"""Reproducible effective-model contexts for essentiality simulations.

The canonical SBML file is never mutated here.  A caller receives an in-memory
model after the declared medium and, optionally, the declared strain overlay
have been applied.  This exact order is shared by screening and acceptance
replay so an FN dossier cannot silently bind a raw model to an overlaid screen.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any
import warnings

import pandas as pd
from cobra.io import read_sbml_model

from .strain_overlay import apply_strain_overlay, load_strain_profile


LEUCINE_EXCHANGE_ID = "R1219"
SIMULATION_CONTEXT_FINGERPRINT_VERSION = "1"
STRAIN_OVERLAY_EFFECT_FINGERPRINT_VERSION = "1"


def canonical_json(value: Any) -> str:
    """Serialize JSON-compatible values deterministically."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_payload(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def load_media(path: Path) -> pd.DataFrame:
    media = pd.read_csv(path, dtype={"exchange": str})
    required = {"exchange", "uptake"}
    missing = required - set(media.columns)
    if missing:
        raise ValueError(f"Media file missing columns {sorted(missing)}: {path}")
    media["exchange"] = media["exchange"].str.strip()
    media["uptake"] = pd.to_numeric(media["uptake"], errors="raise")
    if media["exchange"].duplicated().any():
        repeated = sorted(media.loc[media["exchange"].duplicated(), "exchange"].unique())
        raise ValueError(f"Duplicate exchanges in media file: {repeated}")
    invalid = media[(media["uptake"] < 0) | (media["uptake"] > 1000)]
    if not invalid.empty:
        raise ValueError("Media uptake values must be between 0 and 1000")
    return media


def apply_media(model, media: pd.DataFrame) -> dict[str, float]:
    """Replace all uptake with the explicit medium definition."""
    model_reactions = {reaction.id for reaction in model.reactions}
    missing = sorted(set(media["exchange"]) - model_reactions)
    if missing:
        raise ValueError(f"Media exchange reactions not found in model: {missing}")
    allowed = {
        row.exchange: float(row.uptake)
        for row in media.itertuples(index=False)
        if float(row.uptake) > 0
    }
    model.medium = allowed
    if LEUCINE_EXCHANGE_ID in model.medium:
        raise ValueError(
            f"SD-Leu invariant violated: {LEUCINE_EXCHANGE_ID} is open for uptake"
        )
    return dict(model.medium)


def _excluded_runtime_genes(profile: dict[str, Any] | None) -> tuple[str, ...]:
    if profile is None:
        return ()
    return tuple(
        sorted(
            {
                str(operation["pseudo_gene"])
                for operation in profile.get("operations", [])
                if operation.get("type") == "plasmid_complement"
            }
        )
    )


@dataclass(frozen=True)
class EffectiveSimulationContext:
    """The exact effective model and provenance used by one simulation."""

    model: Any
    active_medium: dict[str, float]
    excluded_runtime_genes: tuple[str, ...]
    model_path: Path
    media_path: Path
    strain_profile_path: Path | None
    strain_profile: dict[str, Any] | None
    strain_overlay_audit: dict[str, Any] | None
    canonical_model_sha256: str
    medium_sha256: str
    strain_profile_sha256: str | None
    strain_overlay_effect_sha256: str | None
    simulation_context_fingerprint: str

    @property
    def strain_overlay_enabled(self) -> bool:
        return self.strain_profile is not None

    @property
    def strain_profile_id(self) -> str | None:
        if self.strain_profile is None:
            return None
        return str(self.strain_profile["profile_id"])

    def provenance(self) -> dict[str, Any]:
        """Path-free binding fields safe to copy into durable evidence."""
        return {
            "simulation_context_fingerprint_version": SIMULATION_CONTEXT_FINGERPRINT_VERSION,
            "simulation_context_fingerprint": self.simulation_context_fingerprint,
            "strain_overlay_enabled": self.strain_overlay_enabled,
            "strain_profile_id": self.strain_profile_id,
            "strain_profile_sha256": self.strain_profile_sha256,
            "strain_overlay_effect_fingerprint_version": (
                STRAIN_OVERLAY_EFFECT_FINGERPRINT_VERSION
                if self.strain_overlay_enabled
                else None
            ),
            "strain_overlay_effect_sha256": self.strain_overlay_effect_sha256,
        }


def load_effective_simulation_context(
    *,
    model_path: str | Path,
    media_path: str | Path,
    strain_profile_path: str | Path | None = None,
) -> EffectiveSimulationContext:
    """Load SBML, then apply medium and optional overlay in that order."""
    model_path = Path(model_path).resolve()
    media_path = Path(media_path).resolve()
    profile_path = Path(strain_profile_path).resolve() if strain_profile_path else None
    if not model_path.is_file() or not media_path.is_file():
        raise ValueError("Simulation context requires existing model and medium files")
    if profile_path is not None and not profile_path.is_file():
        raise ValueError(f"Strain profile does not exist: {profile_path}")

    media = load_media(media_path)
    profile = load_strain_profile(profile_path) if profile_path is not None else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(str(model_path))
    active_medium = apply_media(model, media)
    audit: dict[str, Any] | None = None
    if profile is not None:
        audit = apply_strain_overlay(model, profile, active_medium=active_medium)
        active_medium = dict(model.medium)

    model_sha = sha256_file(model_path)
    medium_sha = sha256_file(media_path)
    profile_sha = sha256_file(profile_path) if profile_path is not None else None
    effect_sha = sha256_payload(audit) if audit is not None else None
    context_payload = {
        "canonical_model_sha256": model_sha,
        "medium_sha256": medium_sha,
        "strain_profile_sha256": profile_sha,
        "strain_overlay_effect_sha256": effect_sha,
    }
    return EffectiveSimulationContext(
        model=model,
        active_medium=active_medium,
        excluded_runtime_genes=_excluded_runtime_genes(profile),
        model_path=model_path,
        media_path=media_path,
        strain_profile_path=profile_path,
        strain_profile=profile,
        strain_overlay_audit=audit,
        canonical_model_sha256=model_sha,
        medium_sha256=medium_sha,
        strain_profile_sha256=profile_sha,
        strain_overlay_effect_sha256=effect_sha,
        simulation_context_fingerprint=sha256_payload(context_payload),
    )


__all__ = [
    "EffectiveSimulationContext",
    "LEUCINE_EXCHANGE_ID",
    "SIMULATION_CONTEXT_FINGERPRINT_VERSION",
    "STRAIN_OVERLAY_EFFECT_FINGERPRINT_VERSION",
    "apply_media",
    "load_effective_simulation_context",
    "load_media",
    "sha256_file",
    "sha256_payload",
]
