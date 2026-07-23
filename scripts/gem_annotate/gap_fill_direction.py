"""Validated direction curation for reactions inserted by gap filling.

The prioritized gap-fill table records an equation but not whether that
equation is biologically reversible.  This module keeps the evidence-backed
direction decision in a separate durable table so generated candidate tables
can be regenerated without erasing manual curation.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from .config import CURATION_DATA_DIR

DEFAULT_GAP_FILL_DIRECTION_TABLE = (
    CURATION_DATA_DIR / "gap_fill_direction_curation.csv"
)

_ALLOWED_STATUSES = {"active", "needs_review"}
_ALLOWED_ACTIONS = {"keep", "reverse"}
_REQUIRED_COLUMNS = {
    "schema_version",
    "status",
    "bigg_reaction",
    "mnxr_id",
    "stoichiometry_action",
    "lower_bound",
    "upper_bound",
    "evidence_url",
    "rationale",
}


@dataclass(frozen=True)
class GapFillDirectionRow:
    """One evidence-backed direction decision for a gap-fill reaction."""

    schema_version: int
    status: str
    bigg_reaction: str
    mnxr_id: str
    stoichiometry_action: str
    lower_bound: float
    upper_bound: float
    evidence_url: str
    rationale: str


def load_gap_fill_direction_curation(
    table_path: str | Path = DEFAULT_GAP_FILL_DIRECTION_TABLE,
) -> dict[str, GapFillDirectionRow]:
    """Load and validate the durable gap-fill direction table."""

    path = Path(table_path)
    if not path.exists():
        raise FileNotFoundError(f"gap-fill direction curation not found: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or ())
        missing = sorted(_REQUIRED_COLUMNS - columns)
        if missing:
            raise ValueError(
                "gap-fill direction curation is missing columns: "
                + ", ".join(missing)
            )

        rows: dict[str, GapFillDirectionRow] = {}
        errors: list[str] = []
        for line_number, raw in enumerate(reader, start=2):
            try:
                schema_version = int(raw["schema_version"])
                if schema_version != 1:
                    raise ValueError("schema_version must be 1")

                status = raw["status"].strip()
                if status not in _ALLOWED_STATUSES:
                    raise ValueError(
                        f"status must be one of {sorted(_ALLOWED_STATUSES)}"
                    )

                bigg_reaction = raw["bigg_reaction"].strip()
                mnxr_id = raw["mnxr_id"].strip()
                if not bigg_reaction:
                    raise ValueError("bigg_reaction must be nonempty")
                if not mnxr_id.startswith("MNXR"):
                    raise ValueError("mnxr_id must start with MNXR")

                action = raw["stoichiometry_action"].strip()
                if action not in _ALLOWED_ACTIONS:
                    raise ValueError(
                        "stoichiometry_action must be 'keep' or 'reverse'"
                    )

                lower_bound = float(raw["lower_bound"])
                upper_bound = float(raw["upper_bound"])
                if not math.isfinite(lower_bound) or not math.isfinite(upper_bound):
                    raise ValueError("bounds must be finite")
                if lower_bound > upper_bound:
                    raise ValueError("lower_bound must not exceed upper_bound")

                evidence_url = raw["evidence_url"].strip()
                rationale = raw["rationale"].strip()
                if not evidence_url.startswith("https://"):
                    raise ValueError("evidence_url must be an https URL")
                if not rationale:
                    raise ValueError("rationale must be nonempty")
                if bigg_reaction in rows:
                    raise ValueError(
                        f"duplicate bigg_reaction {bigg_reaction!r}"
                    )

                rows[bigg_reaction] = GapFillDirectionRow(
                    schema_version=schema_version,
                    status=status,
                    bigg_reaction=bigg_reaction,
                    mnxr_id=mnxr_id,
                    stoichiometry_action=action,
                    lower_bound=lower_bound,
                    upper_bound=upper_bound,
                    evidence_url=evidence_url,
                    rationale=rationale,
                )
            except (TypeError, ValueError) as exc:
                errors.append(f"line {line_number}: {exc}")

    if errors:
        raise ValueError(
            "invalid gap-fill direction curation:\n- " + "\n- ".join(errors)
        )
    return rows


__all__ = [
    "DEFAULT_GAP_FILL_DIRECTION_TABLE",
    "GapFillDirectionRow",
    "load_gap_fill_direction_curation",
]
