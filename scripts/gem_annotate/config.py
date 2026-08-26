"""
config.py — all constants and path definitions for the gem_annotate pipeline.
"""

import os
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RESEARCH_ROOT_ENV = "IYALI26_RESEARCH_ROOT"


@dataclass(frozen=True, slots=True)
class ProjectPaths:
    """Resolved code-repository and external research-workspace paths."""

    repo_root: Path
    research_root: Path
    configured: bool

    @property
    def starting_model(self) -> Path:
        return self.repo_root / "data" / "iyali26.xml"

    @property
    def output_model(self) -> Path:
        return self.repo_root / "model.xml"

    @property
    def metanetx(self) -> Path:
        return self.research_root / "reference" / "metanetx"

    @property
    def cache(self) -> Path:
        return self.research_root / "cache" / "data"

    @property
    def essentiality(self) -> Path:
        return self.research_root / "state" / "essentiality" / "repository"

    @property
    def media(self) -> Path:
        return self.research_root / "state" / "media"

    @property
    def strain_profiles(self) -> Path:
        return self.research_root / "state" / "strain_profiles"

    @property
    def curation_data(self) -> Path:
        return self.research_root / "state" / "curation" / "data"

    @property
    def locus_map(self) -> Path:
        return self.research_root / "reference" / "locus_map"

    @property
    def experiments(self) -> Path:
        return self.research_root / "experiments" / "repository"

    @property
    def results(self) -> Path:
        return self.research_root / "artifacts" / "results"

    @property
    def weekly_briefing(self) -> Path:
        return self.research_root / "artifacts" / "weekly_briefing"

    def resolve_legacy_path(self, path: str | Path) -> Path:
        """Resolve a former repository-relative research path externally."""
        candidate = Path(path)
        if candidate.is_absolute():
            return candidate.resolve()
        parts = candidate.parts
        if parts[:2] == ("data", "essentiality"):
            return (self.essentiality.joinpath(*parts[2:])).resolve()
        if parts[:2] == ("data", "media"):
            return (self.media.joinpath(*parts[2:])).resolve()
        if parts[:2] == ("data", "strain_profiles"):
            return (self.strain_profiles.joinpath(*parts[2:])).resolve()
        if parts[:2] == ("data", "metanetx"):
            return (self.metanetx.joinpath(*parts[2:])).resolve()
        if parts[:2] == ("data", "cache"):
            return (self.cache.joinpath(*parts[2:])).resolve()
        if parts[:2] == ("data", "ncbi"):
            return (self.research_root / "reference" / "ncbi").joinpath(*parts[2:]).resolve()
        if parts[:2] == ("data", "kegg"):
            return (self.research_root / "reference" / "kegg").joinpath(*parts[2:]).resolve()
        if parts[:2] == ("data", "yali1_yali0_map"):
            return (self.locus_map.joinpath(*parts[2:])).resolve()
        if parts == ("data", "iyali26.xml"):
            return (self.repo_root / candidate).resolve()
        if parts and parts[0] == "data" and len(parts) == 2:
            return (self.curation_data / parts[1]).resolve()
        if parts and parts[0] == "experiments":
            return (self.experiments.joinpath(*parts[1:])).resolve()
        if parts and parts[0] == "results":
            return (self.results.joinpath(*parts[1:])).resolve()
        if parts[:2] == ("docs", "weekly_briefing"):
            return (self.weekly_briefing.joinpath(*parts[2:])).resolve()
        return (self.repo_root / candidate).resolve()

    def require(self, *paths: Path) -> "ProjectPaths":
        if not self.configured:
            raise RuntimeError(
                "The external research workspace is not configured. Set "
                f"{RESEARCH_ROOT_ENV} or pass --research-root."
            )
        if not self.research_root.is_dir():
            raise FileNotFoundError(
                f"Research workspace not found: {self.research_root}"
            )
        missing = [path for path in paths if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Required research workspace paths are missing: "
                + ", ".join(str(path) for path in missing)
            )
        return self


def _infer_research_root_from_compatibility_links(repo_root: Path) -> Path | None:
    """Infer the root from compatibility links created by workspace relocation."""
    for candidate in (
        repo_root / "data" / "essentiality",
        repo_root / "data" / "metanetx",
        repo_root / "results",
    ):
        if not candidate.is_symlink():
            continue
        resolved = candidate.resolve()
        for parent in (resolved, *resolved.parents):
            if parent.name in {"state", "reference", "artifacts"}:
                return parent.parent
    return None


def load_project_paths(
    research_root: str | Path | None = None,
    *,
    required: bool = False,
) -> ProjectPaths:
    """Resolve an explicit root, then the environment, then local symlinks."""
    configured = research_root
    if configured is None:
        configured = os.environ.get(RESEARCH_ROOT_ENV)
    if configured is None:
        configured = _infer_research_root_from_compatibility_links(REPO_ROOT)
    if configured is None:
        paths = ProjectPaths(
            repo_root=REPO_ROOT,
            research_root=REPO_ROOT / ".research-unconfigured",
            configured=False,
        )
    else:
        paths = ProjectPaths(
            repo_root=REPO_ROOT,
            research_root=Path(configured).expanduser().resolve(),
            configured=True,
        )
    if required:
        paths.require()
    return paths


PROJECT_PATHS = load_project_paths()
STARTING_MODEL_PATH = PROJECT_PATHS.starting_model
OUTPUT_MODEL_PATH = PROJECT_PATHS.output_model
MNX_DIR = PROJECT_PATHS.metanetx
CACHE_DIR = PROJECT_PATHS.cache
ESSENTIALITY_DIR = PROJECT_PATHS.essentiality
MEDIA_DIR = PROJECT_PATHS.media
CURATION_DATA_DIR = PROJECT_PATHS.curation_data
LOCUS_MAP_DIR = PROJECT_PATHS.locus_map
RESULTS_DIR = PROJECT_PATHS.results

# ── UniProt / NCBI / KEGG URL constants ──────────────────────────────────────
_UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
_KEGG_CONV_URL      = "https://rest.kegg.jp/conv/uniprot/yli"
_TIER_B_LIMIT       = 1100   # per-gene UniProt search cap (~3 min at 0.15 s/call)
# W29/CLIB89 = UP000182444 (Other proteome); CLIB122 reference = UP000001300
_PROTEOME_IDS       = ("UP000182444", "UP000001300")
# NCBI E-utilities — Tier A-prime (bulk NCBI Gene, covers all YALI1 locus tags)
_NCBI_ESEARCH_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EFETCH_URL    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_NCBI_EFETCH_BATCH  = 400    # IDs per efetch request (NCBI recommends ≤500)

# ── Minimal aerobic medium for Y. lipolytica W29 on glucose ──────────────────
# Primary key: BiGG metabolite ID (no compartment suffix).
# Fallback key: lowercased chemical name (before trailing "_FORMULA" in iYali26 names).
# Value: uptake lower bound (mmol / gDW / h).  -1000 = effectively unlimited.
MINIMAL_MEDIUM_BIGG: dict[str, float] = {
    # carbon
    "glc__D":   -10.0,    # D-glucose, sole carbon source
    # nitrogen
    "nh4":      -1000.0,  # ammonium
    # oxygen / aerobic
    "o2":       -1000.0,  # oxygen
    # phosphorus
    "pi":       -1000.0,  # inorganic phosphate
    "h2po4":    -1000.0,  # dihydrogen phosphate (alt BiGG form)
    # sulfur
    "so4":      -1000.0,  # sulphate / sulfate
    # proton and water (exchange reactions are always in extracellular compartment)
    "h":        -1000.0,  # proton
    "h2o":      -1000.0,  # water
    # trace metals
    "fe2":      -1000.0,  # Fe2+
    "fe3":      -1000.0,  # Fe3+
    # macroions
    "k":        -1000.0,  # potassium
    "na1":      -1000.0,  # sodium
    # CO2 (can be re-assimilated)
    "co2":      -1000.0,
}

# Fallback name-based keys for metabolites that lack bigg.metabolite annotation.
# Lowercased chemical name extracted from iYali26 "name_FORMULA" convention.
MINIMAL_MEDIUM_NAMES: dict[str, float] = {
    "d-glucose":      -10.0,
    "glucose":        -10.0,
    "ammonium":       -1000.0,
    "oxygen":         -1000.0,
    "phosphate":      -1000.0,
    "sulphate":       -1000.0,
    "sulfate":        -1000.0,
    "h+":             -1000.0,
    "h2o":            -1000.0,
    "water":          -1000.0,
    "iron":           -1000.0,
    "potassium":      -1000.0,
    "sodium":         -1000.0,
    "carbon dioxide": -1000.0,
    "co2":            -1000.0,
}
