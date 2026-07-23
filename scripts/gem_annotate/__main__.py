import argparse
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_MODEL_PATH = REPO_ROOT / "model.xml"

parser = argparse.ArgumentParser(description="Build the iYali26 model")
parser.add_argument(
    "--research-root",
    type=Path,
    help=(
        "External research workspace. Overrides IYALI26_RESEARCH_ROOT and "
        "is required unless local compatibility links configure the workspace."
    ),
)
parser.add_argument(
    "--provisional-capacity-profile",
    type=Path,
    help=(
        "Build a separately named experimental model with globally active, "
        "marked isozyme-capacity hypotheses"
    ),
)
parser.add_argument(
    "--output-model",
    type=Path,
    default=OUTPUT_MODEL_PATH,
    help="Output SBML path (default: canonical model.xml)",
)
parser.add_argument(
    "--trna-biomass-mode",
    choices=("split",),
    help=(
        "Build a separately named experimental model with 20 independent "
        "AA-tRNA -> tRNA + protein-residue biomass coupling reactions"
    ),
)
args = parser.parse_args()

if args.research_root is not None:
    os.environ["IYALI26_RESEARCH_ROOT"] = str(args.research_root.resolve())

from .config import load_project_paths

project_paths = load_project_paths(args.research_root, required=True)
project_paths.require(
    project_paths.metanetx,
    project_paths.cache,
    project_paths.essentiality,
    project_paths.media,
    project_paths.curation_data,
)

if args.provisional_capacity_profile is not None:
    if not args.provisional_capacity_profile.exists():
        parser.error(
            f"capacity profile not found: {args.provisional_capacity_profile}"
        )
    if args.output_model.resolve() == OUTPUT_MODEL_PATH.resolve():
        parser.error(
            "a provisional capacity profile cannot overwrite canonical model.xml; "
            "provide a separate --output-model"
        )
if (
    args.trna_biomass_mode is not None
    and args.output_model.resolve() == OUTPUT_MODEL_PATH.resolve()
):
    parser.error(
        "an experimental tRNA biomass overlay cannot overwrite canonical "
        "model.xml; provide a separate --output-model"
    )

from .main import main

main(
    provisional_capacity_path=args.provisional_capacity_profile,
    trna_biomass_mode=args.trna_biomass_mode,
    output_model_path=args.output_model,
)
