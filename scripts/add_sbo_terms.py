"""
add_sbo_terms.py — Assign SBO terms to all model objects for Memote compliance.

Rules (only applied when sbo is not already set):
  Metabolites  → SBO:0000247  simple chemical
  Genes        → SBO:0000243  gene
  Reactions:
    exchange        → SBO:0000627
    demand          → SBO:0000628
    sink            → SBO:0000632
    biomass         → SBO:0000629  (id matches BIOMASS/biomass/newBiom/R1372)
    ATP maintenance → SBO:0000630  (id matches MAINTENANCE/ATPM)
    transport       → SBO:0000185  (metabolites span ≥2 compartments)
    otherwise       → SBO:0000176  biochemical reaction

Usage:
  python scripts/add_sbo_terms.py [--model model.xml] [--out model.xml]
"""

import argparse
import logging
import re
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_BIOMASS_RE     = re.compile(r"BIOMASS|biomass|newBiom|R1372")
_MAINTENANCE_RE = re.compile(r"MAINTENANCE|ATPM")


def _set_sbo(obj, term: str) -> bool:
    """Set sbo on obj if not already present. Returns True if set."""
    ann = obj.annotation if isinstance(obj.annotation, dict) else {}
    if "sbo" in ann:
        return False
    merged = dict(ann)
    merged["sbo"] = term
    obj.annotation = merged
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Add SBO terms to GEM objects")
    parser.add_argument("--model", default="model.xml",
                        help="Input SBML model (default: model.xml)")
    parser.add_argument("--out",   default=None,
                        help="Output path (default: overwrite --model)")
    args = parser.parse_args()

    model_path = Path(args.model)
    out_path   = Path(args.out) if args.out else model_path

    if not model_path.exists():
        logger.error(f"Model not found: {model_path}")
        sys.exit(1)

    try:
        from cobra.io import read_sbml_model, write_sbml_model
    except ImportError:
        logger.error("COBRApy not installed: pip install cobra")
        sys.exit(1)

    logger.info(f"Loading model: {model_path}")
    model = read_sbml_model(str(model_path))
    logger.info(f"  {len(model.reactions)} reactions, "
                f"{len(model.metabolites)} metabolites, "
                f"{len(model.genes)} genes")

    # ── Metabolites ───────────────────────────────────────────────────────────
    met_set = met_skip = 0
    for met in model.metabolites:
        if _set_sbo(met, "SBO:0000247"):
            met_set += 1
        else:
            met_skip += 1
    logger.info(f"Metabolites : set={met_set}  already_had_sbo={met_skip}")

    # ── Genes ─────────────────────────────────────────────────────────────────
    gene_set = gene_skip = 0
    for gene in model.genes:
        if _set_sbo(gene, "SBO:0000243"):
            gene_set += 1
        else:
            gene_skip += 1
    logger.info(f"Genes       : set={gene_set}  already_had_sbo={gene_skip}")

    # ── Reactions — build priority sets ──────────────────────────────────────
    exchanges = set(model.exchanges)
    demands   = set(model.demands)
    sinks     = set(model.sinks)

    counts = {
        "exchange":    0,
        "demand":      0,
        "sink":        0,
        "biomass":     0,
        "maintenance": 0,
        "transport":   0,
        "biochemical": 0,
        "skip":        0,
    }

    for rxn in model.reactions:
        # Determine SBO term by priority
        if rxn in exchanges:
            term = "SBO:0000627"
            kind = "exchange"
        elif rxn in demands:
            term = "SBO:0000628"
            kind = "demand"
        elif rxn in sinks:
            term = "SBO:0000632"
            kind = "sink"
        elif _BIOMASS_RE.search(rxn.id):
            term = "SBO:0000629"
            kind = "biomass"
        elif _MAINTENANCE_RE.search(rxn.id):
            term = "SBO:0000630"
            kind = "maintenance"
        elif len({m.compartment for m in rxn.metabolites}) >= 2:
            term = "SBO:0000185"
            kind = "transport"
        else:
            term = "SBO:0000176"
            kind = "biochemical"

        if _set_sbo(rxn, term):
            counts[kind] += 1
        else:
            counts["skip"] += 1

    logger.info(
        f"Reactions   : "
        f"exchange={counts['exchange']}  demand={counts['demand']}  "
        f"sink={counts['sink']}  biomass={counts['biomass']}  "
        f"maintenance={counts['maintenance']}  transport={counts['transport']}  "
        f"biochemical={counts['biochemical']}  already_had_sbo={counts['skip']}"
    )

    logger.info(f"Writing model to: {out_path}")
    write_sbml_model(model, str(out_path))
    logger.info("Done.")


if __name__ == "__main__":
    main()
