# gem_annotate — refactored package

Pure code-motion refactor of `update_model.py` into a proper Python package.
No logic was changed; only file organisation and import structure.

## Package layout

```
scripts/gem_annotate/
├── __init__.py       re-exports main() for backward compatibility
├── __main__.py       enables `python -m gem_annotate`
├── config.py         path constants, medium dicts, URL constants
├── io.py             MetaNetX TSV loaders
├── metabolites.py    metabolite annotation + H+/H2O balancing
├── exchange.py       exchange-bound calibration
├── biomass.py        biomass reaction R1372 diagnosis and fixing
├── reactions.py      reaction annotation via reac_xref
├── genes.py          gene annotation via UniProt REST + NCBI E-utilities
└── gaps.py           gap analysis (FVA), MIS audit, metabolite merging
```

## How to run

From the repository root:

```bash
# As a module (recommended)
python -m scripts.gem_annotate

# Or directly
python scripts/gem_annotate/main.py

# Or from inside the scripts/ directory
cd scripts
python -m gem_annotate
```

## Module responsibilities

| Module | Contents |
|---|---|
| `config.py` | `REPO_ROOT`, `STARTING_MODEL_PATH`, `OUTPUT_MODEL_PATH`, `MNX_DIR`, `MINIMAL_MEDIUM_BIGG`, `MINIMAL_MEDIUM_NAMES`, UniProt/NCBI URL constants |
| `io.py` | `_read_tsv`, `load_chem_xref`, `load_chem_prop`, `load_reac_xref`, `load_reac_prop` |
| `metabolites.py` | `_FORMULA_RE`, `_parse_name_formula`, `_bigg_id`, `annotate_metabolites`, `_build_compartment_met_index`, `fix_proton_water_balance` |
| `exchange.py` | `set_exchange_bounds` |
| `biomass.py` | `_formula_mw`, `fix_biomass_reaction` |
| `reactions.py` | `annotate_reactions` |
| `genes.py` | `_normalise_locus_tag`, `_merge_gene_annotation`, `_parse_uniprot_entry`, `_fetch_proteome`, `_tier_a`, `_tier_b`, `_tier_ncbi`, `_build_kegg_uniprot_index`, `_enrich_kegg_genes`, `annotate_genes` |
| `gaps.py` | `DUPLICATE_PAIRS`, `find_gaps`, `report_gaps`, `audit_mis_reactions`, `merge_duplicate_metabolites` |
| `main.py` | Pipeline orchestration only; calls into all modules above |
