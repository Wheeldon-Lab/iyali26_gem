# Task: Refactor `update_model.py` into a structured Python package

## Context

`update_model.py` is a ~1639-line annotation pipeline for the iYli21 genome-scale metabolic model (GEM) of *Yarrowia lipolytica*. It currently lives as a single flat file. Refactor it into a proper Python package without changing any logic.

## Target package structure

```
gem_annotate/
├── __init__.py          # re-export main() for backward compat
├── main.py              # entry point: orchestrates the pipeline
├── config.py            # all constants and path definitions
├── io.py                # MetaNetX TSV loaders
├── metabolites.py       # metabolite annotation + H+/H2O balancing
├── exchange.py          # exchange bound calibration + medium dicts
├── biomass.py           # biomass reaction diagnosis and fixing
├── reactions.py         # reaction annotation via reac_xref
├── genes.py             # gene annotation via UniProt REST + NCBI
└── gaps.py              # gap analysis, MIS audit, metabolite merging
```

## Module responsibilities

### `config.py`
Move here:
- `REPO_ROOT`, `STARTING_MODEL_PATH`, `OUTPUT_MODEL_PATH`, `MNX_DIR`
- `MINIMAL_MEDIUM_BIGG`, `MINIMAL_MEDIUM_NAMES` dicts
- UniProt/NCBI URL constants: `_UNIPROT_SEARCH_URL`, `_KEGG_CONV_URL`, `_PROTEOME_IDS`, `_NCBI_ESEARCH_URL`, `_NCBI_EFETCH_URL`, `_NCBI_EFETCH_BATCH`, `_TIER_B_LIMIT`

### `io.py`
Move here:
- `_read_tsv()`
- `load_chem_xref()`
- `load_chem_prop()`
- `load_reac_xref()`
- `load_reac_prop()`

### `metabolites.py`
Move here:
- `_FORMULA_RE`, `_parse_name_formula()`
- `_bigg_id()`
- `annotate_metabolites()`
- `_build_compartment_met_index()`
- `fix_proton_water_balance()`

### `exchange.py`
Move here:
- `set_exchange_bounds()`
- (medium dicts stay in `config.py`, imported here)

### `biomass.py`
Move here:
- `_formula_mw()`
- `fix_biomass_reaction()`

### `reactions.py`
Move here:
- `annotate_reactions()`

### `genes.py`
Move here:
- `_normalise_locus_tag()`
- `_merge_gene_annotation()`
- `_parse_uniprot_entry()`
- `_fetch_proteome()`
- `_tier_a()`
- `_tier_b()`
- `_tier_ncbi()`
- `annotate_genes()`

### `gaps.py`
Move here:
- `find_gaps()`
- `report_gaps()`
- `audit_mis_reactions()`
- `merge_duplicate_metabolites()`
- `DUPLICATE_PAIRS` constant (move from `main()` body)

### `main.py`
Keep only the orchestration logic from `main()`. Import everything from the modules above.

## Rules

1. **Zero logic changes.** Move code only. Do not rename functions, change signatures, or alter algorithms.
2. **Fix all imports.** Each module imports only what it needs. Use relative imports within the package (`from .config import ...`).
3. **Shared logger.** Each module creates its own logger: `logger = logging.getLogger(__name__)`. Remove the `logging.basicConfig()` call from module level — move it to `main.py` only.
4. **Backward compatibility.** The package should be runnable as `python -m gem_annotate` or `python gem_annotate/main.py`.
5. **No new dependencies.** Same imports as the original file.
6. **Keep all docstrings and comments intact.**

## Deliverables

- All files in `gem_annotate/` as described above.
- A brief `README_refactor.md` documenting the new structure and how to run it.

## Source file

The original `update_model.py` is provided. Use it as the sole source of truth.
