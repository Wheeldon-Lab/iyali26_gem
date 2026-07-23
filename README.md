[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://mcnaughtonadm.github.io/iyali26)

# iyali26

A history report will be publicly visible at https://wheeldon-lab.github.io/iyali26_gem.


## Usage

The canonical model remains `model.xml`; the immutable pipeline starting model
remains `data/iyali26.xml`. Research inputs, curation state, external databases,
and generated reports live in a separate research workspace.

Configure it explicitly:

```bash
export IYALI26_RESEARCH_ROOT=/absolute/path/to/iyali26_gem_research
python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT"
```

The CLI option takes precedence over `IYALI26_RESEARCH_ROOT`. Data-dependent
commands fail closed when the workspace or a required gate file is missing.
Local compatibility symlinks may preserve old `data/...`, `results/...`, and
`experiments/...` paths, but those links are never committed.

Run the test suite through either supported entry point:

```bash
.venv/bin/pytest -q
.venv/bin/python -m pytest -q
```

## Repository and research boundaries

The Git repository contains code, tests, engineering documentation,
`data/iyali26.xml`, and `model.xml`. The external workspace layout is documented
in `config/research-workspace.example.toml`:

- `raw/`: source PDF and spreadsheet files;
- `reference/`: MetaNetX, KEGG, NCBI, ExPASy, and external models;
- `state/`: essentiality dossiers, ledgers, media, and curation tables;
- `experiments/`: experiment definitions and inputs;
- `artifacts/`: reports, diagnostics, legacy models, and weekly briefings;
- `cache/`: retained download and tool caches;
- `snapshots/`: checksum-verified pre-relocation snapshots.

`scripts/workspace_relocator.py` implements the zero-delete relocation and
verification workflow. `relocation_manifest.csv` maps every original path to its
archived destination and SHA-256.

## MetaNetX

Place downloaded MetaNetX files under
`$IYALI26_RESEARCH_ROOT/reference/metanetx/`:

- https://www.metanetx.org/ftp/latest/chem_prop.tsv
- https://www.metanetx.org/ftp/latest/chem_xref.tsv
- https://www.metanetx.org/ftp/latest/reac_xref.tsv

---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
