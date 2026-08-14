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
- `state/`: essentiality dossiers, ledgers, media, runtime strain profiles, and
  curation tables;
- `experiments/`: experiment definitions and inputs;
- `artifacts/`: reports, diagnostics, legacy models, and weekly briefings;
- `cache/`: retained download and tool caches;
- `snapshots/`: checksum-verified pre-relocation snapshots.

`scripts/workspace_relocator.py` implements the zero-delete relocation and
verification workflow. `relocation_manifest.csv` maps every original path to its
archived destination and SHA-256.

## Model patches and research runs

`scripts.gem_annotate` is the only canonical model build implementation.
`scripts/update_model.py` remains solely as a compatibility forwarding entry
point; its prior monolithic source is retained in the external archive.

For an individual, curated patch, use the shared runner and always choose a
new output path. It refuses both the canonical `model.xml` and an existing
output file:

```bash
python -m scripts.gem_annotate.patch_runner \
  --patch c161-pool-extension \
  --input-model data/iyali26.xml \
  --output-model /tmp/iyali26-c161-experiment.xml
```

The four historical patch scripts remain compatibility wrappers around this
same runner. Their archived source remains available under the external
research workspace's `archive/legacy_code/` directory.

Research runs are recorded append-only in
`$IYALI26_RESEARCH_ROOT/artifacts/run_registry.jsonl`. A matching successful
fingerprint is rejected by default; use `--force-rerun --reproduction-reason
"..."` only for an intentional reproduction. Existing manifests are immutable;
to register historical runs and retained duplicate artifacts, run:

```bash
python -m scripts.gem_annotate.run_registry backfill \
  --research-root "$IYALI26_RESEARCH_ROOT"
```

### B-group aminoacyl-tRNA biomass experiment

The fully split B-group translation representation is now part of the
released canonical `model.xml`. The build pipeline applies 20 independent
`AA-tRNA -> tRNA + protein-residue` reactions, applies the current SD-Leu
medium and the SHA-pinned PO1f runtime profile, then runs the positive-only
1%, 5%, 10%, and 15% essentiality screen:

```bash
python -m scripts.gem_annotate.trna_biomass_pipeline \
  --research-root "$IYALI26_RESEARCH_ROOT"
```

The screening command writes a separately named copy and a provenance manifest
under the external research workspace. It verifies that the copied model
preserves the released canonical B-group reaction/metabolite structure.
The current PO1f regression is 57/265, 63/259, 67/255, and 79/243 TP/FN at
1%, 5%, 10%, and 15%, respectively.

The PO1f profile is stored at
`$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json`. It is applied
only in memory: `R612` is disabled for the `ura3-302` background, while `R45`
is associated with a runtime `PO1f_plasmid_LEU2` pseudo gene to represent
LEU2 complementation by the guide-library plasmid. The pseudo gene is not
treated as a screened genomic target.

The CSM-Leu formulation supplies 20 mg/L uracil (0.178428 mmol/L), but this
concentration does not determine a flux bound in mmol/gDW/h. The profile
therefore retains the formulation and legacy concentration-ratio calculation
as provenance, then applies `R1354=1000` as a permissive static-FBA availability
bound explicitly marked as not experimentally measured. Batch simulations
must instead use 0.178428 mmol/L as an initial extracellular pool.

## MetaNetX

Place downloaded MetaNetX files under
`$IYALI26_RESEARCH_ROOT/reference/metanetx/`:

- https://www.metanetx.org/ftp/latest/chem_prop.tsv
- https://www.metanetx.org/ftp/latest/chem_xref.tsv
- https://www.metanetx.org/ftp/latest/reac_xref.tsv

---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
