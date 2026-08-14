# Quinone multi-agent literature review manifest

- Review date: `2026-08-06`
- Status: `complete — read-only evidence map`
- Governance verdict: `GO for read-only use; no authorization for model, code, curated-data, or Obsidian changes`
- Repository branch: `codex/workspace-cleanup`
- Frozen commit: `35c959b3032b14661653a5bdd8eb2f10c11d5495`

## Frozen inputs

| Input | SHA-256 |
|---|---|
| `model.xml` | `0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415` |
| `data/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| `docs/curation/quinone_branch_cleanup.md` | `b1e030c5548b7abc3e9b9afcc6c2b204fd97993a04ad0f0a9275364b47ad61b4` |
| Obsidian `17-Quinone Biosynthesis.md` | `071ae5ad3a444b0e2cbc2b7117e0083eda0c23a2452b68077c2e348a65ec69a7` |

All four hashes were rechecked after the review and were unchanged.

## Audit summary

- Material claims independently audited: `45/45 (100%)`
- Auditor verdicts: `22 supported`, `7 partially_supported`, `5 unsupported`, `6 contradicted`, `5 unverified`
- High-impact sources independently opened: `23/26 (88.5%)`
- Peer-reviewed high-impact sources independently opened: `21/24 (87.5%)`
- Unopened high-impact records: `SRC-QD-001`, `SRC-QD-003`, `SRC-QD-004`; their limitations are explicit in the ledger and report.
- Final adversarial verdict: `GO for read-only evidence map`.
- Final artifact QA: all TSV schemas, identifiers, foreign keys, counts, locators, source-opened semantics, and overstatement boundaries passed after renaming mixed conflict references to `record_ids_a` / `record_ids_b`.

## Independent work strands

1. Native CoQ homolog, chain length, chemistry, and localization evidence.
2. Net CoQ demand, pool dilution, biomass coupling, and quantitative evidence.
3. COQ synthome identity crosswalk and reaction-specific GPR evidence.
4. Independent source audit performed without treating agent agreement as evidence.
5. Independent adversarial review performed after evidence synthesis.
6. Independent artifact QA performed after the scientific review.

## Artifact hashes

| Artifact | SHA-256 |
|---|---|
| `adversarial_review.md` | `a1a68b95dfd6584a55c15906d1e03fe788277253a9f9324c759d22705a77a1fc` |
| `audit_protocol.md` | `d6943c81178760eeff3216ac8f2d38214c6e934f7e729ff176f055988c88d216` |
| `conflict_matrix.tsv` | `cbcd4e9ca341bc133707fb914d9418b29db4a805edb6015e64e847efc7e0ee1b` |
| `evidence_ledger.tsv` | `9f09a2000639f228fcecbe56375f28624bc035235e2a0675e6a9cd3de0d3f07e` |
| `input_snapshot.md` | `652307385f6db0340acfd55dc015c03224e086316444cf8ad7599fee25c10de2` |
| `local_model_provenance.md` | `e56289f148f8791648fa85cbd98f7a070c8d0a5bdc145dba13afe25cb6bf0414` |
| `quantitative_data.tsv` | `a5e9d6e169cba07eff0f67174b2f1ad5b8cc520dd91aabeee09e956adba71fa1` |
| `review_report.md` | `0018710f3fbaa6672ac79a2b3749b8c589a86de32369625282b4c89fa7a9da60` |
| `search_log.md` | `427805fe7e413369ff857cede6d149c5f9975843af6f06fb3c50e636b5fa1746` |
| `source_audit.tsv` | `6042cc0f00ec3da70200e01e118c43354aa16b21b9dce052a1db6e065ccc4fc9` |
| `source_inventory.tsv` | `6576b5943726ff6811814fc648e17e615cae354260fb95c8996798bed73d534f` |

The manifest is not self-hashed. Any later artifact edit invalidates the corresponding value above and requires a manifest refresh.

## Human gate

This review does not authorize changing CoQ chain length, metabolite formulas,
reaction compartments, GPRs, biomass/dilution demand, sinks, or any other model
content. A separate, explicit human decision is required before designing or
executing a counterfactual or curated patch.
