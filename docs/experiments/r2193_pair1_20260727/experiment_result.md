# R2193 Pair 1 iron-uptake GPR experiment

Date: 2026-07-27  
Status: completed as an isolated evidence-labelled experiment  
Canonical-model change: none  
Canonical `model.xml` SHA-256: `b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48`

## Question

Can the inherited R2193 GPR `YALI1F07747g` be replaced by the Pair 1
iron-uptake candidate on the basis of the 2026 *Yarrowia lipolytica*
experiment?

## Evidence boundary

The 2026 study cloned a `YlFtr1` ORF whose primers map exactly to:

`YALI1D08564g` / `YALI0D06688g` / NCBI Gene `2910500` / UniProt `Q6CA15`.

The experiment co-expressed VHb and YlFtr1 and observed increased
intracellular Fe2+. It did not include an Ftr1-only control,
knockout/complementation, uptake kinetics, or membrane-localization assay.
It therefore provides same-species experimental support for the FtrA
candidate, but not a strict single-gene causal test.

The proposed FetC partner:

`YALI1D08684g` / `YALI0D06754g` / NCBI Gene `2910503` / UniProt `Q6CA12`

was not manipulated in that study. Its multicopper-oxidase identity and its
AND relationship with `YALI1D08564g` remain inferred.

Primary source:
<https://doi.org/10.1016/j.synbio.2026.02.004>

## Implemented experimental variants

Both variants use the latest B-group AA-tRNA biomass model as input
(SHA-256 `b028f8d0cb4cd3b56e274980bd0656aba382376e8dc7d498ad0a7d81669c7d29`).
Neither variant changes R2193 stoichiometry or bounds.

| Variant | R2193 GPR | Evidence status | Model SHA-256 |
|---|---|---|---|
| Primary evidence-bounded model | `YALI1D08564g` | same-species experimental support | `9b8b108610b7193c5706081914db3423908ed984a9373a5e88ce8f16763ac07d` |
| Full Pair 1 sensitivity model | `YALI1D08564g and YALI1D08684g` | inferred Pair 1 hypothesis | `609867529253af0d51c03bed636b543e35fb6e9e511c07088b326937824eb23b` |

The inherited `YALI1F07747g` no longer controls R2193 in either experimental
model. COBRA retains its gene record as an unconnected orphan; the audit
records this explicitly.

## Screen-test result

Both variants produced the same positive-only SD-Leu screen:

| Growth cutoff | TP | FN | Recall |
|---:|---:|---:|---:|
| 1% | 56 | 266 | 17.39% |
| 5% | 62 | 260 | 19.25% |
| 10% | 66 | 256 | 20.50% |
| 15% | 78 | 244 | 24.22% |

WT growth was `1.3993195284 h^-1` in both variants.

At 100% optimal growth:

- R1189 FVA: `0 to 0`
- R2193 FVA: `0 to 0`
- R863 FVA: `0 to 0`
- `YALI1D08564g` knockout growth ratio: `1.0`
- `YALI1D08684g` knockout growth ratio in the Pair 1 AND model: `1.0`

Thus, replacing the GPR corrects the proposed protein identity but does not
activate iron flux. The porphyrin/heme/siroheme branches remain inactive
because the current B-group biomass and downstream reactions do not require
those cofactors.

## Canonical gate

The canonical model was not changed because:

1. direct evidence covers `YALI1D08564g`, not the complete Pair 1 AND rule;
2. R2193 currently uses generic iron metabolites with missing formulas, so
   the schema-v2 balanced-chemistry gate cannot pass;
3. no current evidence dossier, skeptic pass, or exact human-accepted
   `EGC-...` case exists for this GPR replacement.

The isolated experiment allows the candidate to be tested without presenting
the inferred FetC relationship as an accepted canonical annotation.
