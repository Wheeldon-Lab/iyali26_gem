# Cofactor-coupling proteomics prior experiment

Date: 2026-07-25  
Status: completed as an evidence-only experiment  
Canonical-model change: none  
Canonical `model.xml` SHA-256: `b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48`

## Question

This experiment addresses two related questions:

1. Can the current YALI1 genes for sulfite reductase, ferrochelatase, heme synthesis, and respiratory-chain reactions be mapped to YALI0, NCBI Gene, and primary UniProt identifiers and connected to quantitative *Yarrowia lipolytica* proteomics?
2. Can flow matching infer the siroheme occupancy parameter \(\theta\), or a biomass/cofactor-coupling coefficient, from the presently available observations?

## Executive result

The identifier and abundance experiment succeeded for the nuclear genes:

- MET5-like sulfite-reductase beta subunit, YALI1D14058g / Q6C9H5, was quantified in 18/18 samples.
- MET10-like sulfite-reductase alpha subunit, YALI1E19603g / Q6C5P4, was quantified in 18/18 samples.
- Ferrochelatase, YALI1F25924g / Q6C140, was quantified in 18/18 samples.
- Heme O synthase COX10, YALI1F31135g / Q6C0L2, was quantified in 14/18 samples.
- R304/complex IV has quantitative values for only 7 of its 12 current GPR members.
- R305/complex III has quantitative values for 10 of its 11 current GPR members.

The missing respiratory-chain values are not zero-abundance measurements. They are proteins absent from the processed proteomics table, including mitochondrial core subunits. Consequently, this dataset cannot by itself estimate the abundance of an intact complex III or complex IV.

The preliminary siroheme coefficient can be bounded from MET5 and MET10 abundance only after adopting the unverified \(\alpha_2\beta_2\) stoichiometric prior. In all 18 samples, MET5 is the limiting subunit. At full occupancy, \(\theta=1\):

| Physiological phase | Mean coefficient | Sample range |
|---|---:|---:|
| Exponential | \(3.95\times10^{-7}\) mmol siroheme/gDW | \(3.23\times10^{-7}\) to \(4.64\times10^{-7}\) |
| Nitrogen-limited | \(6.90\times10^{-8}\) mmol siroheme/gDW | \(4.19\times10^{-8}\) to \(1.41\times10^{-7}\) |

These values are an order-of-magnitude prior, not a validated iYali26 biomass coefficient. The proteomics experiment used aerobic glycerol mineral medium and engineered strains, not the SD-Leu essentiality condition.

Flow matching cannot presently identify \(\theta\). The available data measure MET5 and MET10 abundance but do not measure siroheme occupancy, siroheme pool size, sulfite-reductase flux, or a condition-matched sulfur phenotype. With no \(\theta\)-dependent observation, the posterior equals the prior. Running FMPE now would give a visually sophisticated but scientifically uninformative result.

## Identifier mapping

The mapping chain was:

\[
\text{current YALI1}
\longleftrightarrow
\text{curated YALI0 crosswalk}
\longleftrightarrow
\text{eciYali5 usage\_prot entry}
\longleftrightarrow
\text{primary UniProt protein}
\]

NCBI Gene identifiers were read from the curated S2 metabolic-gene mapping. Primary UniProt accessions were assigned from the single-gene `usage_prot_*` reactions in eciYali5-GEM, rather than by matching only the terminal digits of YALI0 and YALI1 identifiers.

The complete 29-row mapping and abundance table is in `target_protein_abundance.csv`.

Important mapping qualifications:

- The current model still uses legacy `YALI0F04114g` in R304; the curated crosswalk maps it to YALI1F06244g.
- YALI1M00277g and YALI1M00325g correspond to legacy mitochondrial loci YALIfMp14 and YALIfMp18. The eciYali5 reaction contains NP_075434 and NP_075438, but this experiment did not establish which NP accession belongs to which of the two YALI1 loci. The pair is therefore recorded as unresolved rather than guessed.
- YALI1M00390g / YALIfMp23 corresponds to the mitochondrial complex-III core entry NP_075443 in eciYali5, but it is absent from the processed proteomics table.
- The three R383 gene-to-protein identifier mappings are exact, but their collective biological assignment to the modeled heme O monooxygenase reaction was not independently verified here. Their abundance must not be treated as proof that the R383 GPR is correct.

## Quantitative conversion

The eciYali5 processing code first calculates a normalized total-protein-approach abundance and then applies:

\[
A_i =
\mathrm{TPA}_i
\times 10^3
\times P_{\mathrm{tot}}
\times 0.99
\]

where \(P_{\mathrm{tot}}\) is measured protein content in g protein/gDW. Therefore the values in `data/proteomics.tsv` are interpreted as mg protein/gDW.

For each protein:

\[
E_i\;[\mathrm{mmol/gDW}]
=
\frac{A_i\;[\mathrm{mg/gDW}]}
     {MW_i\;[\mathrm{mg/mmol}]}
\]

The numerical value of molecular weight in Da (g/mol) is also the numerical value in mg/mmol.

Under the temporary \(\alpha_2\beta_2\) prior:

\[
E_{\mathrm{complex}}
=
\min\left(
\frac{E_{\mathrm{MET5}}}{2},
\frac{E_{\mathrm{MET10}}}{2}
\right)
\]

If each complex contains two siroheme molecules and \(\theta\) is the occupied fraction:

\[
c_{\mathrm{siroheme}}
=
2\theta E_{\mathrm{complex}}
=
\theta\min(E_{\mathrm{MET5}},E_{\mathrm{MET10}})
\]

This \(c_{\mathrm{siroheme}}\) is a pool-per-biomass coefficient in mmol/gDW. If it is installed as a fixed biomass coefficient, the dilution flux is:

\[
v_{\mathrm{siroheme\ demand}}
=
\mu\,c_{\mathrm{siroheme}}
\]

where \(\mu\) is the biomass flux in gDW/gDW/h. A fixed coefficient preserves a linear FBA problem. A nonlinear expression appears only if both enzyme abundance and growth are free decision variables and their product is inserted directly as a constraint.

The full sample-level calculation is in `siroheme_prior_by_sample.csv`. Its internal check confirms:

- 18 samples were processed;
- MET5 is limiting in all 18;
- \(c_{\mathrm{siroheme}}=2E_{\mathrm{complex}}\) at \(\theta=1\) for every sample.

## Abundance results

Selected protein-level summaries:

| Target | Exponential mean (mg/gDW) | N-limited mean (mg/gDW) | Detection |
|---|---:|---:|---:|
| MET5 / Q6C9H5 | 0.0620 | 0.0108 | 18/18 |
| MET10 / Q6C5P4 | 0.1396 | 0.0131 | 18/18 |
| Ferrochelatase / Q6C140 | 0.0298 | 0.0114 | 18/18 |
| COX10 / Q6C0L2 | 0.00337 | 0.000536 among detected samples | 14/18 |

The N-limited COX10 mean uses only the five non-missing N-limited measurements. Missing values were not converted to zero.

The large phase dependence is biologically and computationally important. For example, the mean full-occupancy siroheme coefficient inferred from exponential samples is approximately 5.7-fold higher than the N-limited mean. A single universal coefficient would hide this condition dependence.

## Flow-matching identifiability test

### What the existing implementation can do

The present `iyali26_flow` code is a gated two-parameter experiment for:

- `r4_capacity_fraction`
- `r1846_capacity_fraction`

Its simulator, configuration schema, observables, and FMPE network are specific to R4/R1846. It cannot be reused unchanged for siroheme occupancy.

### Why \(\theta\) is not identifiable from the current data

Let the current data be:

\[
D_{\mathrm{prot}} =
\{E_{\mathrm{MET5}},E_{\mathrm{MET10}}\}
\]

These observations determine the maximum complex amount under the assumed stoichiometry, but they do not observe \(\theta\). Therefore:

\[
p(D_{\mathrm{prot}}\mid\theta)
=
p(D_{\mathrm{prot}})
\]

and Bayes' rule gives:

\[
p(\theta\mid D_{\mathrm{prot}})
\propto
p(D_{\mathrm{prot}}\mid\theta)p(\theta)
=
p(D_{\mathrm{prot}})p(\theta)
\propto p(\theta)
\]

Thus the posterior is the prior. A flow-matching model can sample that distribution, but it cannot manufacture information that is absent from the observations.

The five-point sensitivity table in `theta_sensitivity.csv` shows only the deterministic scaling:

\[
c_{\mathrm{siroheme}}(\theta)
=
\theta c_{\mathrm{siroheme}}(1)
\]

It is not a fitted posterior.

### Requirements for a meaningful future flow-matching experiment

A new cofactor-coupling simulator would need:

1. parameters such as \(\theta_{\mathrm{siroheme}}\), possible abundance scaling, and only evidence-supported stoichiometric alternatives;
2. an in-memory model variant that connects siroheme synthesis to R745 or biomass dilution without changing canonical `model.xml`;
3. condition-matched observations that respond to those parameters, preferably absolute siroheme/heme measurements, sulfate or sulfite uptake, sulfide/H2S production, R745 flux, and growth under multiple sulfur conditions;
4. a grid or Sobol sensitivity gate before FMPE training;
5. a stop decision if growth and measured exchanges remain flat across the parameter range.

Only after the sensitivity gate shows informative, non-collinear observables should FMPE be trained. If only WT growth or a single essential/nonessential label is supplied, occupancy, enzyme abundance, \(k_{\mathrm{cat}}\), and biomass dilution remain confounded.

## Decision

This experiment supports the following limited conclusion:

- The Yarrowia proteomics data provide a useful numerical prior for MET5/MET10 abundance and an order-of-magnitude range for a possible siroheme dilution coefficient.
- They do not justify adding the coefficient to canonical biomass yet.
- They do not identify the occupancy parameter \(\theta\).
- They do not support an intact complex-III/IV abundance estimate because required mitochondrial/core subunits are missing.
- They do not validate the current R383 GPR.

No biomass demand, GPR, reaction bound, or `model.xml` content was changed.

## Material Passport

### Input materials

| Material | Role | Locator/version | SHA-256 or stable identifier | Transformation |
|---|---|---|---|---|
| iYali26 `model.xml` | Current YALI1 genes and R302/R304/R305/R383/R384/R745 GPR context | local canonical model | `b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48` | Read only |
| `iyli21_genes_vs_S2.csv` | YALI1-to-YALI0 crosswalk | local curated table | `7c68b3cb244f07848b7ad99d7ed0f95768c9812ddd80311e1cdffa06fff4bae9` | Exact one-to-one mappings only |
| `s2_metabolic_genes.csv` | NCBI Gene identifiers and secondary crosswalk | local curated table | `f7f25bb3e4f3c720d4b0bd534309c5c5f3cfa091bbd22383858d71f8600c1c2e` | Selected target rows |
| `curated_locus_identity_exclusions.csv` | Prevent unsafe cross-assembly aliases | external research workspace via repository compatibility link | `883f1aac417838e923f12a00285545a093b985e90fd87ecfd9e6f29adbeee408` | Applied by safe locus resolver |
| eciYali5-GEM repository | Primary protein IDs, mass, processed abundance, processing method | Git commit `3b97d38d92327a099f3b0f6f9cd92e8aeab78c10` | repository commit | Read only |
| eciYali5 `data/proteomics.tsv` | 18-sample quantitative protein table | same commit | `5d6ab79a3c0c8a46f7286de447f2cb28765b00c8d97046ae1834f69de6511d74` | Target rows selected; missing kept missing |
| eciYali5 `data/uniprot.tsv` | Primary accession and molecular mass | same commit | `68f53654c8403cbd4951d406d5ae93c48fe96e1221e47ee022396a2f5e9e6131` | Target rows selected |
| eciYali5 `model/eciYali5-GEM.txt` | YALI0-to-`usage_prot` mapping | same commit | `13dcc2633060ebb4d87869a79daabfad960cf14f7b1ae7143b580f3e86fcf642` | Single-gene usage reactions parsed |
| `calcTPA.R` | TPA calculation provenance | same commit | `0d0ad3241d751659cf94c85135b73d6c6d40ab4e85e304ff35fffdd63b08a8d4` | Read only |
| `generateProteomics4GECKO.R` | Scaling and units provenance | same commit | `cd8608eae313e574f34ced8c89cad712e86b9368ff41b1da73085f8a2054482a` | Formula inspected |
| `proteinContents.tsv` | Phase-specific total protein content | same commit | `2381de3e02056526a04f1d4df7c45d8b888a559cfff239d85d803a38d1fffc11` | Read only |
| ProteomeXchange PXD072100 | Experimental-condition and acquisition provenance | dataset PXD072100 | PXD072100 | Metadata only |

The eciYali5-GEM repository declares CC-BY-4.0. This experiment stores only target-level derived rows, not a redistributed copy of the full upstream dataset.

### Output materials

| Output | Role | SHA-256 |
|---|---|---|
| `target_protein_abundance.csv` | 29-target mapping and abundance summary | `0840a385119cbcfdd4ec4f144e1ca241055fbf79279f5bc4fb10d1eaca800ed8` |
| `siroheme_prior_by_sample.csv` | 18-sample molar conversion and \(\alpha_2\beta_2\) prior | `1f0b26c0b13b25ccc6dd5b89d1e6ea8cdbbd9780911964de6bb99bf5c7cd90bc` |
| `theta_sensitivity.csv` | Five-point deterministic occupancy sensitivity | `adba62a4a0251825ccc2e0a9ef93d636414281d97ac6b4fe83c3a720136f323e` |

## Limitations

1. PXD072100 used PAR, CAR, and OBE strains in aerobic glycerol mineral medium, with exponential and nitrogen-limited phases. It is not a direct SD-Leu experiment.
2. TPA-derived abundance is model-assisted quantitative proteomics, not a direct purified-standard assay for each target protein.
3. The \(\alpha_2\beta_2\) architecture and two-siroheme-per-complex assumption are temporary cross-species priors, not direct *Y. lipolytica* stoichiometric measurements.
4. Cofactor occupancy \(\theta\) was not measured.
5. Missing proteomics entries mean not quantified, not absent from the cell.
6. The mitochondrial YALI1-to-NP mapping is incomplete for the two complex-IV genes.
7. R383 protein abundance does not validate R383 enzyme identity or GPR logic.
