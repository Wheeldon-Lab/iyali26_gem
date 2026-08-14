# PO1f runtime overlay and uracil-bound audit

Date: 2026-07-29

## Question

The released iYali26 model is a W29/CLIB89 model, whereas the essentiality
screen used PO1f in SD-Leu. This experiment asks whether the PO1f genotype and
the available uracil can be represented at runtime without modifying canonical
`model.xml`.

## Evidence-backed strain layer

Ramesh et al. report PO1f as
`MatA, leu2-270, ura3-302, xpr2-322, axp-2` and cultured the guide-library
transformants in SD-Leu. The library vector was derived from pCRISPRyl, whose
depositor record identifies LEU2 as its selectable marker.

| Gene identity | Protein function and evidence status | Runtime representation |
|---|---|---|
| `YALI1E31685g` (legacy `YALI0E26741g`), **URA3** | Orotidine-5′-phosphate decarboxylase; curated enzyme identity, with `ura3-302` directly reported for PO1f | Disable GPR-less model reaction `R612` |
| `YALI1C00464g` (legacy `YALI0C00407g`), **LEU2** | 3-isopropylmalate dehydrogenase; curated enzyme identity, with `leu2-270` directly reported for PO1f | Replace the genomic `R45` GPR at runtime with `PO1f_plasmid_LEU2` to represent plasmid complementation |
| `YALI0F31889g`, **XPR2** | Alkaline extracellular protease; experimentally established protein function | Genotype provenance only; no current small-molecule GEM reaction |
| `YALI1B07696g` (legacy `YALI0B05654g`), **AXP1** | Acid extracellular protease; curated annotation | Genotype provenance only; no current small-molecule GEM reaction |

Cas9/Cas12a integration at A08 is assay provenance only. KU70 is intact in the
fitness screen; the ΔKU70 strains were used for guide cutting-score calibration
and are not part of this metabolic overlay.

## Uracil bound audit

The SD-Leu recipe uses 0.69 g/L Sunrise CSM-Leu. The manufacturer lists
20 mg/L uracil at that formulation:

\[
C_{\mathrm{uracil}}
=
\frac{20\ \mathrm{mg/L}}{112.09\ \mathrm{mg/mmol}}
=
0.178428\ \mathrm{mmol/L}.
\]

The former `R1354=0.01607` value was obtained from:

\[
10\ \frac{\mathrm{mmol}}{\mathrm{gDW\,h}}
\times
\frac{0.178428\ \mathrm{mmol/L\ uracil}}
{111.015\ \mathrm{mmol/L\ glucose}}
=
0.0160724\ \frac{\mathrm{mmol}}{\mathrm{gDW\,h}}.
\]

The arithmetic is correct, but the calculation assumes uracil and glucose have
the same concentration-normalized clearance. No PO1f uracil \(V_{\max}\),
\(K_m\), specific uptake rate, residual-uracil time course, or biomass integral
was found. Therefore 0.01607 is a concentration-ratio surrogate, not a
validated kinetic upper bound.

For static essentiality FBA, the overlay now sets `R1354=1000` with the explicit
meaning “available and not artificially rate-limiting; not experimentally
measured.” The physical 0.178428 mmol/L remains formulation metadata for a
future batch/dFBA extracellular pool.

## In-memory results

With `R612` disabled:

| R1354 uptake bound | Growth, h⁻¹ |
|---:|---:|
| 0 | 0 |
| 0.01607 | 0.116990 |
| 0.2 | 1.414631 |
| 1 | 1.414631 |
| 1000 | 1.414631 |

At bound 1000, the optimized model uses only
`0.1943166 mmol·gDW⁻¹·h⁻¹` uracil. Across 1074 model genes, bounds 0.2, 1, and
1000 produced identical knockout growth ratios to numerical precision
(maximum absolute difference below \(2.4\times10^{-14}\)).

The PO1f four-cutoff screen produced:

| Growth cutoff | TP | FN | Recall |
|---:|---:|---:|---:|
| 1% | 57 | 265 | 17.70% |
| 5% | 63 | 259 | 19.57% |
| 10% | 67 | 255 | 20.81% |
| 15% | 79 | 243 | 24.53% |

The one additional TP at every cutoff is `YALI1D07232g`, an uncharacterized
multi-pass membrane protein assigned by the model to `R935` uracil transport.
That GPR is model-dependent and is not independently validated by the improved
essentiality match. The current optimal salvage route also uses GPR-less
`R1308` in reverse plus `R1927`; no salvage GPR was changed in this experiment.

## Implementation and safety

- Canonical `model.xml` is not written or changed.
- The external JSON profile and its SHA-256 are included in the run key,
  summary, and manifest.
- The summary records a composite simulation-context fingerprint.
- The runtime plasmid pseudo gene is excluded from genomic deletion screening.
- Formal FN dossier generation is blocked for runtime overlays until the
  durable evidence workflow can replay the strain-profile fingerprint.

## Sources

- [Ramesh et al. 2023, acCRISPR methods](https://www.nature.com/articles/s42003-023-04996-8)
- [Sunrise Science CSM-Leu formulation](https://sunrisescience.com/shop/growth-media/amino-acid-supplement-mixtures/csm-formulations/csm-leu-powder-10-grams/)
- [Addgene pCRISPRyl #70007](https://www.addgene.org/70007/)
