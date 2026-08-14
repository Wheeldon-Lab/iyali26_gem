## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: run
- Origin Date: 2026-08-06
- Verification Status: UNVERIFIED (deterministic replay required)
- Version Label: exp_result_v1

# Runtime-only Quinone counterfactual experiment

## Experiment Result

- **ID**: `quinone-runtime-counterfactual-20260806`
- **Type**: simulation
- **Status**: completed
- **Base context fingerprint**: `a85218c0be14e14f1b82f4e9c3aec57322080f9d4e44a80df0f4754f66fa53db`
- **Canonical model SHA-256**: `0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415`
- **Baseline B–PO1f growth**: `1.46507605681 h^-1`

## Main diagnosis

The frozen route is inactive because it has no effective terminal net demand:
all route FVA ranges are `[0,0]` without a demand.  An irreversible diagnostic
demand activates the complete route, including the previously omitted connector
`R969` in reverse.

This diagnosis survives two atom-balanced terminal variants:

- Legacy Q6 zero-growth stoichiometric reachability maximum:
  `0.934496499` mmol/gDW/h,
  but legacy `R385` is chemically imbalanced.
- Atom-balanced oxidized-terminal Q6 zero-growth reachability maximum:
  `0.934280331` mmol/gDW/h
  after making `R385` produce ubiquinone-6.
- Atom-balanced oxidized-terminal Q9 zero-growth reachability maximum:
  `0.673485721` mmol/gDW/h
  after the explicit four-IPP nonaprenyl counterfactual.
- Closing `R385` reduces all three diagnostic demands to zero, so the result is
  not supplied by an alternate modeled source.

Therefore **missing net demand is the immediate cause of inactivity in this
simulation context**. It is not the only model defect: the legacy endpoint is
imbalanced, the biologically appropriate terminal redox microspecies is
unresolved, the chain length is CoQ6 rather than the evidence-supported CoQ9
counterfactual, and the native localization/GPR dependency remains unresolved.

## Growth-dilution sensitivity

The following coefficients are numerical sensitivity points, not fitted
physiological parameters:

| cQ (mmol/gDW) | Growth (h⁻¹) | Q9 dilution flux |
|---:|---:|---:|
| 0 | 1.46507606 | 0 |
| 1e-06 | 1.46507334 | 1.46507334e-06 |
| 1e-05 | 1.46504891 | 1.46504891e-05 |
| 0.0001 | 1.46480464 | 0.000146480464 |
| 0.001 | 1.46236646 | 0.00146236646 |
| 0.01 | 1.43842368 | 0.0143842368 |

## Essentiality attribution

The full positive-only B–PO1f screen was run for the frozen baseline and for a
Q9 dilution coefficient of `0.001 mmol/gDW` under two
GPR interpretations. Recall is TP / `322` experimental
positive genes that map into the screened model. The source positive list has
`1612` genes in total; genes absent from the model are not
counted as TP or FN.

| Scenario | Cutoff | TP | FN | Recall |
|---|---:|---:|---:|---:|
| Q0_B_PO1f | 1% | 57 | 265 | 17.70% |
| Q0_B_PO1f | 5% | 63 | 259 | 19.57% |
| Q0_B_PO1f | 10% | 67 | 255 | 20.81% |
| Q0_B_PO1f | 15% | 79 | 243 | 24.53% |
| Q9_current_repeated_AND | 1% | 66 | 256 | 20.50% |
| Q9_current_repeated_AND | 5% | 72 | 250 | 22.36% |
| Q9_current_repeated_AND | 10% | 77 | 245 | 23.91% |
| Q9_current_repeated_AND | 15% | 88 | 234 | 27.33% |
| Q9_step_specific_GPR | 1% | 64 | 258 | 19.88% |
| Q9_step_specific_GPR | 5% | 70 | 252 | 21.74% |
| Q9_step_specific_GPR | 10% | 75 | 247 | 23.29% |
| Q9_step_specific_GPR | 15% | 86 | 236 | 26.71% |

- Positive-list call changes caused by adding balanced Q9 chemistry and demand:
  `37` gene-threshold rows.
- Positive-list call changes caused by route-wide step-specific GPR remapping,
  including adding candidate GPRs to previously ungated `R39/R40` and replacing
  the repeated seven-gene `AND` assignments:
  `8` gene-threshold rows.

These counterfactual call changes are causal model diagnostics, not evidence
that the proposed GPRs or coefficient are biologically correct.

## Chemistry boundary

- Q9 is built by reinterpreting the existing chain-specific IDs in memory and
  adding `C15H24` to each Q6 intermediate.
- Counterfactual `R763` is
  `4 IPP + pentaprenyl-PP -> 4 PPi + nonaprenyl-PP`.
- The oxidized-terminal counterfactual `R385` is
  `SAM + 3-demethylubiquinone-9 -> SAH + ubiquinone-9`.
- Every reaction in the Q9 route passes elemental and charge balance.
- [KEGG R08781](https://www.kegg.jp/entry/R08781) and
  [Rhea 81218](https://www.rhea-db.org/rhea/81218) represent generic oxidized
  quinone equations, whereas [Rhea 44381](https://www.rhea-db.org/rhea/44381)
  and Q9-specific [Rhea 17049](https://www.rhea-db.org/rhea/17049) represent a
  reduced ubiquinol equation with proton production. These are
  database-compatible alternatives; they do not resolve the native Yarrowia
  terminal redox state or supply a PO1f biomass coefficient.

## Anomalies and limitations

1. The current canonical `R385` has residual `H:+1, charge:-1`.
2. The native terminal substrate/product redox microspecies remains unresolved.
3. No W29/PO1f SD-Leu whole-cell Q9 mmol/gDW value is available.
4. `R39/R969/R808` localization and transport remain unverified in native
   *Yarrowia*.
5. The step-specific GPR is the minimal direct-catalyst interpretation; native
   COQ8/COQ9 accessory necessity remains unresolved.
6. The `0...1000` demand is used only as a reachability objective and is not a
   biological demand coefficient.

## No-write verification

The experiment wrote reports only.  It did not write SBML, curated tables,
formal GPRs/bounds, the durable FN ledger, or Obsidian.  See the frozen evidence
review in `docs/research/quinone_review_2026-08-06/review_report.md`.
