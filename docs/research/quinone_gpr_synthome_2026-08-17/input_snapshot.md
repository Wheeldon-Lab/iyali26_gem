# Quinone GPR and synthome review: frozen input snapshot

Date: 2026-08-17

This review is read-only with respect to the canonical model and the durable
essentiality workflow. It may support a disposable runtime counterfactual, but
it does not authorize changes to `model.xml`, canonical GPRs, reaction bounds,
curated patch tables, or the FN dossier.

## Simulation context

- Model: `model.xml`
- Model SHA-256: `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`
- Model size: 2,313 reactions; 1,877 metabolites; 1,074 genes
- Medium: PO1f SD-Leu runtime medium
- Medium SHA-256: `ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85`
- Strain profile: `po1f_sd_leu_accrispr_v1`
- Profile SHA-256: `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383`
- Overlay effect SHA-256: `d15acbde9438f5d2391c4da23705a34a3585833062d616517d5af052088606c2`
- Composite simulation-context fingerprint: `c243b23e7344e3f1e2b4962be25f0f2a38980990c6fe88ee32d3aa4f7af90e30`
- Essential-positive reference SHA-256: `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`
- Cas9/Cas12a fitness table SHA-256: `97ef559651cd99ce63144ffe08ae64e983ecaec60d568f859da8926c14c7d9ff`
- Solver for model-side baseline: Gurobi
- Wild-type growth: `1.4650760568106092 h^-1`

## Current CoQ9 route

The exact route inspected in the frozen model is:

`R763 -> R407 -> R969 (reverse, mitochondrion to cytosol) -> R39 -> R808 -> R715 -> R40 -> R19 -> R18 -> R695 -> R385`

All 11 route reactions pass the current formula/charge balance check. `R39`
and `R19` have no GPR. `YALI1B20527g` (COQ8 candidate) and
`YALI1F34675g` (COQ9 candidate) are retained gene objects but have no reaction
associations.

## Current model-side essentiality observation

With no CoQ9 net demand in the canonical model, all nine reviewed COQ-gene
knockouts have `KO/WT = 1.0` (within solver tolerance) and are nonessential at
the 1%, 5%, 10%, and 15% growth cutoffs. This is a statement about the current
model representation, not evidence that the biological genes are dispensable.

## Human and evidence gates

- Cross-species evidence must be labeled as such.
- Fitness score and model growth ratio are separate scales and must not be
  converted into one another.
- A runtime dependency requires direct perturbation or omission plus a
  diagnostic intermediate, rescue, or equivalent reaction-specific evidence.
- Rate enhancement or general synthome stabilization is not encoded as a
  binary GPR.
- No formal EGC candidate may exceed `awaiting_human`; formal implementation
  requires a later explicit `接受 EGC-...` command.
