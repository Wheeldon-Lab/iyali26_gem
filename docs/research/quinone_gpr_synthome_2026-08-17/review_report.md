# iYali26 Quinone GPR and Synthome Evidence Review — 2026-08-17

## Executive decision

This review reached the prespecified stop conditions. It does not authorize a
new Quinone GPR, an accessory `AND`, a CoQ9 abundance or degradation
coefficient, a runtime reference mapping, or an EGC candidate.

The decisive findings are:

1. `YALI1A08781g — no established gene name (COQ6 candidate) — FAD-dependent
   CoQ monooxygenase` (**curated/computational annotation; native Yarrowia
   function unverified**) is compatible with the comparative C5 regiochemistry
   represented by `R39`, but the exact carboxylated-substrate order remains
   unresolved and the current reaction is cytosolic while direct
   yeast experiments place Coq6 on the matrix side of the mitochondrial inner
   membrane. The modeled `R969/R808` hydrophobic-intermediate shuttle and the
   Yarrowia electron-transfer partners are unverified. `R39` is therefore
   chemistry-supported but topology/electron-bookkeeping blocked.
2. The exact `R19` product is not the direct product made by COQ6 in the
   reconstructed system: COQ6 produces a hydroquinone, while `R19` produces an
   oxidized quinone. Native yeast C1 modification can also persist without
   Coq6. The current reaction lumps hydroxylation and oxidation and cannot
   receive a binary COQ6 GPR.
3. `YALI1B20527g — no established gene name (COQ8 candidate) — ADCK-family
   ATPase/kinase-like CoQ-synthome accessory` (**curated annotation; native
   Yarrowia role unverified**) has pathway-level regulation and enhancement
   evidence, but no reviewed reaction-specific absolute dependency.
4. `YALI1F34675g — no established gene name (COQ9 candidate) — lipid-binding
   and substrate-presenting CoQ-synthome accessory` (**curated annotation;
   native Yarrowia role unverified**) has strong comparative accessory evidence,
   but no default binary reaction dependency. The strongest candidate,
   COQ9–COQ7/`R695`, is directly conflicted across systems.
5. Whole-cell total CoQ9+CoQ9H2 in W29/CLIB89 or PO1f under SD-Leu or a matched
   growth condition was not located in mmol/gDCW. A Yarrowia CoQ9 half-life or
   molecular degradation rate was also not located. No `c_Q` or `k_deg` value
   is authorized.

`NOT LOCATED` means that the bounded review did not find qualifying evidence;
it does not mean that the parameter or biological process is zero.

## Frozen context

- Canonical model SHA-256:
  `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`
- PO1f SD-Leu simulation-context fingerprint:
  `c243b23e7344e3f1e2b4962be25f0f2a38980990c6fe88ee32d3aa4f7af90e30`
- Solver: Gurobi
- Wild-type growth: `1.4650760568106092 h^-1`
- Current route:
  `R763 -> R407 -> R969(reverse) -> R39 -> R808 -> R715 -> R40 -> R19 -> R18 -> R695 -> R385`

All 11 route reactions pass the current formula/charge balance check. `R39`
and `R19` have empty GPRs. COQ8 and COQ9 candidate gene objects have no reaction
associations. The canonical model has no CoQ9 net demand.

## Evidence gate

A new runtime binary dependency required direct perturbation or omission plus a
diagnostic intermediate, rescue, or equivalent reaction-specific evidence.
General complex membership, protein stabilization, or rate enhancement was not
treated as absolute reaction incapacity. Cross-species evidence was allowed,
but its species, construct, substrate and condition limits were retained.

### Gate outcome

| Candidate mapping | Chemistry | Topology / exact product | Absolute dependency | Decision |
|---|---|---|---|---|
| COQ6 -> `R39` | Comparative C5 regiochemistry supported; exact carboxylated substrate sequence unresolved | Current cytosolic location and shuttle unsupported; electron partners absent | Native Yarrowia dependency not located | Stop; separate chemistry/topology case |
| COQ6 -> `R19` | Reconstructed animal COQ6 can hydroxylate an R19-type substrate | Direct product is hydroquinone, not current oxidized quinone; native yeast C1 dependency is contradicted | Not established | Stop; separate chemistry/redox case |
| COQ8 -> any reviewed late step | Pathway accessory activity supported | No atom-changing step assigned | No reaction-specific absolute dependency | No mapping |
| COQ9 -> `R39` or `R19` | COQ6 cooperation/stability supported | Effect is indirect and system-dependent | No absolute dependency | No mapping |
| COQ9 -> `R695` | COQ7 cooperation and substrate presentation supported | Cell and purified systems disagree on necessity | Conflicted | No default mapping |
| COQ9 -> `R715/R40/R18/R385` | Synthome association only | No step-specific perturbation | Unsupported | No mapping |

## COQ6 chemistry and topology

### R39

The C5 regiochemistry represented in `R39` is compatible with comparative
COQ6 evidence, although the exact carboxylated-substrate order in native yeast
remains unresolved. *S. cerevisiae* genetics, diagnostic intermediates and
precursor rescue identify Coq6 as the C5 enzyme
([Ozeir et al. 2011](https://doi.org/10.1016/j.chembiol.2011.07.008),
Fig. 2, 3 and 5), while direct substrate conversion is demonstrated in a
reconstructed short-chain system. Direct fractionation and
protease-protection experiments place yeast Coq6 on the matrix side of the
mitochondrial inner membrane ([Gin et al. 2003](https://doi.org/10.1074/jbc.M303234200),
Fig. 7b–e).

The current model instead places `R39` in the cytosol and uses two GPR-less
transport reactions, `R969` and `R808`, to move highly hydrophobic prenylated
intermediates out of and back into mitochondria. No biological support for
this Yarrowia shuttle was located.

Coq6 also does not behave as an enzyme reduced by NAD(P)H alone. Bound FAD and
the need for a ferredoxin/reductase chain are supported by protein studies
([Ismail et al. 2016](https://doi.org/10.1371/journal.pcbi.1004690), Fig. 1 and
Supplementary Fig. S2; [Gonzalez et al. 2024](https://doi.org/10.1002/cbic.202300738),
Fig. 7A–D). The reconstructed COQ system likewise required FDXR/FDX2, FAD and
NADPH for the C5 conversion ([Nicoll et al. 2024](https://doi.org/10.1038/s41929-023-01087-z),
Fig. 2). The current balanced `0.5 O2` equation is a formal net representation,
not the demonstrated enzyme mechanism.

### R19

The reconstructed ancestral tetrapod system shows that COQ6 can hydroxylate an
R19-type substrate, but its product is a hydroquinone. The model product is an
oxidized quinone, so `R19` combines COQ6 hydroxylation with an additional
two-electron oxidation ([Nicoll et al. 2024](https://doi.org/10.1038/s41929-023-01087-z),
Fig. 4). The oxidation may be enzymatic, oxygen-driven or spontaneous; it does
not automatically inherit a COQ6 GPR.

There is also lineage counterevidence. In native *S. cerevisiae*, C1
decarboxylation and hydroxylation can proceed in Coq6-deficient backgrounds,
and the C1 hydroxylase was unresolved ([Ozeir et al. 2015](https://doi.org/10.1074/jbc.M115.675744),
Fig. 1–4 and Discussion). Bacterial pathways use different specialist or
multifunctional FMOs depending on lineage, so family membership alone cannot
fix the Yarrowia C1 assignment.

## COQ8 and COQ9 accessory dependence

COQ8 acute inhibition can nearly stop new CoQ formation and cause early
intermediate accumulation, but it does not resolve a single atom-changing late
step ([Reidenbach et al. 2018](https://doi.org/10.1016/j.chembiol.2017.11.001),
Fig. 6–7). In the reconstructed pathway, COQ8B plus ATP enhanced overall output
by roughly five-fold while the individual COQ3, COQ6 and COQ7:COQ9 reactions
were not strictly dependent on it ([Nicoll et al. 2024](https://doi.org/10.1038/s41929-023-01087-z),
Fig. 6 and Extended Data Fig. 8). This supports a regulatory/organizational
accessory, not a direct GPR for every step.

COQ9 binds hydrophobic CoQ intermediates and interacts with COQ7 and other
synthome proteins. Co-immunoprecipitation, lipid-binding, structures and
interface mutants support substrate presentation and complex stability, not an
independent atom-changing reaction. The most important conflict concerns the
COQ7 step:

- In a *S. cerevisiae* `delta-coq9` background with COQ8 overexpression,
  precursor bypass can leave diagnostic DMQ accumulation consistent with a
  strong cellular Coq7-step block ([Xie et al. 2012](https://doi.org/10.1074/jbc.M112.360354),
  Fig. 8).
- Purified reconstructed COQ7 remains active without COQ9, and COQ9 increases
  activity by only about 1.5-fold ([Nicoll et al. 2024](https://doi.org/10.1038/s41929-023-01087-z),
  Fig. 5–6).

The difference can reflect substrate delivery and complex stability in intact
cells rather than intrinsic catalytic necessity. Without Yarrowia knockout
intermediates and rescue, it cannot be encoded as a default binary `AND`.

No native Yarrowia COQ8/COQ9 AP-MS, co-IP, proximity-labeling,
reaction-specific knockout metabolomics, diagnostic-intermediate or rescue
study was located. This is a direct-species evidence gap, not evidence that the
accessories are biologically unimportant.

## CoQ9 pool and turnover

The literature supports Q9 identity and a nonzero physical pool in Yarrowia,
but none of the located numbers meets the target coefficient definition:

| Evidence | Reported quantity | Why it is not `c_Q` or `k_deg` |
|---|---:|---|
| Historical *Endomycopsis lipolytica* | Q9 system identified | No absolute amount or dry-weight denominator |
| ATCC 20362 patent | 0.2–0.3% of extracted oil = 0.002515–0.003773 mmol/g oil | Oil denominator, unmatched culture, no recovery or uncertainty |
| Purified PIPO complex I | 0.2 and 0.4 Q9/complex | Protein-complex occupancy, not whole-cell pool |
| Purified GB30 delta-ST1 complex I | 1.9 Q9/complex | Preparation-specific occupancy, not whole-cell pool |
| *S. cerevisiae* whole-cell measurements | Dry-weight-normalized Q6/Q6H2 values | Different species, homolog, medium and growth phase |
| Plant turnover experiment | approximately 30 h UQ9/UQ10 half-life | Different kingdom and physiology |

The target values remain:

- `c_Q` for W29/CLIB89 or PO1f, total CoQ9+CoQ9H2 in mmol/gDCW under SD-Leu
  or a defensibly matched growth condition: **NOT LOCATED**.
- Yarrowia CoQ9 molecular turnover, half-life or degradation rate separated
  from biomass dilution: **NOT LOCATED**.

The runtime H-Q9-1 parameters therefore remain hypothetical sensitivity
parameters. They are not estimates.

## Fitness score and model cutoffs

For the CRISPR experiments,

\[
2^{FS}=\frac{\text{normalized guide abundance after selection}}
{\text{normalized guide abundance in the control}}.
\]

This is not a growth-rate ratio. The paper's essential calls use assay-specific
q values. In the normalized dataset, the empirical boundary is approximately
`-2.918` for Cas9 and `-1.855` for Cas12a, but those values describe these two
assays only.

The model uses a different quantity:

\[
r_g=\frac{\mu_{KO}}{\mu_{WT}},
\]

with strict 1%, 5%, 10% and 15% growth cutoffs. In the current canonical
PO1f/B-group/Gurobi context, all nine reviewed genes have \(r_g\approx1\)
because the CoQ9 pool has no net growth demand and three candidates are not
connected to reactions. The complete values are in
`fitness_threshold_matrix.tsv`.

The Cas9/Cas12a pattern is heterogeneous: two genes are essential in both
assays (COQ7 and COQ8 candidates), two are nonessential in both (COQ3 and COQ5
candidates), and five are discordant. These phenotypes prioritize unresolved
biology; they do not identify a reaction, compartment, pool coefficient or
binary accessory rule.

## Runtime decision

The evidence-approved mapping set is empty. Therefore the planned new
nine-gene mapped counterfactual was **not run by prespecified stop condition**.
Running an empty mapping would duplicate the canonical baseline, while running
the conflicted COQ9–R695 hypothesis would violate the rule that only supported
absolute dependencies enter the reference runtime experiment.

The previously archived H-Q9-1 study remains unchanged and may be cited only
as `runtime_only` and `sensitivity_only_not_calibrated`. Its result that six
currently connected step-specific genes show conditional complete blocks while
COQ6/COQ8/COQ9 remain unchanged follows from the current GPR representation;
it is not independent validation of that representation.

No new runtime runner or test was needed. A general accessory toolkit was not
created: native COBRA GPR `AND` would already be the minimal representation if
a future absolute dependency passes the evidence gate.

## Formal workflow decision

- No EGC candidate was created.
- No `awaiting_human` record was created.
- `model.xml`, formal GPRs, reaction bounds, biomass, curated data and the FN
  dossier remain unchanged.
- A future formal patch still requires an independently audited chemistry and
  identity review plus an explicit user command naming `接受 EGC-...`.

## Reopen conditions

Reopen the COQ6 case only after evidence resolves the mitochondrial reaction
location, `R969/R808` topology, Yarrowia ferredoxin/reductase identities, and
the separate hydroquinone-to-quinone oxidation near `R19`.

Reopen the accessory case only after Yarrowia reaction-specific perturbation
plus a diagnostic intermediate, rescue, or equivalent omission experiment.

Reopen parameterization only after a target-strain measurement reports total
CoQ9+CoQ9H2 with matched DCW, growth phase, internal standard/recovery,
biological replication and uncertainty, or a target-strain time course
separates dilution from molecular degradation.

## Audit artifacts

- `input_snapshot.md`
- `source_inventory.tsv`
- `evidence_ledger.tsv`
- `conflict_matrix.tsv`
- `gene_evidence_matrix.tsv`
- `fitness_threshold_matrix.tsv`
- `source_audit.tsv`
- `research_manifest.json`

The independent source audit reopened all 20 prespecified primary/grey sources
and audited all 19 material atomic claims. Its verdict distribution was eight
supported, four partially supported, five contradicted and two unverified.
Kogan 1985 was available only as a bibliographic record/abstract and was not
misrepresented as a full-text audit. No direct Yarrowia reaction-level
COQ6/COQ8/COQ9 experiment was located. Full corrective wording and locators
are recorded in `source_audit.tsv`.

The permanent human gate remains active: fitness/recall, automated database
annotation and cross-species synthome membership cannot create a GPR,
accessory `AND`, CoQ9 coefficient or EGC. Any future formal patch still needs
current-fingerprint chemistry, identity and topology reviews, an independent
source audit, and an explicit user command naming `接受 EGC-...`.
