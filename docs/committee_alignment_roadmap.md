# Committee-aligned research roadmap

## Material Passport

- **Origin skill:** `academic-research-suite / experiment-agent`
- **Origin mode:** plan
- **Origin date:** 2026-07-20
- **Verification status:** BASELINED — repository evidence checked; reactor and economic hypotheses not yet experimentally verified
- **Version label:** `committee_roadmap_v1`
- **Project baseline:** Git `694e6ac`; `model.xml` SHA-256 `dee22ccb2febe39282ac4fae240b5d56c95a23ebaa2b262702c7b7de0c3c012b`
- **Primary inputs:** committee feedback supplied by the researcher; current GEM repository; existing 13-paper comparison matrix at `/Users/david/.codex/.chatgpt-projects/g-p-6a5c155487c08191b6713fd3bb505787/.research/literature_matrix.md`
- **Known missing material:** calibrated fermentation time series, feedstock specifications and prices, downstream recovery data, and reactor geometry/operating limits

## Research objective

Build and defend a multiscale explanation of *Yarrowia lipolytica* lipid fermentation that connects:

> gene and enzyme complex → reaction capacity and control point → intracellular flux → online reactor signal → control action → yield, titer, productivity, and process viability

The target is not merely a larger or better annotated GEM. The final research contribution should show when the model is biologically credible, how reactor measurements constrain its dynamic state, which control action improves an explicit process objective, and whether that improvement survives whole-process constraints.

## Baseline decision

The current repository is strong enough for curation and hypothesis generation, but not yet reliable enough to make high-confidence real-time control decisions.

| Committee expectation | Current asset | Evidence gap | Immediate implication |
|---|---|---|---|
| Explain enzymes, subunits, control points, and reaction rates | 2,300 reactions, 1,073 genes, explicit GPRs, and essentiality diagnostics | No `kcat`, proteomics, enzyme-capacity, or rate-parameter layer was found in code or documentation | Build a mechanism-and-capacity map before claiming rate control |
| Connect metabolic genetics to lipid production | Key lipid reactions and several complex GPRs are present | Essentiality recall at the primary cutoff is 0.186; 646 reactions lack GPRs; bypass and isozyme redundancy remain substantial | Rank claims by evidence rather than treating every GPR as equally certain |
| Use pH, RQ, CTR, OTR, and feed for real-time control | No reactor/off-gas/control implementation or data schema was found | No time-aligned signals, gas-balance equations, actuator limits, controller, or validation run | Establish the measurement dictionary and dynamic balances first |
| Take a whole-process view | Medium and gene-essentiality inputs are tracked | No feedstock, utilities, downstream, waste, or economic system boundary | Add a process inventory and minimum economic model |
| Support quantitative conclusions | Model QA and curation logs are extensive | 837 reactions are universally blocked (36.38%); mass/charge issues, stoichiometric cycles, biomass inconsistency, and no identified NGAM remain | Control simulations must carry a model-validity warning and sensitivity analysis |

The numerical baseline comes from `data/memote_summary.csv`, `results/essentiality/essentiality_summary.json`, and the matching model fingerprint. These limitations are research tasks, not merely software-cleanup items: each one constrains which biological or process conclusion can be defended.

## Core research questions and hypotheses

1. **Mechanism:** Which enzyme complexes and cofactor/compartment constraints control the transition from carbon uptake to lipid accumulation?
   - Working hypothesis: acetyl-CoA supply, malonyl-CoA formation, fatty-acid synthase capacity, NADPH supply, TAG assembly, and competing β-oxidation form a small set of regime-dependent control points.
   - Falsifier: measured or reconstructed flux changes cannot be explained by the proposed nodes, or perturbing them does not change lipid rate under the relevant reactor condition.

2. **State estimation:** Can online gas, pH, and feed signals identify metabolic phase changes and provide quantitative bounds on substrate uptake and respiratory state?
   - Working hypothesis: reconciled O2/CO2 balances and feed history constrain dynamic fluxes more tightly than static medium bounds alone.
   - Falsifier: inferred rates are not identifiable within useful uncertainty, or phase calls are not reproducible across runs.

3. **Control:** Does a model-informed feed/oxygen strategy improve lipid productivity without violating oxygen-transfer, substrate, or actuator constraints?
   - Working hypothesis: a constrained controller using estimated respiratory state outperforms a fixed feed recipe under feed and oxygen disturbances.
   - Falsifier: improvement disappears under realistic sensor lag, noise, parameter uncertainty, or scale-dependent oxygen-transfer limits.

4. **Whole process:** Is the improvement commercially meaningful after feedstock, aeration/agitation, downstream recovery, and waste handling are included?
   - Working hypothesis: feedstock cost, oxygen-transfer power, titer/productivity, and recovery yield dominate economic sensitivity.
   - Falsifier: the proposed biological/control gain has negligible effect on the minimum selling price or moves burden to an unmodeled unit operation.

## Workstream A — biochemical and metabolic-genetic mechanism

### A1. Required evidence table

Create `docs/mechanistic_control_map.md`. One row represents one mechanistic claim, not just one reaction.

| Field | Required content |
|---|---|
| Biological function | Role in carbon assimilation, precursor supply, redox, lipid synthesis, degradation, or regulation |
| Gene(s) and protein | Current and legacy locus IDs, protein name, evidence source |
| Enzyme/subunit logic | Homomer, heteromer, isozyme, accessory enzyme, or uncertain; expected `AND`/`OR` GPR |
| Compartment | Cytosol, mitochondrion, peroxisome, ER, lipid particle, or transport step |
| Model mapping | Reaction ID, equation, direction, bounds, current GPR, and blocked status |
| Kinetic/capacity evidence | `kcat`, `Km`, enzyme abundance, saturation assumption, temperature/pH context, and uncertainty |
| Control mechanism | Transcriptional, allosteric, substrate/cofactor availability, thermodynamic, or enzyme-capacity control |
| Observable consequence | Expected change in OUR/CTR/RQ, substrate uptake, growth, citrate, lipid rate, or product profile |
| Perturbation evidence | Knockout, overexpression, isotope flux, proteomics, metabolomics, or fermentation phenotype |
| Confidence | E0–E5 evidence level and unresolved contradictions |

Evidence levels:

- **E0:** present only as a model reaction.
- **E1:** gene/protein function is annotated.
- **E2:** subunit/isozyme/compartment logic is supported.
- **E3:** biochemical kinetics or enzyme abundance is available in a relevant context.
- **E4:** in vivo perturbation or flux evidence supports the mechanism.
- **E5:** the mechanism predicts a reactor-scale rate or signal across independent runs.

### A2. First-priority audit nodes

| Node | Current model anchor | Question that must be resolved |
|---|---|---|
| ATP-citrate lyase | `R1894`; GPR currently `YALI1E41315g or YALI1D32268g` | Are these isozymes or required subunits? Does the GPR encode the verified complex architecture? |
| Acetyl-CoA carboxylase/biotin activation | `R88`, `R175`; `R87` is a mitochondrial reaction | Which protein supplies catalytic ACC activity, which activates biotin, and how should their dependence be represented? |
| Fatty-acid synthase | `R1392` and `R1393`; FAS subunits `YALI1B19844g and YALI1B25427g` | Are overall chain-length reactions, added gene requirements, stoichiometry, and reversibility biologically justified? |
| Cytosolic NADPH supply | PPP reactions plus mitochondrial NADP malic enzyme `R538` | Which source supplies lipid NADPH in vivo under the target condition, and is compartmental transfer represented? |
| Glycerol backbone | `R348`, `R350–R353` | Which steps limit glycerol-3-phosphate supply and acylation under nitrogen limitation? |
| TAG assembly | `R1771–R1772` and `R1728–R1729` | Are localization and direction correct? Why is `R1772` disabled while `R1771` is reversible? |
| β-oxidation | `R97–R102`, `R1487–R1498` in the peroxisome | Which POX isozymes act on each chain length, and when does degradation compete with accumulation? |
| Maintenance/respiration | no identified NGAM; respiratory network contains complexes | What ATP and redox demand is needed before respiratory signals can constrain flux? |

### A3. Acceptance gate

Pass **G1** when every first-priority node has an evidence-ranked gene–enzyme–reaction–rate entry, all complex GPRs have an explicit subunit justification, and disputed mappings are labeled rather than silently fixed. Pass **G2** when the corresponding pathway is mass/charge checked, feasible under the target medium, and tested for alternate optima and bypasses.

## Workstream B — reactor measurements, estimation, and control

### B1. Signal dictionary and units

Create `docs/reactor_signal_dictionary.md` before fitting any controller.

| Quantity | Minimum definition | Control use | Main caveat |
|---|---|---|---|
| pH and acid/base addition | calibrated pH; acid/base flow and concentration | constraint control; indirect acid/base production signal | base addition is not a metabolic rate unless buffering and chemistry are modeled |
| DO | liquid dissolved-O2 concentration or % air saturation | oxygen constraint; inner-loop controlled variable | depends on probe dynamics, pressure, temperature, and saturation calibration |
| OTR | `kLa (C*O2 - CL)` on a volumetric molar basis | oxygen supply/capacity limit | is not automatically equal to OUR during transients or gas/liquid accumulation |
| CTR/CER | CO2 transfer/evolution from inlet–outlet gas balance | respiratory and carbon-rate estimate | correct for dry/wet gas, pressure, flow, and gas holdup |
| RQ | `CER / OUR`, or approximately `CTR / OTR` only under stated quasi-steady assumptions and common molar units | metabolic-state indicator | sign convention, dynamics, and non-respiratory CO2 chemistry must be explicit |
| Feed | mass or volume flow, substrate concentration, density, and cumulative feed | principal manipulated input and substrate balance | pump calibration, evaporation, sampling, and changing reactor volume matter |
| Agitation/air/O2 enrichment | rpm, gas flow, inlet O2 fraction, pressure | DO-cascade actuators | actuator saturation and scale-dependent `kLa` must be recorded |

Minimum dynamic balances for a fed-batch baseline are:

```text
dV/dt       = F - Fsample - Fevap
d(XV)/dt    = μ X V
d(SV)/dt    = F Sf - qS X V - S Fsample
d(CL V)/dt  = OTR V - OUR V - CL Fsample
OTR         = kLa (C*O2 - CL)
RQ          = CER / OUR
```

An exponential-feed expression may initialize a feed-forward policy,

```text
F(t) = μset X0 V0 exp(μset t) / (YX/S Sf)
```

but it is not a complete controller. Maintenance, product formation, changing yield, volume, feed density, pump bounds, and oxygen capacity must be included or bounded before use.

### B2. Minimum data contract

For at least one complete fermentation run, retain raw and calibrated values with timestamps:

- reactor volume, sampling and evaporation estimates;
- substrate feed flow, concentration, density, and cumulative mass;
- DCW/biomass, residual substrates, citrate/byproducts, lipid concentration and composition;
- pH, acid/base flow and concentration, temperature, DO;
- agitation, inlet gas flow, inlet O2 fraction, pressure;
- inlet and outlet O2/CO2, off-gas flow or the correction used to infer it;
- sensor calibration, sample frequency, filter settings, missing intervals, and known time delays.

No derived column should overwrite a raw signal. Every calculated rate must store its equation, units, smoothing window, and uncertainty.

### B3. Controller progression

1. **Replay and reconcile:** align timestamps, close volume/gas/carbon balances, quantify residuals, and mark physiological phases.
2. **Estimate:** infer OUR, CER/CTR, RQ, `qS`, `qO2`, growth rate, and their uncertainty; test identifiability.
3. **Baseline:** document the existing DO cascade and feed recipe, including actuator priority and saturation.
4. **Feedback:** test a bounded PID/feed-feedback policy with filtering, delay, anti-windup, and safety interlocks.
5. **Model-based control:** compare feed-forward, feedback, and a constrained MPC/state-estimation design in simulation before reactor use.
6. **Robustness:** perturb yield, `kLa`, sensor lag/noise, feed concentration, and initial biomass; report failures, not only mean improvement.

Primary process metrics are lipid productivity, carbon yield, titer, time under oxygen limitation, substrate accumulation, actuator saturation, tracking error, and run-to-run robustness. A controller is not accepted solely because it tracks RQ or DO; it must improve a process objective without violating constraints.

Pass **G3** when reconstructed rates and phase changes are reproducible and balance residuals/uncertainty are reported. Pass **G4** when the proposed controller outperforms the documented baseline under nominal and disturbed simulations and remains within actuator and biological limits. Any wet-lab deployment requires a separate approved operating protocol.

## Workstream C — whole-process and commercial boundary

Create `docs/whole_process_boundary.md` with this initial boundary:

```text
feedstock production/acquisition
→ transport, storage, and pretreatment
→ medium and sterile-feed preparation
→ inoculum and aerobic fed-batch fermentation
→ cell/product recovery and lipid extraction
→ product finishing
→ wastewater, off-gas, residual biomass, and coproduct handling
```

Inventory each stage using a common basis, initially **1 kg saleable lipid product**. At minimum track:

- feedstock composition, impurities/inhibitors, variability, price, pretreatment loss, and sterilization demand;
- seed and production volume, batch time, yield, titer, productivity, contamination/failure allowance;
- air/O2, agitation and cooling energy, water, nitrogen source, minerals, acid/base, and antifoam;
- cells, residual carbon/nitrogen, citrate/other byproducts, CO2, heat, and wastewater load;
- recovery/extraction yield, solvent and energy demand, product specification, and coproduct credit;
- scale-dependent oxygen-transfer and heat-removal ceilings.

Use a transparent screening equation before a detailed TEA:

```text
minimum selling price
= (annualized capital + annual operating cost - coproduct credit)
  / annual saleable product
```

Report sensitivity to feedstock price, fermentation yield/titer/productivity, oxygen and power demand, batch turnaround, and recovery yield. Pass **G5** when the top cost/scale bottlenecks remain identifiable across plausible low/base/high assumptions and the system boundary prevents burden shifting.

## Workstream D — paper and defense deliverables

The paper should contain, at minimum:

1. **Mechanism figure:** compartments, enzymes/subunits, cofactors, verified control points, and model reaction IDs.
2. **Evidence table:** each central claim linked to biochemical, genetic, modeling, and fermentation evidence with confidence and contradictions.
3. **Reactor/control figure:** sensor → estimator → controller → actuator → reactor, with sample times and constraints.
4. **Dynamic results figure:** feed, DO, OTR/OUR, CTR/CER, RQ, biomass, substrate, and lipid rate on synchronized time axes with phase annotations.
5. **Model-validation table:** balances, blocked reactions, essentiality, bypasses, parameter sources, uncertainty, and domain of validity.
6. **Whole-process figure/table:** material/energy boundary, key inputs/outputs, bottlenecks, and economic sensitivity.
7. **Claim-to-evidence appendix:** enough detail for another researcher to reproduce every derived rate and model constraint.

Pass **G6** when each main conclusion is supported by at least two independent evidence types where possible, model-only claims are labeled as predictions, and the conclusion remains true under the reported uncertainty/sensitivity analysis.

## Staged execution plan

| Stage | Deliverable | Evidence gate | Stop/go decision |
|---|---|---|---|
| 0. Baseline | This roadmap; fingerprints; limitation register | Repository and model state are reproducible | Do not use an untracked model state in later results |
| 1. Mechanism | First-priority mechanism-and-capacity map | G1–G2 | Do not build control claims on unresolved central GPR/stoichiometry |
| 2. Reactor replay | Signal dictionary, data inventory, reconciled rates, phase map | G3 | If rates are not identifiable, improve measurements before controller complexity |
| 3. Control simulation | Baseline, feed-forward, PID, and constrained model-based comparison | G4 | Advance only if gains survive delay, noise, uncertainty, and actuator limits |
| 4. Whole process | Mass/energy inventory and screening economics | G5 | Reject improvements that are immaterial or shift the bottleneck downstream |
| 5. Manuscript | Integrated figures, claim ledger, methods, limitations, defense question bank | G6 | Submit only claims within the verified domain of validity |

## First executable experiment plan

### Experiment Overview

- **Title:** Reconstruct metabolic state from one fed-batch fermentation and test whether online signals can constrain the GEM
- **Objective:** produce a time-resolved, unit-consistent set of uptake, respiratory, growth, and lipid-production rates; identify physiological phases; translate measured rates into dynamic model bounds
- **Hypothesis:** reconciled feed and off-gas measurements reduce flux uncertainty and reveal regime changes that a static FBA condition misses
- **Type:** data reconciliation + dynamic simulation

### Inputs

| Input | Status | Success criterion |
|---|---|---|
| Raw reactor historian export | Not yet supplied | Original timestamps and unfiltered values retained |
| Feed and sampling log | Not yet supplied | Volume and substrate mass can be reconstructed |
| Off-gas calibration/flow basis | Not yet supplied | O2 and CO2 molar rates are computable with documented corrections |
| Offline biomass, substrate, byproduct, and lipid data | Not yet supplied | At least enough points to estimate phase-level rates and uncertainty |
| Reactor geometry and operating limits | Not yet supplied | OTR and actuator constraints can be bounded |
| GEM and target medium | Available | Exact model hash and medium are recorded for every run |

### Expected outputs

| Output | Planned location | Success criterion |
|---|---|---|
| Signal dictionary | `docs/reactor_signal_dictionary.md` | definition, unit, calibration, frequency, and lag for every signal |
| Clean analysis table | `data/fermentation/derived/<run_id>_rates.csv` | raw-source lineage and units for every derived column |
| Balance and phase report | `results/process_control/<run_id>_reconciliation.md` | residuals and uncertainty shown; no forced closure hidden |
| Dynamic constraint file | `data/fermentation/derived/<run_id>_dfba_bounds.csv` | time-indexed bounds with sign convention and confidence |
| Baseline simulation | `results/process_control/<run_id>_baseline/` | reproduces stated inputs and reports mismatch, saturation, and uncertainty |

### Analysis and decision rules

- Primary metric: uncertainty and residuals in reconstructed `qS`, OUR, CER/CTR, RQ, growth, and lipid-production rates.
- Secondary metric: reduction in feasible flux ranges relative to static medium bounds without contradicting offline measurements.
- Unit audit: 100% of raw and derived variables must have explicit units and sign conventions.
- No balance will be declared closed if unmeasured products or evaporation can explain the residual; these become bounded unknowns.
- Before analysis, set numerical residual thresholds from sensor accuracy and sampling frequency; do not choose them after seeing controller performance.
- Compare phase calls against at least one independent offline measurement.
- A failed identifiability result is actionable: it defines the next sensor, sample, or calibration requirement.

## Immediate next actions

1. Complete the eight-node mechanistic map, beginning with ACL, ACC, FAS, NADPH supply, and TAG assembly.
2. Audit the highlighted GPR/direction questions against primary biochemical and genetic evidence; make no silent model edit.
3. Obtain one complete fermentation historian export plus the accompanying feed, sample, calibration, and offline assay records.
4. Freeze a raw-data manifest and write the signal dictionary before deriving OTR/CTR/RQ.
5. Reconcile one run and establish the existing feed/DO-control baseline before proposing a more advanced controller.
6. Start the whole-process inventory in parallel with explicit unknowns; use it to decide which biological improvements are economically worth optimizing.
