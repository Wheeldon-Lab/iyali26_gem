# Reactor signal and control dictionary

## Material Passport

- **Date:** 2026-07-20
- **Version:** `reactor_signals_v0.1`
- **Verification status:** DESIGN SPECIFICATION — equations and conventions are defined; no historian export has yet been reconciled
- **Target process:** aerobic *Yarrowia lipolytica* batch/fed-batch lipid production
- **Required next input:** raw historian, feed/sample log, gas and probe calibrations, reactor geometry, and offline biomass/substrate/lipid measurements

## Control architecture

```text
raw sensors
  pH, DO, temperature, pressure, gas flow, O2/CO2, mass/level, pump state
        ↓ calibration, time alignment, filtering, delay correction
derived rates / state estimator
  OTR, OUR, CTR/CER, RQ, qS, μ, carbon balance, metabolic phase
        ↓
outer-loop optimizer or supervisory controller
  growth/lipid objective, feed target, oxygen-demand target, constraints
        ↓
inner regulatory loops
  feed pump; pH acid/base; agitation; air; O2 enrichment; pressure
        ↓
reactor and cells ────────────────────────────────────────────────┘
```

The inner loops stabilize physical variables. The outer loop should optimize a biological/process objective. Keeping DO or RQ at a setpoint is not, by itself, proof of improved lipid production.

## Units and sign conventions

Use one internal unit system:

- time: `h`;
- liquid volume: `L`;
- biomass: `gDCW/L`;
- liquid concentrations: `mol/L` or `g/L`, never mixed without molecular-weight conversion;
- total gas rates: `mol/h`; volumetric rates: `mmol/L/h`; specific rates: `mmol/gDCW/h`;
- gas composition: dry mole fraction unless explicitly tagged `wet`;
- process OUR and substrate uptake are reported as positive consumption rates;
- process CER/CTR and product rates are reported as positive production rates;
- GEM mapping follows COBRA convention: substrate and O2 uptake are negative exchange fluxes; CO2/product secretion is positive.

Every column name should carry a machine-readable unit, or its metadata must do so. Never store `OTR`, `CTR`, `feed`, or `air` without specifying whether it is total, volumetric, specific, mass, or molar.

## Raw online signals

| Signal | Canonical name | Unit | Calibration/context required | Typical role |
|---|---|---:|---|---|
| Timestamp | `time_h` plus original timestamp | h | clock, timezone, historian resampling | common independent variable |
| pH measurement/setpoint | `pH`, `pH_sp` | 1 | two-point calibration, temperature compensation, probe age | controlled constraint |
| Acid/base addition | `acid_L_h`, `base_L_h`, cumulative equivalents | L/h; mol H+/h | reagent identity/concentration, pump calibration | actuator and possible acid/base soft signal |
| Dissolved O2 | `do_pct_air_sat`, preferably `cO2_mmol_L` | % or mmol/L | zero/span, temperature, pressure, salinity, probe time constant | controlled variable; oxygen-balance state |
| Temperature | `temp_C`, `temp_sp_C` | °C | sensor calibration | affects kinetics, gas solubility, `kLa` |
| Agitation | `rpm`, `rpm_sp` | min⁻¹ | impeller geometry, torque/power if available | DO-cascade actuator |
| Gas flows | `air_sL_min`, `o2_sL_min`, `n2_sL_min` | standard L/min | definition of standard T/P, mass-flow calibration | gas-balance input and actuator |
| Inlet/outlet O2 and CO2 | `yin_O2`, `yout_O2`, `yin_CO2`, `yout_CO2` | mol/mol dry | analyzer span/zero, water removal, sample-line lag | OUR/CTR calculation |
| Pressure | `p_abs_bar` | bar abs | absolute, not gauge | gas conversion and O2 saturation |
| Feed pump | `feed_L_h`, `feed_kg_h`, cumulative feed | L/h or kg/h | gravimetric pump calibration, density, concentration | principal carbon actuator/input |
| Reactor mass/volume | `mass_kg`, `volume_L` | kg; L | tare, density, evaporation and sampling | dilution and all volumetric rates |
| Antifoam and samples | event log plus amount | mL or g | exact time and composition | volume/property disturbance |
| Controller states | mode, output, saturation, alarms | categorical; % | manual/auto/cascade state and output limits | reconstruct actual control action |

## Offline and at-line measurements

At minimum record DCW, viable biomass if relevant, residual carbon substrates, nitrogen/ammonium, citrate and other major organic acids, lipid concentration, lipid fraction of DCW, fatty-acid composition, and sample volume. Store assay uncertainty and replicate information.

Separate three lipid quantities:

```text
lipid content       = g lipid / gDCW
lipid titer         = g lipid / L broth
lipid productivity  = Δ(g lipid) / (L · h)
```

An increase in content caused by loss of non-lipid biomass is not automatically an increase in titer or productivity.

## Gas and oxygen rates

### Gas-phase balance

If inlet and outlet dry molar flow are measured:

```text
OURtotal = Fin · yO2,in  - Fout · yO2,out
CTRtotal = Fout · yCO2,out - Fin · yCO2,in
OURvol   = OURtotal / V
CTRvol   = CTRtotal / V
```

If outlet flow is not measured and an inert-gas balance is valid:

```text
Fout = Fin · yinert,in / yinert,out
```

The inert fraction must be recomputed when O2 enrichment or N2 addition is used. Correct all flows to the same dry/wet and temperature/pressure basis. Account for analyzer transport delay before combining gas signals with feed or DO.

### Liquid oxygen balance

```text
OTR = kLa · (C*O2 - CL)

d(CL V)/dt
  = OTR · V - OUR · V
    + F · C_O2,feed - Fsample · CL
```

At a stable DO setpoint with small dissolved-O2 accumulation, `OTR ≈ OUR`. During a DO transient, gas-flow step, rapid feed change, or scale-up gradient, the two are not interchangeable.

`C*O2` depends on temperature, pressure, broth composition, and inlet O2 fraction. `kLa` depends on agitation, gas flow, viscosity, antifoam, volume, geometry, and cell/broth state; a water calibration at one operating point is not a universal process constant.

## CO2 and respiratory quotient

Define two ratios:

```text
biological RQ = CER / OUR
transfer quotient = CTR / OTR
```

They approach one another only when gas/liquid accumulation is negligible and all rates share a molar basis. Dissolved CO2/bicarbonate can accumulate or be released when pH changes. Therefore off-gas `CTR/OTR` must not be labeled a true biological RQ during rapid pH or feed transients without a dissolved-inorganic-carbon correction.

This distinction is especially important when on/off acid/base dosing makes CTR noisy. Preserve the raw pH-controller output so that apparent RQ events can be tested against dosing events.

Do not import a numerical RQ setpoint from *S. cerevisiae* into *Y. lipolytica*. Establish a strain-, substrate-, phase-, and objective-specific range using carbon balances and offline measurements.

## Feed, volume, substrate, and biomass rates

For a single-substrate fed batch:

```text
dV/dt      = F - Fsample - Fevap
d(XV)/dt   = μ X V
d(SV)/dt   = F Sf - qS X V - Fsample S

qS = [F Sf - d(SV)/dt - Fsample S] / (X V)
μ  = d ln(XV) / dt
```

Add every carbon-containing cosubstrate and base/antifoam contribution when material. If residual substrate is sparsely sampled, report interval-average `qS` with uncertainty rather than an apparently continuous exact rate.

An exponential feed-forward starting point is:

```text
F(t) = μset X0 V0 exp(μset t) / (YX/S Sf)
```

It assumes constant yield, negligible maintenance/product formation, accurate initial biomass, and no oxygen constraint. These assumptions usually degrade during the growth-to-lipid transition, so the policy needs estimator feedback or scheduled parameters.

## Carbon and electron reconciliation

Use carbon-moles, not only grams:

```text
carbon in feed
= carbon in new biomass + lipid + measured byproducts + CO2
 + residual carbon accumulation + sampled carbon + unmeasured residual
```

Report the unmeasured residual and its uncertainty. Do not force closure by adjusting one measured signal without an independent calibration reason.

Use degree-of-reduction/electron balances to distinguish a carbon-balanced but redox-inconsistent reconstruction. This is the bridge between off-gas rates and the NADPH/respiration mechanisms in `mechanistic_control_map.md`.

## Phase labels and evidence

| Phase/state | Candidate online signatures | Independent confirmation required |
|---|---|---|
| Exponential growth | rising OUR/CTR; stable specific rates; feed/growth relation | DCW and substrate slope |
| Nitrogen limitation onset | ammonium approaches limit; growth slope changes; gas-rate inflection | ammonium assay and biomass composition |
| Lipid accumulation | lipid rate increases; growth decouples from carbon uptake; phase-specific RQ/OUR pattern | lipid titer/content and fatty-acid analysis |
| Carbon overfeed | residual substrate rises; OUR/feed decouple; possible citrate/byproduct shift | residual substrate and citrate |
| Oxygen-transfer limitation | DO at lower bound; cascade saturated; OTR ceiling; feed increase no longer raises OUR | `kLa`/gas balance and actuator saturation |
| Carbon depletion | sharp OUR/CTR/feed-response change | residual substrate and cumulative mass balance |

Online signatures are hypotheses until confirmed. A threshold selected on one run must be tested on independent runs.

## Regulatory loops

### pH loop

- controlled variable: pH;
- manipulated variables: acid/base pumps;
- mandatory constraints: pump limits, deadband, minimum pulse, mixing delay;
- metabolic use: cumulative equivalents may indicate changing acid/base production only after buffer, CO2 chemistry, feed, and sampling are included.

### DO cascade

A common priority is agitation → air flow → O2 enrichment → pressure, with each actuator clamped to validated bounds. The exact order must reflect equipment, shear, foam, power, and cost constraints. Record which actuator is active and saturated.

The DO controller can hide a rising OUR by increasing oxygen supply. Therefore interpret DO together with controller output, OTR/OUR, and saturation state.

### Feed outer loop

Candidate policies, in order of complexity:

1. fixed or scheduled feed;
2. exponential feed-forward;
3. DO-stat or RQ-stat feedback with explicit guard conditions;
4. estimated-`qS` or estimated-growth feedback;
5. constrained model-predictive control using the dynamic balance/GEM state.

All policies require hard bounds on feed, volume, residual substrate, DO, gas flows, agitation, pressure, and any product-quality constraint.

## Controller mathematics and validation ladder

### 1. Identify dynamics

Apply safe small input steps or use existing excitation to estimate gain, time constant, delay, and interaction for feed→OUR/RQ/substrate and gas/agitation→DO/OTR. Report identifiability and uncertainty.

### 2. Regulatory PID

```text
e(t) = ysp(t) - y(t)
u(t) = Kp [e(t) + (1/Ti)∫e(t)dt + Td de(t)/dt]
```

Use derivative filtering, anti-windup, output/rate limits, bumpless manual/auto transfer, and alarm/interlock behavior. A clean nominal trace without saturation/noise tests is insufficient.

### 3. Supervisory/MPC layer

A generic finite-horizon objective is:

```text
min Σ [wtrack(y-ysp)^2 + wmove(Δu)^2 - wprod·lipid_rate]
subject to dynamic balances, GEM/enzyme constraints,
           state bounds, actuator bounds, and rate-of-change bounds
```

State estimation may progress from filtered algebraic rates to an extended/unscented Kalman filter or moving-horizon estimator. The estimator must be validated separately from the controller.

### 4. Required comparisons

Compare open-loop recipe, documented current control, feed-forward, feedback, and model-based control under the same initial states and disturbances. Test:

- ± uncertainty in yield, maintenance, biomass, `kLa`, and kinetic/capacity parameters;
- gas-analyzer and DO-probe delay/noise/drift;
- feed-concentration and pump-calibration error;
- actuator saturation, lost signal, sample removal, and antifoam events;
- scale-dependent mixing/oxygen-transfer limits.

Primary outcomes: lipid titer, yield, productivity, citrate/byproducts, time under oxygen limitation, residual substrate, actuator saturation, power/O2 use, and run-to-run robustness.

## Data-table contract

Use a long-form raw table or a lossless wide export, plus a separate variable dictionary. A derived table should minimally contain:

```text
run_id, time_h, phase, volume_L, biomass_g_L,
feed_kg_h, substrate_g_L, qS_mmol_gDCW_h,
do_pct, kla_h_1, otr_mmol_L_h, our_mmol_L_h,
ctr_mmol_L_h, cer_mmol_L_h, rq_bio, tq_transfer,
pH, acid_mol_h, base_mol_h, rpm, air_sL_min, o2_sL_min,
citrate_g_L, lipid_g_L, lipid_g_gDCW,
quality_flag, uncertainty_note, source_columns
```

Use `NA`, never zero, for unavailable measurements. Quality flags must identify calibration periods, manual sampling, analyzer cleaning, controller manual mode, saturation, and interpolation.

## Acceptance gates

- **S1 — traceability:** every derived value maps to raw columns, calibration, equation, unit, and software version.
- **S2 — timing:** sensor and actuator delays are estimated or bounded; raw and corrected times are retained.
- **S3 — balances:** volume, gas, carbon, and electron residuals are reported with uncertainty and unmeasured pools.
- **S4 — state validity:** phase/rate estimates agree with at least one independent offline measurement.
- **S5 — control validity:** gains persist under noise, delay, uncertainty, and actuator constraints and improve lipid titer/yield/productivity rather than only tracking.
- **S6 — deployment safety:** any reactor implementation has an independently reviewed operating envelope, fallback mode, and abort criteria.

## Sources used

- Kavšček et al. (2015), [GEM-guided feed and aeration strategy in *Y. lipolytica*](https://pmc.ncbi.nlm.nih.gov/articles/PMC4623914/).
- De Biaggi et al. (2026), [phase-resolved bioreactor exchange rates and enzyme-constrained modeling](https://research.chalmers.se/publication/551434/file/551434_Fulltext.pdf).
- Royce and Thornhill (1992), [why CTR/OTR can differ from biological RQ during pH and CO2 transients](https://pubmed.ncbi.nlm.nih.gov/18601064/).
- Löser et al. (2021), [online OUR/CER and momentary RQ calculation in an aerobic bioreactor](https://pmc.ncbi.nlm.nih.gov/articles/PMC7923609/).
- Craven et al. (2014), [oxygen-transfer characterization and DO-control capacity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3790518/).
