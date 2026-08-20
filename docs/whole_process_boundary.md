# Whole-process boundary and commercial-practicality framework

## Material Passport

- **Date:** 2026-07-20
- **Version:** `whole_process_v0.1`
- **Verification status:** SCREENING FRAMEWORK — system boundary and equations are defined; project-specific inventory and prices are not yet supplied
- **Default functional unit:** 1 kg saleable lipid product at the stated specification and plant gate
- **Scenario status:** no feedstock, product grade, annual capacity, geography, scale, or downstream route has yet been selected

## Decision question

Can a defined feedstock be converted by the target *Yarrowia lipolytica* strain and control strategy into a defined lipid product at a yield, titer, productivity, recovery, energy demand, and risk profile that is competitive for the intended market?

The words **defined feedstock** and **defined product** are mandatory. “Microbial lipid” is not a single commercial product: fuel-range TAG, food fat, oleochemical intermediate, and high-value structured lipid have different purity, regulation, recovery, price, and scale requirements.

## System boundary

```text
feedstock source
  production or waste generation
  collection, transport, storage
        ↓
pretreatment and conditioning
  separation, hydrolysis, detoxification, concentration, sterilization
        ↓
medium and feed preparation ───── seed train
        ↓                            ↓
production fermentation
  carbon/nitrogen feed, air/O2, agitation, cooling, pH reagents, antifoam
        ↓
broth harvest and concentration
        ↓
cell disruption → lipid extraction → solvent recovery → product finishing
        ↓
saleable lipid product

side streams:
  CO2/off-gas, citrate/polyols, residual carbon/nitrogen,
  aqueous waste, cell debris, spent solvent, heat, and off-spec product
```

The initial boundary ends at the plant gate. Feedstock upstream burdens and coproduct/waste treatment remain inside the inventory even if performed off-site. Transport to the final customer can be added as a separate scenario.

## Functional unit and product specification

Before comparing scenarios, complete:

| Field | Required definition |
|---|---|
| Product | TAG mixture, specific fatty acid profile, biodiesel, food fat, oleochemical, or other |
| Purity/composition | total lipid purity; TAG/FFA/sterol fractions; target acyl-chain distribution; contaminants |
| Physical form | crude oil, refined oil, esterified fuel, powder/encapsulated product |
| Regulatory market | fuel, feed, food, cosmetic, chemical, research reagent |
| Functional unit | default 1 kg saleable product; add energy or functional-performance basis if needed |
| Annual capacity | kg or tonne product/year with uptime basis |
| Geography/year | currency year, utility prices, feedstock availability, regulation, and logistics |

No minimum selling price should be quoted without these fields.

## Feedstock scenario matrix

At least four candidate classes should be screened, even if only one advances:

| Feedstock class | Potential advantage | Required penalties/risks to quantify | Key characterization |
|---|---|---|---|
| Refined glucose/sugar | consistent and easy to sterilize/model | purchase cost; land/upstream burden; competing use | sugar spectrum, concentration, price, seasonal supply |
| Purified or crude glycerol | biodiesel coproduct; established *Yarrowia* use | methanol, salts, soaps, metals, water, lot variability; purification loss | carbon purity, ash, methanol, fatty matter, density, viscosity |
| Waste hydrolysate/molasses | low-cost carbon and circularity potential | inhibitors, variable C/N, color/solids, detoxification and wastewater | fermentable carbon, inhibitors, COD, N/P, suspended solids |
| Volatile fatty acids | integration with anaerobic food-waste processing | acid composition, pH/base demand, toxicity, upstream digester yield | acetate/propionate/butyrate, salts, COD, variability |
| Lipid/oil substrate | high carbon conversion to lipid may be possible | emulsification, oxygen/mixing, cost, imported fatty-acid profile, attribution of de novo lipid | TAG/FFA profile, impurities, emulsifier and transfer behavior |

For each lot, distinguish **nominal carbon concentration** from **bioavailable carbon**. A low or negative gate price does not imply a low delivered fermentation cost if concentration, detoxification, sterilization, or wastewater treatment is expensive.

## Block-level inventory

### 1. Feedstock acquisition and conditioning

Inputs: raw feedstock, transport fuel, water, heat/electricity, chemicals, filters/adsorbents, storage.

Outputs: conditioned carbon feed, rejected solids/aqueous fraction, inhibitor-rich waste, losses.

Required measurements: delivered composition distribution, carbon recovery, dilution/concentration factor, sterility strategy, storage loss, and per-lot fermentation performance.

### 2. Medium, feed, and seed train

Inputs: conditioned carbon, nitrogen source, minerals/vitamins, water, pH reagent, steam or alternative sterilization, inoculum vessels.

Outputs: sterile feeds, inoculum, cleaning/sterilization waste.

Required measurements: medium mass per batch, seed-to-production ratio, seed time, contamination/failure probability, and clean-in-place/steam-in-place demand.

### 3. Production fermentation

Inputs: inoculum, feeds, air/O2, agitation power, cooling, acid/base, antifoam, water.

Outputs: biomass/lipid broth, CO2, heat, citrate/polyols/other byproducts, residual nutrients.

Required dynamic inventory: integrated feed mass, O2 consumed, CO2 evolved, agitator power, O2 enrichment, cooling duty, batch/turnaround time, samples and evaporative loss.

The reactor result is summarized by:

```text
YP/S = kg lipid formed / kg bioavailable substrate consumed
titer = kg lipid / m3 final broth
Qp = kg lipid / (m3 · h process time)
lipid content = kg lipid / kg DCW
```

Use full process time in productivity, including fill, sterilization, seed/lag if allocated, harvest, cleaning, and turnaround when estimating annual output.

### 4. Harvest, disruption, extraction, and finishing

Inputs: centrifugation/filtration power, wash water, disruption energy, solvent or supercritical fluid, heat, refining reagents.

Outputs: recovered lipid, wet/dry cell debris, solvent loss, aqueous waste, off-spec fractions.

Required measurements: cell recovery, solids concentration, disruption yield, extraction yield, solvent recovery, product loss, purity, and number of recycle stages.

Because the lipid is intracellular, fermentation titer cannot be treated as final product recovery. Overall recovery is multiplicative:

```text
ηoverall = ηharvest · ηdisruption · ηextraction · ηrefining
saleable product = fermentation lipid · ηoverall
```

### 5. Coproducts and waste

Do not credit citrate, polyols, cell debris, or recovered salts until a sellable specification, recovery operation, market capacity, and price are defined. Otherwise report them as zero-credit sensitivity cases.

CO2, high-COD water, residual nitrogen/phosphorus, solvents, salts, and genetically modified biomass handling must have treatment routes and costs.

## Minimum mass and energy model

### Annual production

```text
batches_per_year
= operating_hours_per_year / (fermentation_time + turnaround_time)

annual_saleable_product
= working_volume · final_lipid_titer · ηoverall · batches_per_year · success_fraction
```

For continuous operation, replace the batch expression with steady volumetric productivity, operating hours, and availability.

### Feedstock requirement

```text
kg raw feedstock / kg product
= 1 / (bioavailable_carbon_fraction · YP/S · ηpretreatment · ηoverall)
```

Use an elemental carbon balance when comparing chemically different substrates.

### Oxygen and power

```text
mol O2 / batch = ∫ OURtotal(t) dt
mol CO2 / batch = ∫ CERtotal(t) dt
agitation energy = ∫ Pshaft(t) dt
```

Separate compressed-air electricity, pure-O2 purchase/generation, agitation, cooling, sterilization, evaporation/concentration, harvest, disruption, extraction, and solvent recovery. Oxygen-transfer and heat-removal capacity are scale constraints, not only operating costs.

## Screening economics

Use a transparent annual cash-cost model before a detailed discounted-cash-flow TEA:

```text
annualized capital = installed capital · capital-recovery factor

minimum selling price
= (annualized capital + fixed OPEX + variable OPEX
   + waste treatment - defensible coproduct credits)
  / annual saleable product
```

Variable OPEX must include feedstock delivered cost, nutrients, water, acid/base, antifoam, air/O2, electricity, steam/heat, cooling, solvent/reagent makeup, waste treatment, and consumables. Fixed OPEX should include labor, maintenance, quality, insurance/overhead, and facility-dependent allocations.

Report low/base/high cases and a break-even value for each uncertain variable. Do not combine optimistic assumptions for every parameter into one “best case” without showing joint probability or scenario consistency.

## Commercial-practicality gates

| Gate | Go criterion | Evidence required | Typical no-go signal |
|---|---|---|---|
| P0 Product/market | product specification and market are defined | customer/specification/regulatory basis | only “microbial lipid” is specified |
| P1 Feedstock | sufficient delivered quantity and acceptable lot variability | composition distribution, logistics, storage, price/credit | cheap average price but inhibitor/seasonal failure |
| P2 Biological performance | yield, titer, productivity, and robustness meet scenario targets | replicated controlled fermentations with full balances | high lipid content but low titer/productivity |
| P3 Scale transfer | OTR, mixing, heat removal, and control actions fit scalable equipment | `kLa`, power, gas, cooling, and gradient analysis | bench result requires infeasible rpm/O2/power |
| P4 Recovery | product can be recovered to specification at defensible yield | integrated harvest/disruption/extraction data | fermentation gain is lost during recovery |
| P5 Economics | MSP is competitive with target price under uncertainty | itemized CAPEX/OPEX, sensitivity and break-even analysis | viability requires simultaneous optimistic assumptions |
| P6 Environment/regulation | burdens and handling routes are acceptable | preliminary LCA inventory, wastewater/emissions, GMO/product pathway | burden shifts to energy, solvent, water, or waste |

## Initial bottleneck hypotheses

These are hypotheses to test, not conclusions:

1. **Commodity versus specialty economics:** a fuel/commodity product requires far lower cost and much larger scale than a high-value composition-specific lipid.
2. **Yield–titer–productivity trade-off:** improving lipid fraction without increasing recovered lipid per reactor volume and time may not reduce cost.
3. **Aerobic scale-up:** oxygen transfer, agitation/compression power, cooling, foam, and mixing gradients may dominate the feasible operating window.
4. **Intracellular recovery:** concentration, disruption, and extraction can erase gains achieved in fermentation.
5. **Feedstock variability:** crude glycerol, molasses, hydrolysates, and VFAs may reduce raw-material cost while increasing conditioning, control, assay, and waste-treatment costs.
6. **Citrate/polyol branching:** coproducts can add value only if their formation, separation, specification, and market are jointly designed; otherwise they reduce carbon yield and add separation load.
7. **Model domain:** a steady GEM without enzyme capacity or validated dynamic exchange rates cannot alone predict plant throughput, oxygen equipment, or control performance.

## Sensitivity and break-even table

Populate this table before making a commercial claim:

| Variable | Low | Base | High | Distribution/source | MSP elasticity or rank |
|---|---:|---:|---:|---|---:|
| Delivered feedstock cost/credit | TBD | TBD | TBD | supplier/market/year | TBD |
| Bioavailable carbon fraction | TBD | TBD | TBD | lot analysis | TBD |
| Fermentation `YP/S` | TBD | TBD | TBD | balanced runs | TBD |
| Final lipid titer | TBD | TBD | TBD | balanced runs | TBD |
| Full-cycle productivity | TBD | TBD | TBD | historian + turnaround | TBD |
| O2 and agitation energy | TBD | TBD | TBD | OUR, gas, power | TBD |
| Overall recovery | TBD | TBD | TBD | integrated downstream run | TBD |
| Batch success/availability | TBD | TBD | TBD | pilot/operational assumption | TBD |
| Installed equipment cost | TBD | TBD | TBD | scaled quotes/correlation | TBD |
| Coproduct credit | 0 | TBD | TBD | specification + market cap | TBD |

Use tornado/spearman/global sensitivity to identify the top three economic levers. Feed those levers back to strain and controller objectives. For example, if oxygen/power dominates, maximizing carbon yield at any OUR is the wrong objective; if recovery dominates, increasing titer or changing lipid localization may matter more than percentage content.

## Integration with the GEM and controller

| Process quantity | Model/control mapping |
|---|---|
| Feed composition and rate | time-varying exchange bounds with impurity/inhibitor scenarios |
| OUR/CER and RQ | respiratory exchange constraints and redox/carbon validation |
| Lipid composition | strain/phase-specific biomass and product equations |
| Enzyme/protein state | ecGEM capacity constraints and phase-specific targets |
| OTR ceiling | dynamic upper bound on O2 uptake; scale-dependent controller constraint |
| Product/recovery value | economic objective or weighting, not a biological flux alone |
| Power/O2/feed costs | controller stage cost and TEA variable OPEX |

The controller should optimize an economically and physically meaningful objective such as recovered lipid productivity minus feed/O2/power penalties, subject to biological and equipment constraints. It should not maximize a lipid flux in isolation.

## Immediate data request

To establish the project base case, obtain:

1. target lipid product and specification;
2. intended feedstock(s), supplier/source, composition data, delivered price or tipping credit, and annual availability;
3. complete fermentation mass/gas/power history and offline product/byproduct data;
4. proposed working volume, annual production, operating days, and scale-up path;
5. actual or candidate harvest, disruption, extraction, solvent-recovery, and refining steps;
6. measured recovery/purity or literature analogs with uncertainty;
7. local utility, labor, waste, and equipment-cost basis with currency year.

## Sources used

- Karamerou et al. (2021), [techno-economic minimum-cost analysis for a microbial palm-oil substitute](https://doi.org/10.1186/s13068-021-01911-3).
- Kumar, Tyagi, and Drogui (2023), [economic analysis of *Y. lipolytica* lipid/citric-acid coproduction on purified crude glycerol](https://doi.org/10.1007/s13399-021-01772-8).
- Rakicka et al. (2015), [industrial byproducts, oxygenation, yield, and productivity in *Y. lipolytica* lipid production](https://pmc.ncbi.nlm.nih.gov/articles/PMC4513389/).
- Pereira, Lopes, and Belo (2023), [crude glycerol/VFA fed-batch and coproduct formation](https://doi.org/10.1039/D3SE00682D).
- Kavšček et al. (2015), [model-guided feed/aeration trade-offs and lipid remobilization](https://pmc.ncbi.nlm.nih.gov/articles/PMC4623914/).
