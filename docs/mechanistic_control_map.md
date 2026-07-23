# Mechanistic control map for lipid production

## Material Passport

- **Origin:** committee-aligned mechanism audit
- **Date:** 2026-07-20
- **Version:** `mechanistic_map_v0.1`
- **Verification status:** PARTIALLY VERIFIED — current SBML mappings were inspected directly; external statements below are tied to primary papers where available; kinetic parameters have not been reproduced locally
- **Model baseline:** `M_iYali26`, 2,300 reactions, 1,073 genes; `model.xml` SHA-256 `dee22ccb2febe39282ac4fae240b5d56c95a23ebaa2b262702c7b7de0c3c012b`
- **Naming convention:** `YALI1*` denotes the current W29/CLIB89 mapping in this repository; `YALI0*` denotes the legacy CLIB122 loci commonly used in the literature
- **Change policy:** this document records hypotheses and conflicts. It does **not** authorize silent edits to the GEM.

## The mechanistic story to be able to defend

```text
carbon uptake
→ glycolysis or glycerol assimilation
→ mitochondrial citrate production
→ citrate export
→ ACL1–ACL2 supplies cytosolic acetyl-CoA
→ ACC1 supplies malonyl-CoA
→ FAS1–FAS2 consumes acetyl-CoA, malonyl-CoA, and NADPH to form acyl-CoA
→ GPD1/SCT1/SLC1/PAH1 builds the glycerolipid backbone
→ DGA1/DGA2 or LRO1 converts DAG to TAG
→ TAG is stored in lipid particles

competing processes:
  citrate secretion; amino-acid and biomass synthesis; sterol/phospholipid synthesis;
  TAG lipolysis; peroxisomal POX1–POX6 β-oxidation; ATP/redox maintenance
```

Nitrogen limitation changes allocation among these processes; it is not sufficient to say that “nitrogen limitation turns lipid synthesis on.” In a chemostat multi-omics study, the lipid phenotype was associated with carbon reallocation away from amino-acid biosynthesis rather than a general transcriptional upregulation of lipid-synthesis genes. In the 2026 enzyme-constrained study, exponential and nitrogen-limited phases produced different enzyme-capacity targets. The mechanism is therefore condition- and phase-specific.

## Four meanings of “control”

These must not be used interchangeably in a paper or defense:

1. **Structural requirement:** all necessary subunits/accessory proteins are present; represented by a justified GPR.
2. **Reaction capacity:** the maximum enzyme-supported rate, approximately

   ```text
   vi ≤ kcat,i · Ei
   ```

   after unit, oligomer, saturation, and isozyme corrections.
3. **Local kinetic response:** substrate, product, cofactor, allosteric state, pH, and temperature change the fraction of capacity used.
4. **System-level flux control:** changing one enzyme changes the pathway output. A local approximation is

   ```text
   CEi^J = ∂ln(J) / ∂ln(Ei)
   ```

   A large flux or abundant protein does not by itself imply a large flux-control coefficient.

## First-pass evidence matrix

| Node | Confirmed biological role | Current model representation | Current conflict or uncertainty | Reactor-scale consequence to test | Evidence status |
|---|---|---|---|---|---|
| ATP-citrate lyase (ACL1–ACL2) | Heteromeric enzyme supplies cytosolic acetyl-CoA from citrate; both subunits form the catalytic architecture | `R1894`; `YALI1E41315g or YALI1D32268g` | The biological subunit logic indicates `AND`, not the present `OR`. The 2026 ecGEM also treats the two legacy genes together and predicts high lipid flux control, but that FCC is model-derived | A capacity limit should divert citrate toward secretion and reduce lipid synthesis; compare citrate, lipid rate, and respiratory demand | Subunit architecture E2; capacity-control claim provisional E3-model |
| Acetyl-CoA carboxylase (ACC1) and biotin activation | ACC1 converts cytosolic acetyl-CoA to malonyl-CoA; BPL1/holocarboxylase synthetase enables biotin-dependent carboxylase activity | `R88`: `YALI1C15991g and YALI1E35955g`; `R175` represents biotin-AMP formation; `R87` is a mitochondrial ACC-like reaction | Need distinguish catalytic ACC1 from the upstream BPL1 dependency and verify whether the two-reaction representation preserves biotin catalysis without making BPL1 an artificial stoichiometric subunit | A malonyl-CoA limit should reduce fatty-acid synthesis even if acetyl-CoA/citrate is available; biotin availability becomes a medium/control covariate | Function E2; current kinetic capacity E0 |
| FAS1–FAS2 fatty-acid synthase | Native type-I FAS contains six FAS2 α and six FAS1 β subunits; consumes malonyl-CoA and NADPH to form fatty acyl products | `R1392`: `YALI1B19844g and YALI1B25427g`; `R1393` additionally requires `YALI1C15991g and YALI1E27279g` | Core FAS `AND` is supported. `R1393` uniquely adds ACC1 plus acyl-CoA-binding protein 2 even though malonyl-CoA is already a substrate; C16/C18 GPR asymmetry requires mechanistic justification | FAS capacity should affect lipid-production rate and acyl-chain profile; the 2026 ecGEM predicts the C18 step has the largest lipid FCC | Oligomer E2; FCC provisional E3-model; current `kcat` absent |
| Cytosolic NADPH supply | On glucose, 13C-MFA supports oxidative PPP as the principal lipogenic NADPH source; malic-enzyme flux did not rise with the doubled fatty-acid flux | `R325` ZWF1/G6PDH (`YALI1E26811g`), `R639` GND1 (`YALI1B20462g`); `R538` is a mitochondrial NADP malic enzyme | Do not cite mitochondrial `R538` as the default cytosolic lipogenic NADPH source. Quantify PPP, isocitrate-dehydrogenase, and any modeled redox-transfer routes separately, especially on glycerol | PPP/redox limitation may change CO2 production, oxygen demand, growth/lipid yield, and oxidative-stress response; relationships depend on substrate | Glucose 13C-MFA E4; glycerol/other substrates require new evidence |
| Glycerol-3-phosphate supply and shuttle | Cytosolic GPD1 controls glycerol-3-phosphate backbone supply; GUT2 and the shuttle couple this pool to redox metabolism; perturbations affect TAG accumulation | `R348` GPD1 (`YALI1B04433g`); `R347/R349` mitochondrial dehydrogenases; `R1142` is a reversible G3P “transport” with no GPR | `R1142` may be a bookkeeping representation of the redox shuttle rather than a literal transport step. Its thermodynamic and compartment logic must be tested before assigning a transporter or kinetic bound | GPD1/GUT2 activity can couple NADH balance, respiration, and TAG backbone supply; test glycerol-3-P, RQ, OUR, and lipid rate together | Genetic perturbation E4; transport representation E0 |
| Glycerolipid assembly (SCT1/SLC1/PAH1) | Sequential acylation produces lysophosphatidate and phosphatidate; PAH1 supplies DAG for TAG formation | `R350–R353` map one acyltransferase to ER/lipid-particle reactions; downstream SLC1/PAH1 reactions need a complete reaction-ID audit | `R352` is disabled while related chain/compartment variants are active; chain specificity, compartment transfer, and aggregate lipid pools may create artificial bottlenecks | Backbone limitation should accumulate acyl-CoA or PA intermediates and alter TAG/phospholipid partitioning, not necessarily total fatty-acid synthesis | Annotation E1–E2; pathway capacity not established |
| TAG synthases DGA1, DGA2, LRO1 | DGA1/DGA2 use acyl-CoA; LRO1 uses phospholipid as acyl donor. Deletion/heterologous-expression studies confirm distinct contributions and carbon-source dependence | DGA1 `R1771–R1772` (`YALI1E38810g`); DGA2 `R2165/R2168` (`YALI1D10264g`); LRO1 `R1728–R1729` (`YALI1E20049g`) | DGA1 and LRO1 reactions are reversible in places despite their storage-lipid role; `R1772` is disabled; DGA2 has both a mechanistic ER reaction and a lumped cytosolic TAG synthesis reaction, creating possible duplication | TAG assembly capacity should change fatty-acid-to-TAG partitioning, lipid-droplet formation, and free-acyl stress; distinguish lipid content from de novo fatty-acid rate | Enzyme function/perturbation E4; current direction/duplication E0–E1 |
| Peroxisomal β-oxidation (POX1–POX6) | Six acyl-CoA oxidases catalyze the first, rate-limiting β-oxidation step with different chain preferences; POX2 and POX3 are prominent long- and short-chain activities | `R97–R102`, `R1487–R1498`; current GPR set uses four identifiers, including legacy `YALI0D27654g`, across many chain lengths | Reconcile all six current loci, chain specificity, and induction. A single enzyme is currently assigned across C10–C26 in `R97–R102`, while several unsaturated/long-chain reactions use broad OR rules | β-oxidation competes with storage and changes O2 demand and CO2 production; the effect depends strongly on whether carbon enters as sugar/glycerol or fatty acid | Isozyme/substrate evidence E3–E4; present GPR coverage incomplete |

## Model conflict register

| ID | Priority | Finding | Required resolution | Acceptance evidence |
|---|---:|---|---|---|
| MC-ACL-01 | Critical | `R1894` uses ACL1 `OR` ACL2 although ACL is heteromeric | Compare iYali4/iYli21/eciYali5 GPRs and direct subunit evidence; stage an `AND` correction only after model feasibility and essentiality tests | Primary source + cross-model provenance + before/after growth, flux, and essentiality report |
| MC-FAS-01 | Critical | `R1393` makes ACC1 and acyl-CoA-binding protein 2 obligatory for only the C18 FAS reaction | Separate catalytic subunits, precursor-producing reactions, and accessory binding proteins; verify why C16 and C18 differ | Reaction mechanism, original reconstruction provenance, and no artificial gene lethality |
| MC-TAG-01 | High | DGA1/LRO1 reversibility, disabled variants, and DGA2 lumping may duplicate or reverse TAG synthesis | Reconstruct a compartment- and chain-aware DAG→TAG route with explicit physiological direction | Element/charge balance, thermodynamic rationale, lipid-production feasibility, no duplicate free route |
| MC-NADPH-01 | High | Mitochondrial malic enzyme is easy to misinterpret as cytosolic lipogenic NADPH supply | Produce compartment-specific NADPH production/consumption tables under glucose and glycerol | Redox balance plus 13C/proteome/exchange-data consistency |
| MC-G3P-01 | Medium | `R1142` directly transports G3P and has no GPR | Decide whether it is a bookkeeping shuttle, a transporter, or should be represented as coupled redox reactions | Compartment mechanism and no free redox/energy cycle |
| MC-POX-01 | High | POX chain specificity is compressed into incomplete/broad GPRs and one legacy locus | Map POX1–POX6 current loci and substrate ranges; distinguish verified activity from permissive model coverage | Locus provenance, biochemical substrate evidence, and fatty-acid-growth validation |

## Initial capacity and control claims

| Claim | What the evidence actually supports | What it does not yet support |
|---|---|---|
| FAS is a central lipid control point | In eciYali5, a 1% `kcat` perturbation analysis ranked the C18 FAS1–FAS2 step highest for TAG FCC under the stated glycerol-bound optimization; FAS oligomer structure is independently supported | That FAS always controls lipid production in this repository, every strain, or every process phase |
| ACL has substantial lipid control | The same ecGEM analysis assigned high FCC to ACL, consistent with precursor-supply logic | That either subunit alone is sufficient, or that ACL overexpression will outperform feed/oxygen interventions |
| PPP supplies lipogenic NADPH | 13C-MFA on glucose found oxidative-PPP NADPH production tracked the increase in fatty-acid synthesis while malic-enzyme flux did not | That the exact same fractional contribution holds on glycerol, fatty acids, or all nitrogen-limited states |
| DGA1/LRO1 control TAG storage | Gene deletions, heterologous complementation, and enzyme assays support distinct acyl donors and condition-dependent contributions | That reversible GEM reactions or a single lumped TAG equation reproduce that biology |
| POX controls competing degradation | Genetic and biochemical studies support multiple chain-selective Aox isozymes and a rate-limiting first β-oxidation step | The present broad OR rules or exact quantitative POX capacity in the target reactor condition |

## Measurements needed to move from E2/E3 to E5

For each candidate control point, collect or derive a synchronized evidence tuple:

```text
[gene/protein abundance]
+ [substrate/cofactor state]
+ [estimated intracellular or exchange flux]
+ [online reactor state]
+ [lipid rate/composition]
+ [defined perturbation]
```

Minimum useful comparisons are exponential growth versus nitrogen-limited lipid accumulation, and parental versus high-lipid strain, on the actual target feedstock. A single end-point lipid percentage cannot distinguish precursor, redox, FAS, TAG-assembly, or degradation control.

The first quantitative tests should be:

1. Map current-locus proteins to the absolute proteomics and `kcat` assignments in eciYali5.
2. Recompute enzyme-usage and FCC rankings with the exact target substrate, phase, exchange rates, and biomass/lipid composition.
3. Perform ±10× sensitivity on non-organism-specific or manually relaxed `kcat` values.
4. Compare predicted rate changes against time-resolved glycerol/glucose, O2, CO2, biomass, citrate, and lipid data.
5. Treat disagreement as a localization/GPR/kinetic hypothesis, not as a reason to tune a parameter until the observation fits.

## Sources used in this audit

- De Biaggi et al. (2026), [enzyme-constrained model and growth-phase-specific targets](https://research.chalmers.se/publication/551434/file/551434_Fulltext.pdf).
- Domenzain et al. (2022), [GECKO 2.0 and the limits of cross-species/low-specificity kinetic parameters](https://www.nature.com/articles/s41467-022-31421-1).
- Kerkhoven et al. (2016), [amino-acid regulation and carbon reallocation during lipid accumulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC5516929/).
- Wasylenko, Ahn, and Stephanopoulos (2015), [13C-MFA evidence for oxidative PPP as the primary lipogenic NADPH source on glucose](https://doi.org/10.1016/j.ymben.2015.02.007).
- Athenstaedt (2011), [DGA1 and LRO1 biochemical/genetic characterization](https://pmc.ncbi.nlm.nih.gov/articles/PMC3161177/).
- Dulermo and Nicaud (2011), [G3P shuttle and β-oxidation perturbations](https://pubmed.ncbi.nlm.nih.gov/21620992/).
- Wang et al. (1999), [Aox isozyme functional evaluation](https://pubmed.ncbi.nlm.nih.gov/10464181/).
