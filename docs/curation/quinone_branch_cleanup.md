# Quinone branch cleanup, CoQ9 main-chain correction, and reviewed GPRs

Initial cleanup: 2026-08-04

CoQ9 follow-up: 2026-08-12

Step-specific GPR follow-up: 2026-08-14

## Decision

Remove `R189`, `R2242`, `R2247`, `R2248`, `R2249`, and `R2250` as one
atomic curated patch. Keep the mitochondrial `R407` COQ2 route and the generic
SAM/SAH transport reactions `R2243`-`R2246`.

This is an identity and topology correction. It does not add a quinone demand,
open a bound, or attempt to improve essentiality recall.

The six genes that become reaction-orphaned are deliberately retained in the
GEM. Five occur in the positive-only essentiality reference; deleting the gene
objects would move them outside the evaluation intersection and create an
artificial recall increase. They therefore remain explicit unresolved FNs with
zero trusted metabolic reactions.

## Follow-up: balanced CoQ9 main chain

The formal build pipeline now replaces the inherited CoQ6 identities with the
CoQ9 homologous series before FVA. The reviewed route is:

`R763 -> R407 -> R969 (mitochondria to cytosol) -> R39 -> R808 -> R715 -> R40 -> R19 -> R18 -> R695 -> R385`.

Twelve chain-specific metabolites were changed from the C30/CoQ6 series to the
C45/CoQ9 series. Every homolog gains `C15H24`; charge is unchanged. Only two
reaction equations change:

- `R763`: `4 IPP + pentaprenyl diphosphate -> 4 PPi + nonaprenyl diphosphate`.
- `R385`: `SAM + 3-demethylubiquinone-9 -> SAH + ubiquinone-9`.

All 11 reactions in the synthesis route are now element- and charge-balanced.
The four-IPP lump in `R763` is balanced model bookkeeping, and the oxidized
`R385` endpoint is a declared model convention; neither is claimed as a direct
native Yarrowia reaction measurement. The transport steps and their
compartments also remain model assignments.

The CoQ9 chemistry patch itself changes no GPR, bound, compartment, biomass
coefficient, demand, sink, or model object count. It therefore corrects
chemical identity without claiming that the pathway is active.

## Follow-up: step-specific GPR decomposition

The formal pipeline now replaces the inherited seven-gene `AND` repeated on
five downstream reactions with the following reviewed, step-specific rules:

| Reaction | Formal GPR | Interpreted step | Evidence boundary |
|---|---|---|---|
| `R715` | `YALI1B20835g` | COQ3 first O-methylation | Cross-species biochemical/structural evidence plus compatible predicted fold; native Yarrowia biochemistry unverified. |
| `R385` | `YALI1B20835g` | COQ3 terminal O-methylation | Same boundary; the native substrate redox state remains unresolved. |
| `R18` | `YALI1C25352g` | COQ5 C-methylation | Cross-species biochemical/structural evidence plus compatible predicted fold; native Yarrowia biochemistry unverified. |
| `R695` | `YALI1E18269g` | COQ7 hydroxylation | Cross-species biochemical/structural evidence plus compatible predicted fold; native Yarrowia biochemistry unverified. |
| `R40` | `YALI1F34625g` | COQ4 C1 decarboxylation | Cross-species experiments support the step, but oxidative versus sequential decarboxylation/hydroxylation remains unresolved. |
| `R19` | empty | unresolved COQ6-associated monooxygenation | No exact GPR is asserted because regioselectivity, product redox state, and ferredoxin/reductase coupling do not yet match the encoded reaction. |

An empty `R19` GPR means that the catalyst identity is unknown; it does not
mean that the reaction is spontaneous. The COQ6, COQ8, and COQ9 candidates
remain in the model as gene objects. COQ6 is not assigned to the unresolved
`R19`; COQ8 and COQ9 are no longer represented as direct atom-changing
catalysts. Their possible accessory or synthome-coupling roles remain
unmodeled. Fitness phenotypes were used as a consistency check, not as proof
of a reaction assignment.

This follow-up changes only the six GPRs and their evidence notes. It does not
add a CoQ9 demand, open a bound, alter reaction chemistry, or activate the
pathway.

## Why the removed reactions are not a valid Yarrowia pathway

| Reaction | Reviewed problem | Steady-state consequence |
|---|---|---|
| `R189` | The equation is a chemically plausible Q9-like 4-hydroxybenzoate prenylation skeleton, but its name/GPR describe a CAAX protein farnesyltransferase, it is placed in cytosol, and both branch-specific endpoints are disconnected. A future mitochondrial CoQ9 reconstruction may reuse this chemistry, but not this reaction as currently encoded. | `v_R189 = 0` |
| `R2250` | An octaprenyl-labelled substrate is joined to a chemically hexaprenyl product. Its GPR consists of protein prenyltransferase/dolichol-pathway proteins rather than COQ2. | Its substrate has no producer, so `v_R2250 = 0`. |
| `R2247`-`R2249` | These reactions form a cytosolic-to-nuclear continuation of the invalid `R2250` branch. They are not an independently reachable route. | Mass balance propagates `v = 0` downstream. |
| `R2242` | The reaction is a bacterial-style 2-polyprenyl-6-hydroxyphenol methylation in the nucleus, but its GPR is `YALI1E01159g`/DIM1, the 18S-rRNA dimethyltransferase EC 2.1.1.183. Eukaryotic COQ3 uses different carboxylated intermediates in mitochondria. | Its product has no consumer, so `v_R2242 = 0`. |

`R2250` and the retained mitochondrial `R407` share the same MetaNetX reaction
mapping. Merely moving `R2250` into mitochondria would therefore create a
duplicate rather than repair the model.

## Gene identity audit

| YALI1 gene | Established symbol | Protein function | Evidence status | Model decision |
|---|---|---|---|---|
| `YALI1B21088g` | no established Yarrowia symbol | Protein farnesyltransferase/geranylgeranyltransferase type-I alpha subunit | computational/homology annotation | Remove from quinone GPRs. |
| `YALI1D17983g` | no established Yarrowia symbol | Protein farnesyltransferase beta subunit | computational/homology annotation | Remove from quinone GPRs. |
| `YALI1E11415g` | no established Yarrowia symbol | Geranylgeranyltransferase type-I beta subunit | computational/homology annotation | Remove with `R2250`. |
| `YALI1E16694g` | no established Yarrowia symbol | Rab geranylgeranyltransferase alpha subunit | computational/homology annotation | Remove with `R2250`. |
| `YALI1E33302g` | no established Yarrowia symbol | Rab geranylgeranyltransferase beta subunit | computational/homology annotation | Remove with `R2250`. |
| `YALI1E01159g` | `DIM1` | 18S-rRNA dimethyltransferase and ribosome-maturation factor | curated annotation | Remove from the quinone methylase reaction. |
| `YALI1C26017g` | no established Yarrowia symbol (`COQ1` candidate) | Long-chain trans-prenyl diphosphate synthase assigned to the CoQ side-chain step | heterologous functional support for Yarrowia CoQ9 synthesis; native mitochondrial localization not directly verified | Keep on `R763`; retain the four-IPP lump as a model convention. |
| `YALI1F08349g` | `COQ2` | Mitochondrial 4-hydroxybenzoate polyprenyltransferase | homology-supported curated annotation; no direct Yarrowia locus experiment found | Keep on `R407`. |
| `YALI1B20835g` | `COQ3` candidate | CoQ O-methyltransferase | Cross-species biochemical/structural evidence and compatible predicted fold; native Yarrowia biochemistry unverified | Assign to `R715` and `R385`. |
| `YALI1C25352g` | `COQ5` candidate | CoQ-ring C-methyltransferase | Cross-species biochemical/structural evidence and compatible predicted fold; native Yarrowia biochemistry unverified | Assign to `R18`. |
| `YALI1E18269g` | `COQ7` candidate | Demethoxyubiquinone hydroxylase | Cross-species biochemical/structural evidence and compatible predicted fold; native Yarrowia biochemistry unverified | Assign to `R695`. |
| `YALI1F34625g` | `COQ4` | CoQ-ring C1 decarboxylase and synthome-organising protein | Reviewed homology annotation plus cross-species experiments; exact product mechanism unresolved | Assign provisionally to `R40`. |
| `YALI1A08781g` | no established Yarrowia symbol (`COQ6` candidate) | FAD-dependent CoQ monooxygenase candidate | Family/fold compatible, but the exact `R19` regioselectivity, redox state, and electron partners conflict with available mechanisms | Do not assign to `R19`; leave its GPR unresolved. |
| `YALI1B20527g` | no established Yarrowia symbol (`COQ8` candidate) | ADCK-family ATPase/kinase-like accessory protein | Cross-species evidence supports pathway regulation and complex stability, not direct ring chemistry | Remove from repeated direct-catalyst `AND`; retain the gene object. |
| `YALI1F34675g` | no established Yarrowia symbol (`COQ9` candidate) | Lipid-binding, substrate-presenting accessory protein | Cross-species evidence supports COQ7 cooperation, not direct atom transfer | Remove from repeated direct-catalyst `AND`; retain the gene object. |

Cross-assembly mappings used here are
`YALI1F08349g <-> YALI0F05610g <-> YALI2_F00880g` and
`YALI1B20835g <-> YALI0B15884g <-> YALI2_C00448g`.

## Evidence

- [UniProt Q6C2S2 and A0A1H6PM88 (COQ2 cross-assembly records)](https://www.uniprot.org/uniprotkb/Q6C2S2/entry)
- [UniProt Q6CEG2 and A0A1D8N802 (COQ3 cross-assembly records)](https://www.uniprot.org/uniprotkb/Q6CEG2/entry)
- [UniProt Q6C7H6 (YALI0E00770g/DIM1)](https://www.uniprot.org/uniprotkb/Q6C7H6/entry)
- [ENZYME EC 2.1.1.183: 18S-rRNA dimethyltransferase](https://enzyme.expasy.org/EC/2.1.1.183)
- [ENZYME EC 2.1.1.222: bacterial 2-polyprenyl-6-hydroxyphenol methylase](https://enzyme.expasy.org/EC/2.1.1.222)
- [Rhea 44504: 4-hydroxybenzoate polyprenyltransferase](https://www.rhea-db.org/rhea/44504)
- [Yarrowia mitochondrial complex-I work reporting endogenous Q9](https://doi.org/10.1016/S0005-2728(02)00307-9)
- [Yarrowia CoQ1 functional study](https://journals.asm.org/doi/10.1128/mbio.00342-24)
- [KEGG R08781: balanced oxidized terminal methylation convention](https://www.kegg.jp/entry/R08781)
- [UbiG/COQ3 structure 4KDC](https://www.rcsb.org/structure/4KDC)
- [Saccharomyces Coq5 structure 4OBW](https://www.rcsb.org/structure/4OBW)
- [Human COQ7:COQ9 structure 7SSS](https://www.rcsb.org/structure/7SSS)
- [Pelosi et al.: Coq4 oxidative decarboxylation](https://doi.org/10.1016/j.molcel.2024.01.003)
- [Nicoll et al.: reconstructed eukaryotic CoQ synthome steps](https://doi.org/10.1038/s41929-023-01087-z)

## Explicitly not fixed here

- The exact `R19` catalyst remains unresolved; an empty GPR records this
  uncertainty and does not assert a spontaneous reaction.
- COQ8/COQ9 accessory dependence and other synthome-level coupling remain
  unmodeled even though their unsupported direct-catalyst assignments were
  removed.
- The model still lacks a validated net CoQ pool dilution/demand term.
- Consequently the corrected CoQ9 synthesis route remains structurally blocked
  (`FVA = [0, 0]`) rather than being activated by an arbitrary sink.
- `R1889`, `R305`, and `R570` retain pre-existing proton/charge residuals in the
  wider Q/QH2 respiratory component; this patch validates the synthesis route,
  not the entire respiratory connected component.
- `R2243`-`R2246` remain as a separate inert nuclear SAM/SAH transport artifact;
  they were not folded into this narrowly scoped cleanup.

Those are the separate `18-Ubiquinone and other terpenoid-quinone biosynthesis`
curation problem and must not be inferred solved by this branch cleanup.

## Validation

- Final canonical model SHA-256:
  `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`.
- A separately generated candidate established the same CoQ9 stoichiometry;
  the final artifact additionally carries the canonical B-group metadata and
  the independently audited `R385` reaction-level EC correction.
- All 11 CoQ9 synthesis-route reactions pass element and charge balance.
- Model counts remain 2,313 reactions, 1,877 metabolites, and 1,074 genes. The
  chemistry follow-up changes only `R763` and `R385` stoichiometry; the later
  GPR follow-up changes only `R715`, `R385`, `R18`, `R695`, `R40`, and `R19`
  GPRs and evidence notes.
- Default growth changed only by floating-point noise:
  `1.3875914754870884 -> 1.387591475487086`.
- CoQ demand reactions remain absent and the only demand is the pre-existing
  `R1373`; the model contains no sink reactions.
- All six removed reactions had FVA range `[0, 0]` before deletion.
- Default growth was `1.3875914754870913` before and
  `1.3875914754870884` after deletion (floating-point noise only).
- The rebuilt model contains 20 independent B-group AA-tRNA biomass reactions.
- The five affected experimental-essential loci remain inside the screen
  denominator; this cleanup does not receive credit for solving their biology.
- PO1f SD-Leu positive-only TP/FN at 1%, 5%, 10%, and 15% are respectively
  `57/265`, `63/259`, `67/255`, and `79/243`; the TP series is unchanged.
- The formal cleanup accepts only the exact reviewed GPR, stoichiometry, bounds,
  reversibility, and compartment signatures. The raw-source and canonical R189
  proton variants are enumerated separately; any other future change fails
  closed rather than being silently deleted.
