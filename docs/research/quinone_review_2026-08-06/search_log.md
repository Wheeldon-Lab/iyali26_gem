# Reproducible search log

Date: 2026-08-06  
Scope: root-agent searches supplementing the three independent evidence lines  
Mode: read-only; search hits are candidates until claim-level verification

## Search families

The searches were framed around competing explanations rather than around a
desired model edit. Query families included:

1. `Yarrowia lipolytica ubiquinone CoQ9 COQ1 chain length`
2. `Yarrowia lipolytica Q9 complex I HPLC`
3. `Yarrowia lipolytica coenzyme Q9 content dry weight HPLC`
4. `Yarrowia lipolytica COQ2 COQ3 COQ4 COQ5 COQ6 COQ7 COQ8 COQ9`
5. `Yarrowia lipolytica ubiquinone biosynthesis mitochondria localization`
6. `Yarrowia lipolytica coenzyme Q turnover dilution biomass model`
7. `Yarrowia lipolytica quinone genome scale model CoQ6 hexaprenyl`
8. `Yarrowia lipolytica CoQ synthome ER mitochondria contact`

For each potentially material hit, the review attempts to trace the exact
primary article, stable DOI/PMID/PMCID, strain, growth condition, analytical
method, and an inspectable passage or figure. Reviews and databases are used
to locate or cross-check primary sources, not to replace them.

## Root-verified candidate sources

The `SRC-R*` labels below were temporary root-search IDs. Their final inventory
crosswalk is: `SRC-R001 → SRC-QC-003`, `SRC-R002 → SRC-QC-002`,
`SRC-R003 → SRC-QC-004`, `SRC-R004 → SRC-DB-001`, and
`SRC-R005 → SRC-QC-008`. Only the final IDs are used by the evidence ledger.

### SRC-R001 — Saeed et al. (2024), mBio

- DOI: `10.1128/mbio.00342-24`
- PMID: `38747615`
- PMCID: `PMC11237637`
- Stable full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11237637/>
- Search role: direct test of the product chain length encoded by the
  *Y. lipolytica* `COQ1` coding region.
- Verified observation: a strain of *S. cerevisiae* whose `COQ1` locus was
  replaced with the *Y. lipolytica* coding sequence produced CoQ9 by HPLC.
- Material limitation: the construct used the *S. cerevisiae* mitochondrial
  import signal. It therefore supports a chain-length function of the coding
  sequence but does not directly validate the native *Yarrowia* targeting
  peptide or localization.
- Additional limitation: this is heterologous complementation, not a native
  *Yarrowia* locus perturbation.
- Candidate current-model mapping:
  `YALI1C26017g — COQ1 — solanesyl/nonaprenyl-diphosphate synthase candidate`
  (direct heterologous functional support for a CoQ9-producing side chain;
  native localization remains unverified).

### SRC-R002 — Dröse, Zwicker & Brandt (2002), BBA Bioenergetics

- DOI: `10.1016/S0005-2728(02)00307-9`
- Stable landing page:
  <https://www.sciencedirect.com/science/article/pii/S0005272802003079>
- Search role: endogenous quinone associated with purified *Y. lipolytica*
  mitochondrial complex I.
- Verified context: strain PIPO; modified YPD with 2.5% glucose; 27°C;
  purified mitochondrial complex I; quinone measured by HPLC.
- Verified observation: the extracted quinone was identified as Q9. Two
  preparations contained approximately 0.2 and 0.4 mol Q9 per mol complex I
  after the reported purification/proteolysis workflow.
- Material limitation: this is not a whole-cell CoQ abundance measurement and
  cannot directly supply a biomass or dilution coefficient.
- Audit resolution: the older direct classification source was traced to
  Yamada et al. 1976 (`SRC-QC-001`); Table 2 identifies *E. lipolytica*
  `CBS 6124` as a Q9 system. This still does not establish W29/PO1f identity or
  a target-condition abundance.

### SRC-R003 — US20090142322A1 / US8815567B2

- Stable text: <https://patents.google.com/patent/US20090142322A1/en>
- Evidence class: grey literature/patent; not peer reviewed.
- Search role: possible quantitative lower-tier evidence.
- Verified context: wild-type *Y. lipolytica* ATCC 20362; extracted oil;
  HPLC-DAD with a CoQ9 standard.
- Reported observation: CoQ9 in extracted oil was approximately 0.2–0.3%.
- Material limitations: the denominator is extracted oil, not dry cell weight;
  industrial conflicts and incomplete uncertainty reporting apply. A separate
  statement that this strain may contain roughly 40% lipid per dry cell weight
  cannot be multiplied into a biomass coefficient without proving that the
  values refer to matched samples and recovery conditions.
- Patent-family correction: `US20090142322A1` and its grant
  `US8815567B2` are one family. `WO2008073367A1` is a different carotenoid
  patent family that is cited only as engineering-background provenance; it is
  neither a family member nor an independent replication of the ATCC 20362
  CoQ9-in-oil measurement.

### SRC-R004 — KEGG pathway `yli00130`

- Stable entry: <https://www.kegg.jp/entry/yli00130>
- Evidence class: database annotation; identity aid only.
- Search role: candidate locus mapping for the COQ pathway.
- Candidate annotations observed include COQ2, COQ3, COQ4, COQ5, COQ6 and a
  COQ7-like hydroxylase. These assignments require source-level and
  protein-family auditing before they can establish a GPR.

### SRC-R005 — conserved CoQ biosynthetic domains at ER–mitochondria contacts

- PMCID: `PMC6446851`
- Stable full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6446851/>
- Evidence class: primary work in *S. cerevisiae*, not *Yarrowia*.
- Search role: conserved-mechanism context for a multi-protein CoQ synthome.
- Material limitation: it may motivate hypotheses about organization or
  localization, but it cannot establish a *Yarrowia* locus-specific GPR.

## Explicitly unresolved during root search

- Direct native localization of *Y. lipolytica* COQ1/Coq1.
- A native *Yarrowia* gene perturbation that measures CoQ9 production for each
  current candidate COQ locus.
- A condition-matched whole-cell CoQ9 pool in mol or mmol per gDW for W29 or
  PO1f.
- CoQ9 turnover or growth-dilution measurements suitable for a numerical GEM
  coupling coefficient.
- Experimental evidence that the seven proteins in the inherited current GPR
  are jointly required for every one of five distinct catalytic steps.

These are recorded as unresolved search outcomes, not as evidence that the
corresponding biology is absent.

## Independent-audit correction log

- The root search could not reliably recover the historical isolate from the
  OCR text of Yamada et al. (1976). The independent auditor subsequently
  opened Table 2 and read the entry as *E. lipolytica* `CBS 6124`, assigned to
  Q9. The audited table reading supersedes the earlier unresolved note; any
  relationship to modern W29/PO1f backgrounds still requires a separate
  strain-provenance check and is not assumed here.
