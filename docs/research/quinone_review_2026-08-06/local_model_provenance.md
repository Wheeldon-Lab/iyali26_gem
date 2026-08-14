# Local model provenance audit

Date: 2026-08-06  
Scope: read-only comparison of locally available Yarrowia GEM files  
Evidence class: model provenance only; not biological validation

## Models inspected

| Label | Local file identity | Size | Quinone observation |
|---|---|---:|---|
| iYali4 local copy | SHA-256 `c0b54165301bfba7edc7083efee1106309f049494f3a3f35fb436306d4fc86ae` | 1924 reactions | Hexaprenyl/CoQ6 route and repeated seven-gene `AND` GPR |
| Duplicate `iYali.xml` local copy | same SHA as the inspected iYali4 XML | 1924 reactions | Byte-identical to the inspected iYali4 XML |
| iYali5 | SHA-256 `b73e99d88d2de954a0abe49bd326e8c1acf11c2c579c5a4c5629c78cbf5f0019` | 1986 reactions | Retains the same hexaprenyl chemistry and the same seven YALI0 loci on five steps |
| iYli21 | SHA-256 `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631` | 2285 reactions | Retains the same route after YALI0→YALI1 ID conversion |
| Current iYali26 | SHA-256 `0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415` | 2313 reactions | Retains the iYli21 main route; the separately reviewed spurious branches were removed |

The local iYali5 README explicitly identifies iYali4/iYali4_corr as its model
source. Its citation field is a template, so it is used only as local lineage
metadata. Independent publication provenance is recorded separately as
`SRC-PR-001` (Kerkhoven et al. 2016, iYali4) and `SRC-PR-002` (Guo et al. 2022,
iYli21); neither publication is treated as independent validation of the copied
quinone chemistry or GPR.

## Repeated GPR lineage

The same seven-locus conjunction appears on the five downstream reactions in
iYali4 and iYali5:

`YALI0F27247g and YALI0B15664g and YALI0A09042g and YALI0F27313g and YALI0C18205g and YALI0B15884g and YALI0E15224g`

iYli21 converts it to:

`YALI1F34625g and YALI1B20527g and YALI1A08781g and YALI1F34675g and YALI1C25352g and YALI1B20835g and YALI1E18269g`

The reaction correspondences visible in the local files are:

| iYali4/iYali5 role | iYli21/iYali26 reaction | GPR behavior |
|---|---|---|
| 2-hexaprenyl-6-methoxy-1,4-benzoquinone methyltransferase | `R18` | same seven-locus `AND` after identifier conversion |
| 2-hexaprenyl-6-methoxyphenol monooxygenase | `R19` | same |
| terminal methyltransferase | `R385` | same |
| quinone oxidoreductase/hydroxylation step | `R695` | same |
| early O-methyltransferase | `R715` | same |

`YALI0F05610g`, mapped in the current model to `YALI1F08349g`, is separately
assigned to the hydroxybenzoate polyprenyltransferase reaction (`R407` in
iYli21/iYali26). The local files alone do not establish the correctness of that
identity or the substrate chain length.

## Provenance conclusion

Agreement among these four GEM generations is not four independent pieces of
evidence. The exact chemistry and repeated GPR are visibly inherited through a
model lineage. The completed source review evaluates each locus and reaction
separately in `evidence_ledger.tsv`; regardless of those biological verdicts,
the four-model agreement remains one copied modeling assertion with several
descendants.

This audit does not show that the inherited rule is false. It shows that model
replication cannot be counted as experimental replication.
