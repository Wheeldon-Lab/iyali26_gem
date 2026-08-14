# UMP microspecies component review

The pH 7.3 Rhea identity for UMP with charge -2 is `C9H11N2O9P`. The current
model carries five UMP compartment copies as a single deferred component:
`m1002[C_nu]`, `m1120[C_mi]`, `m149[C_cy]`, `m1967[C_er]`, and `m838[C_go]`.

Activating this family as one change repairs the hydrogen residual in `R612`,
but a prospective balance audit regresses previously balanced reactions `R58`,
`R781`, and `R782` (and leaves `R804` unresolved). Therefore this implementation
adds the supported `R612 = YALI1E31685g` GPR but does not activate the UMP formula
change. The existing curation row remains `component_review` until the connected
pyrimidine/nucleotide chemistry component can be repaired and tested together.
