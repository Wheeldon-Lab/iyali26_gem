# Independent source-audit protocol

Date: 2026-08-06  
Status: fixed before the independent audit

## Unit of audit

The audit unit is an atomic, material claim. A claim is material when changing
its verdict could change any of the following conclusions:

- native CoQ chain length or chemical species;
- current reaction chemistry or compartment interpretation;
- a gene identity, catalytic role, or proposed GPR decomposition;
- existence or numerical magnitude of a net CoQ requirement;
- whether a source can justify a model representation rather than only a
  hypothesis or feasibility experiment.

## Required checks

For every audited claim, the auditor must independently open the cited source
where access permits and check:

1. source identity (title, authors, year, DOI/PMID/URL);
2. exact locator and whether the cited passage/figure supports the claim;
3. organism and strain;
4. medium, carbon source, growth phase and sample type when relevant;
5. measurement or perturbation method;
6. directness and whether the wording overgeneralizes the experiment;
7. relevant counterevidence or alternative interpretation;
8. independence from copied database or GEM assertions.

Permitted verdicts are `supported`, `partially_supported`, `unsupported`,
`contradicted`, and `unverified`. `Unverified` is used for inaccessible or
insufficiently identifiable evidence; it is not treated as support.

## Audit coverage

Two coverage measures are reported:

\[
C_{claims}=\frac{N_{material\ claims\ independently\ audited}}
{N_{material\ claims\ in\ the\ synthesis}}
\]

and

\[
C_{sources}=\frac{N_{high\text{-}impact\ primary\ sources\ independently\ opened}}
{N_{high\text{-}impact\ primary\ sources\ cited}}
\]

A source is high impact if it supplies the only or strongest evidence for a
material claim, a numerical value, a gene-function assignment, or a claimed
conflict. Database records and inherited GEM files are audited for provenance
but are not counted as independent primary studies.

## Independence rules

- The auditor receives the neutral question, ledger, and sources but no target
  verdict or desired edit.
- Agreement between agents is not counted as evidence.
- Multiple GEM generations with copied chemistry/GPRs count as one provenance
  lineage.
- A review that cites a primary study does not count as an independent
  replication of that study.
- “No evidence found” and “evidence of absence” remain separate outcomes.

## Human gate

The audit may classify evidence and recommend what is or is not presently
defensible. It cannot authorize a chemistry, GPR, bound, demand, biomass,
pipeline, model, or Obsidian change.
