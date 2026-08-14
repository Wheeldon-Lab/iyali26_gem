# Reproducibility verification

- Status: **VERIFIED — exact replay**
- Verification date: 2026-08-06
- Primary run: `docs/experiments/quinone_runtime_counterfactual_20260806_final`
- Independent replay: `/tmp/quinone_runtime_counterfactual_final_replay_20260806`
- Comparison: all 11 original run artifacts matched byte-for-byte (`diff -qr` returned no differences)
- Primary manifest SHA-256: `f14ec84357ebe1ea04f3a1892ed80124aa14460d8c9f7d638af104ae5d3851b5`
- Replay manifest SHA-256: `f14ec84357ebe1ea04f3a1892ed80124aa14460d8c9f7d638af104ae5d3851b5`
- Experiment script SHA-256: `ad43f632329b365b967970b6adbd10b251bce28e976c607a6875b1c2b0638003`
- Experiment-design SHA-256: `5c13517f6e08abacdae3a14ef7e2aebca6920976ed409933b442e55dbeda09b9`

## Independent scientific audit

The final wording incorporates the independent audit boundaries:

- maximum-demand values are zero-growth stoichiometric reachability maxima, not physiological capacities;
- oxidized and reduced terminal redox representations are both database-compatible, so the native *Yarrowia* redox microspecies remains unresolved;
- the GPR contrast is route-wide and includes candidate assignments for previously ungated `R39/R40`;
- none of these qualifications changes the immediate modeled diagnosis of missing terminal net demand.

## No-write audit

The following protected inputs retained their frozen hashes after the primary run and replay:

- `model.xml`: `0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415`
- `data/iyali26.xml`: `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d`
- Obsidian `17-Quinone Biosynthesis.md`: `4dfa2e30190fbd835b7520146c76eca60fc4dbadc2078d6a8540bbb5416baad3`

This post-run attestation is intentionally outside the original manifest that it verifies.
