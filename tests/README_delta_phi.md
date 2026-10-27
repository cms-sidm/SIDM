# Synthetic angular regression tests

With the repository's analysis dependencies installed, run from the repository root:

```sh
python -m unittest discover -s tests -p 'test_delta_phi.py' -v
```

These tests build all vectors and matches in memory. They exercise the real
histogram definitions, weighted histogram fills, vector behavior and DSA
cross-cleaning predicates without opening event files or running a processor.

## Angular contract and audit

`sidm.tools.utilities.delta_phi(phi1, phi2)` returns the signed difference in
`[-pi, pi)` radians; its absolute value is in `[0, pi]`. It calculates in float64,
including for float32 inputs, so the rounded float32 representation of pi cannot
put an antipodal pair above the histogram boundary. The corrected full-range
axes end one float64 step above pi so exact pi falls in the last regular bin.

The following definitions now use the wrapped absolute separation:

- `lj_lj_absdphi` and its `invmass`, `absdR`, `absdeta` and `ptRatio` correlations.
- `genE_genE_dPhi` and `genMu_genMu_dPhi`; their existing 0–1 radian zoom ranges
  remain, so larger separations still intentionally enter overflow.
- `dsa_pair_phi_diff`, `back_to_back_dsa_pair_phi_diff` and
  `parallel_dsa_pair_phi_diff`. Pair selection still uses the existing 3D angle.

`dR_outer` also wraps its azimuthal component. This affects `segmatch veto + dR
outer`, `old_all`, `all` and `all + charge` in the DSA muon cuts. Synthetic tests
check that a segment-matched pair across the branch cut is recognized as close,
and preserve the segment-count, charge and missing-match exceptions.

The LJ and PF-muon event cuts and `genAs_genAs_absdphi` already use the vector
library's wrapped `delta_phi`; their definitions remain unchanged. The cosmic
distance expressions also remain unchanged: for canonical phi coordinates,
`(pi - abs(phi1 - phi2))**2` is invariant under ordinary angular wrapping.
Individual phi-coordinate plots retain their coordinate definitions.

`muon_muon_dphi` has an ambiguous name: both its current calculation and label
specify Delta R. This patch preserves that behavior; deciding whether to rename
it or change its observable requires clarification of the plot's intent. The
existing axis identifiers are preserved, including the misnamed
`back to back dsa_pair_eta_diff` axis in the corresponding phi histogram.

## Compatibility and limits

Full-range angular axes change from 0–2pi (LJ) or 0–7 (DSA) to 0–pi. Do not
combine old and new histogram outputs as if their axes or definitions matched.
The outer-track fix can also change selected candidates and downstream yields;
rebinned old histograms cannot establish those effects.

The suite covers float32/float64 endpoints, periodicity, pair swapping,
jagged/empty/missing inputs, leading-pair choice, histogram/cut agreement,
weighted integrals and variances, and the affected cross-cleaning predicates.
It establishes these code behaviors, not their impact on event samples or
scientific validation of a selection. Candidate-collection conventions and cut
thresholds are unchanged.
