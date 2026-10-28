# Back-to-back DSA phi histogram

The histogram key remains `back_to_back_dsa_pair_phi_diff`. Its angular axis is
named `back to back dsa_pair_phi_diff`, so a booked or loaded histogram made
with the corrected definition can be projected by name:

```python
phi_projection = phi_hist.project("back to back dsa_pair_phi_diff")
```

This name correction preserves the display label, wrapped absolute delta-phi
values, 100-bin range including exact pi, back-to-back pair selection and weights.

Previously saved phi histograms retain the incorrect axis name
`back to back dsa_pair_eta_diff`; loading them does not rename it. Check
`phi_hist.axes.name` and use the matching name for projections or named slicing.
Do not combine old and new outputs until their axis names have been reconciled
explicitly. Matching names alone does not make older outputs with different
angular definitions or binning compatible. No automatic conversion is performed.

The separate `back_to_back_dsa_pair_eta_diff` histogram still uses
`back to back dsa_pair_eta_diff`. Do not rename that genuine eta axis.
