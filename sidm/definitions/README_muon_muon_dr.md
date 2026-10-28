# Muon-pair separation

`muon_muon_dr` measures $\Delta R = \sqrt{(\Delta\eta)^2 + (\Delta\phi)^2}$
for the first two muons in the supplied `muons` collection, using the vector
method's wrapped azimuthal separation. Events with fewer than two muons do not
contribute. The processor orders muons by descending $p_T$ before selections;
the histogram itself does not reorder the supplied collection. It keeps 50 bins from 0 to 6.28, with weighted storage
and underflow/overflow bins; values at or above 6.28 enter overflow.

The `ttjets_base` histogram collection includes it. Custom histogram collections
should request `muon_muon_dr`. Its output key and observable axis both use that
name; for example, with a processor result dictionary `output`:

```python
histogram = output[sample]["hists"]["muon_muon_dr"]
separation = histogram.project("muon_muon_dr")
```

Projection sums any other axes, including channels. Slice those axes first when
a separate channel or reconstruction choice is needed.

## Existing outputs

The former key and axis name `muon_muon_dphi` already held **DeltaR**, despite its
name. The calculation, label, muon ordering, event mask, binning and weights are
unchanged. Newly produced outputs use `muon_muon_dr`; the old name is not an alias.
Update custom configurations and scripts that request the old key or axis.

Saved outputs retain their original `muon_muon_dphi` key and axis. Read those
outputs with their stored names. Changing a dictionary key alone does not rename
the axis. Do not blindly merge old and new output dictionaries: the different
keys can remain separate, and histograms with different axis names are not
interchangeable. Regenerate outputs consistently, or explicitly align both names
in a separate conversion after checking binning, selections and normalization.
No automatic output conversion is provided.
