# Muon-pair DeltaR regression tests

From the repository root, with the analysis dependencies installed:

```sh
python -m unittest discover -s tests -p 'test_muon_muon_dr.py' -v
```

The tests use in-memory Awkward/Coffea vectors and production histogram
booking/filling. They cover DeltaR with nonzero DeltaEta and wrapped phi,
first-two input ordering, multiplicity masks and empty inputs, signed-weight
bin contents and variances, named-axis projection, and the preserved 50-bin
range `[0, 6.28]` with flow bins. They also resolve and book `ttjets_base`
through `SidmProcessor` to check the renamed key and axis.

No event files, notebook execution, or processor event loop are used. These
tests establish synthetic software behavior, not physics validation or
compatibility of previously saved outputs.
