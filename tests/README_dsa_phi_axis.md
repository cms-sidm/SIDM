# Back-to-back DSA phi-axis tests

From the repository root, with the analysis dependencies installed:

```sh
python -m unittest discover -s tests -p 'test_dsa_phi_axis.py' -v
python -m unittest discover -s tests -p 'test_delta_phi.py' -v
```

The focused tests exercise actual weighted filling, named projection and slicing,
bin contents and variances across the axis rename, exact-pi inclusion, existing
3D pair selection, empty collections and preservation of the real eta histogram.
The angular suite checks the surrounding wrapped-angle behavior. All inputs are
synthetic in-memory arrays; no datasets or analysis notebooks are used.
