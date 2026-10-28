# Synthetic angular regression tests

For usage, see [the angular helper documentation](../sidm/tools/README_delta_phi.md).

Run the synthetic regression tests from the repository root with the analysis
dependencies installed:

```sh
python -m unittest discover -s tests -p 'test_delta_phi.py' -v
```

The tests create their inputs in memory; no event files are needed.
