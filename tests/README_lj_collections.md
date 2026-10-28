# Synthetic lepton-jet collection tests

For usage, see [the collection documentation](../sidm/definitions/README.md).

Run the tests from the repository root in the SIDM analysis Python environment,
with the [analysis dependencies](../README.md#3-install-requirements) installed:

```sh
PYTHONPATH=. python tests/test_lj_collections.py
```

The tests create Awkward records in memory and use the configured selections;
no event files are needed. They check pair selection and ordering, stable
identity, parent/child synchronization, channel exclusivity, empty collections,
and agreement between pair cuts and histogram values.
