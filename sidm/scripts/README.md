# Read the configured cuts

`sidm/scripts/print_cuts.py` resolves YAML anchors/merges and flattens inherited
cut lists, preserving their order and repeated cuts. It requires only PyYAML
(`python3 -m pip install PyYAML`), with no Coffea installation or event inputs.
From the repository root:

```bash
python3 sidm/scripts/print_cuts.py                         # all selections
python3 sidm/scripts/print_cuts.py --list                  # selection names
python3 sidm/scripts/print_cuts.py -s 4mu_sr               # one selection
python3 sidm/scripts/print_cuts.py -s '*_sr'               # matching selections
python3 sidm/scripts/print_cuts.py -s 4mu -s 2mu2e         # several selections
python3 sidm/scripts/print_cuts.py > cuts.txt              # save the full report
python3 sidm/scripts/print_cuts.py path/to/selections.yaml --json
python3 sidm/scripts/print_cuts.py --compare base_sr 4mu_sr # compare two selections
python3 sidm/scripts/print_cuts.py -s 4mu_sr --color never  # disable colors
python3 sidm/scripts/print_cuts.py --compare base_sr 4mu_sr --color always | less -R
```

`--compare LEFT RIGHT` accepts two exact selection names, for example denominator
then numerator. It shows shared cuts (`=`), left removals (`-`, red), and right
additions (`+`, green), grouped by section and object collection. Each row shows
the original `[left:right]` cut positions. Separate summaries check equal cut
labels (including repeats) and equal order within each group; reordered cuts
appear as removal/addition pairs. Missing and empty groups both mean no cuts.
This compares configuration labels, not equivalent function behavior or a proof
that one event selection is a subset of another. Add `--json` for a structured
comparison report. Comparison cannot be combined with `-s` or `--list`.

Colors are automatic in terminals and disabled for redirected/piped output,
`NO_COLOR`, or `TERM=dumb`. Use `--color always` to force them or `--color never`
to disable them; explicit modes override auto detection. JSON never has colors.

The default config is located relative to the script, so the script also works
from another directory. Output separates object, LJ/post-LJ, and event cuts.
These are the configured cut labels; their Python implementations live in
`sidm/definitions/cuts.py`. The script does not evaluate cut functions or prove
that a runtime cut succeeded. YAMLs for sample lists, cross sections, and
histograms do not contain selection cuts and are reported as such.

Run the feature tests with only PyYAML installed:

```bash
python3 -m unittest discover -s tests -p test_print_cuts.py
```

