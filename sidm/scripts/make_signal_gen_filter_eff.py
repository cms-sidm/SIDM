#!/usr/bin/env python3
"""Build sidm/configs/signal_gen_filter_eff.yaml, the per-point generator-filter
efficiency of the 2018 v10 signal samples.

Inputs (repo files):
  - sidm/studies/truth_kinematics_forAN/central_genFilterEfficiencies.yml:
    efficiency of the central production's filter (pT > 5 GeV), per point.
  - sidm/studies/truth_kinematics_forAN/v10_filter_repass.json: counts of stored
    v10 events passing the 5 GeV cut string (genFilterEmulation channel).
  - sidm/configs/ntuples/signal_{2mu2e,4mu}_v10.yaml: the sample list.

The v10 samples ran the same filter at pT > 4 GeV, so
eff_v10 = eff_central / repass_5gev. Central keys are matched to our sample
names by the (channel, mXX, mA, ctau) values, not by string. The script writes
nothing if any v10 sample lacks a central efficiency or re-pass count.

Run from the repo root:
    python sidm/scripts/make_signal_gen_filter_eff.py
"""

import argparse
import hashlib
import json
import math
import os
import re
import sys

import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
STUDY = "sidm/studies/truth_kinematics_forAN"
CENTRAL = f"{STUDY}/central_genFilterEfficiencies.yml"
REPASS = f"{STUDY}/v10_filter_repass.json"
NTUPLES = ["sidm/configs/ntuples/signal_2mu2e_v10.yaml", "sidm/configs/ntuples/signal_4mu_v10.yaml"]
OUTPUT = "sidm/configs/signal_gen_filter_eff.yaml"
REPASS_CHANNEL = "genFilterEmulation"

CENTRAL_RE = re.compile(r"^XXTo2ATo(2Mu2E|4Mu)_mXX-([0-9p]+)_mA-([0-9p]+)_lxy-[0-9p]+_ctau-([0-9p]+)$")
SAMPLE_RE = re.compile(r"^(2Mu2E|4Mu)_([0-9p]+)GeV_([0-9p]+)GeV_([0-9p]+)mm$")

HEADER = """\
# Generator-filter efficiency of the 2018 signal productions, per sample.
# The v10 samples were produced with a generator filter (>= 4 prompt final-state
# e/mu with pT > 4 GeV, |eta| < 2.4, production vertex rho < 740 cm, |z| < 960 cm;
# cut string "... && pt>4. && isLastCopy() && isPromptFinalState() &&
# fromHardProcessFinalState()" in the GEN process PSet stored in the v10 files, read
# back in sidm/studies/run3_signal_validation/gen_requirement.py on branch
# run3-signal-validation). Failing events are not stored, so the normalization needs
# eff_filter: utilities.get_lumixs_weight applies it, choosing the field by the
# production tag make_fileset derives from the location config (v10 ntuples -> v10,
# ffntuples.yaml -> central, anything else -> no factor, with a warning).
# Fields:
#   central:     {eff, err} of the central production's filter, per point, from
#                phylsix/Firefighter ffConfig/python/production/Autumn18/sigmc/central/,
#                copied in sidm/studies/truth_kinematics_forAN/central_genFilterEfficiencies.yml
#   repass_5gev: {n_all, n_pass, value, err} stored v10 events before/after the
#                central (5 GeV) cut string, the genFilterEmulation channel, from
#                sidm/studies/truth_kinematics_forAN/v10_filter_repass.json; binomial err,
#                which is 0 where value = 1
#   v10:         {eff, err} = central.eff / repass_5gev.value, relative errors in quadrature
# The central production (McM EXO-RunIIFall18wmLHEGS-02157) ran the same filter string
# at pT > 5 GeV. Every event passing the 5 GeV string passes the 4 GeV one, so
# eff_v10 = eff_central / repass_5gev; the re-pass deficit is the 4-5 GeV band.
# The re-pass uses raw event counts while the normalization divides by the sum of
# generator weights; the two agree because the v10 signal genWeight is constant (1.0 for
# every event in the checked files, with Runs genEventSumw = genEventCount).
# Values rounded to 4 significant figures.
"""


def _num(token):
    """'1p2' -> 1.2"""
    return float(token.replace("p", "."))


def _sig4(x):
    return float(f"{x:.4g}")


def central_by_point(path):
    """Map (channel, mXX, mA, ctau) -> (eff, err) from the Firefighter table."""
    with open(path, encoding="utf8") as f:
        raw = yaml.safe_load(f)
    table = {}
    for key, val in raw.items():
        m = CENTRAL_RE.match(key)
        if m is None:
            raise ValueError(f"Unparsed central key: {key}")
        point = (m.group(1),) + tuple(_num(g) for g in m.groups()[1:])
        if point in table:
            raise ValueError(f"Duplicate central point {point}: {key}")
        eff, err = (float(v) for v in str(val).split("+-"))
        table[point] = (eff, err)
    return table


def sample_point(name):
    """Map our sample name to (channel, mXX, mA, ctau)."""
    m = SAMPLE_RE.match(name)
    if m is None:
        raise ValueError(f"Unparsed sample name: {name}")
    return (m.group(1),) + tuple(_num(g) for g in m.groups()[1:])


def v10_samples(paths):
    """Sorted sample names of the v10 ntuple yamls."""
    names = []
    for path in paths:
        with open(path, encoding="utf8") as f:
            names += list(yaml.safe_load(f)["llpNanoAOD_v2"]["samples"])
    return sorted(names)


def table_entry(central_eff, central_err, n_all, n_pass):
    """One sample's entry: central, binomial re-pass, and v10 = central / re-pass."""
    value = n_pass / n_all
    value_err = math.sqrt(value * (1 - value) / n_all)
    v10_eff = central_eff / value
    v10_err = v10_eff * math.hypot(central_err / central_eff, value_err / value)
    return {
        "central": {"eff": _sig4(central_eff), "err": _sig4(central_err)},
        "repass_5gev": {"n_all": n_all, "n_pass": n_pass,
                        "value": _sig4(value), "err": _sig4(value_err)},
        "v10": {"eff": _sig4(v10_eff), "err": _sig4(v10_err)},
    }


def build_table(samples, central, repass):
    """Table keyed by sample name, and the list of samples lacking an input."""
    table, missing = {}, []
    for name in samples:
        point = sample_point(name)
        if point not in central or name not in repass:
            missing.append(name)
            continue
        counts = repass[name][REPASS_CHANNEL]
        table[name] = table_entry(*central[point], counts["n_all"], counts["n_pass"])
    return table, missing


def build_from_inputs(repo=REPO):
    """Table and unmatched samples from the committed inputs under repo."""
    samples = v10_samples(os.path.join(repo, p) for p in NTUPLES)
    central = central_by_point(os.path.join(repo, CENTRAL))
    with open(os.path.join(repo, REPASS), encoding="utf8") as f:
        repass = json.load(f)
    return build_table(samples, central, repass)


def input_ids(repo=REPO):
    """Header lines naming each input by its git blob id (git hash-object <file>)."""
    lines = []
    for path in [CENTRAL, REPASS] + NTUPLES:
        with open(os.path.join(repo, path), "rb") as f:
            content = f.read()
        blob = hashlib.sha1(b"blob %d\0" % len(content) + content).hexdigest()
        lines.append(f"#   {blob[:12]} {path}\n")
    return lines


def main():
    """Build the table and write it, or exit without writing if a sample is unmatched."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    parser.add_argument("-o", "--output", default=OUTPUT, help="output yaml, relative to the repo root")
    args = parser.parse_args()

    table, missing = build_from_inputs()
    if missing:
        sys.exit(f"Not writing {args.output}: {len(missing)} v10 samples unmatched: {missing}")

    command = "python sidm/scripts/make_signal_gen_filter_eff.py"
    if args.output != OUTPUT:
        command += f" -o {args.output}"
    with open(os.path.join(REPO, args.output), "w", encoding="utf8") as f:
        f.write(HEADER)
        f.write(f"# Made by: {command}, from the inputs (git blob ids)\n")
        f.writelines(input_ids())
        yaml.safe_dump(table, f, sort_keys=True, default_flow_style=None)
    print(f"Wrote {len(table)} samples to {args.output}")


if __name__ == "__main__":
    main()
