"""Sidecar metadata for .coffea outputs.

``write_run_metadata`` produces a ``.meta.yaml`` sidecar next to a ``.coffea``
recording the selection definitions, hist-collection contents, input ROOT
file list, per-sample cross section, SIDM git commit, coffea version, schema,
chunksize, ``unweighted_hist`` flag, and a UTC timestamp. ``load_run_metadata``
reads the sidecar back.

Used by both the dask-from-notebook path
(``sidm/test_notebooks/lpc_dask_example.ipynb``) and the Condor batch path
(``sidm/scripts/merge_coffea_chunks_eos.py``).
"""

import contextlib
import io
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import coffea
import yaml

from sidm.tools.utilities import get_xs, load_yaml


def write_run_metadata(
    coffea_path,
    *,
    fileset,
    selections,
    hist_collections,
    schema=None,
    chunksize=None,
    unweighted_hist=None,
    extra=None,
    sidm_root=None,
):
    """Write a ``.meta.yaml`` sidecar describing what produced ``coffea_path``.

    Args:
        coffea_path: path to the ``.coffea`` output. The sidecar path is derived
            from this: ``foo.coffea`` -> ``foo.meta.yaml``.
        fileset: coffea-style ``{sample_name: {"files": [...], "metadata":
            {...}}}`` dict, the same shape passed to ``processor.Runner.run``.
            Recorded under ``samples`` in the sidecar.
        selections: list of selection names from
            ``sidm/configs/selections.yaml``. The full nested definition of each
            is copied into the sidecar.
        hist_collections: list of histogram-collection names from
            ``sidm/configs/hist_collections.yaml``. The full list of histogram
            names in each collection is expanded into the sidecar.
        schema: optional name of the awkward/coffea schema
            (e.g. ``"LLPNanoAODSchema"``).
        chunksize: optional chunk size used by ``processor.Runner``.
        unweighted_hist: optional bool — whether ``SidmProcessor`` was
            constructed with ``unweighted_hist=True``.
        extra: optional dict of additional fields to merge into the sidecar
            (e.g. condor cluster ID, dask cluster dashboard URL).
        sidm_root: optional repo-root override; default walks up from this file.

    Returns:
        Path to the written sidecar.
    """
    sidm_root = Path(sidm_root) if sidm_root else _find_sidm_root()

    selections_cfg = load_yaml(str(sidm_root / "sidm" / "configs" / "selections.yaml"))
    hists_cfg = load_yaml(str(sidm_root / "sidm" / "configs" / "hist_collections.yaml"))

    selections_used = [
        {"name": s, "definition": _yaml_safe(selections_cfg.get(s, "<not found>"))}
        for s in selections
    ]
    hist_collections_used = [
        {"name": hc, "hists": _yaml_safe(hists_cfg.get(hc, []))}
        for hc in hist_collections
    ]

    samples = []
    for name, info in fileset.items():
        files = list(info.get("files", []))
        # Cross-section from sidm/configs/cross_sections.yaml (pb).
        # utilities.get_xs falls back to 0.001 pb for 2Mu2E/4Mu signals not
        # explicitly listed, and raises for unknown bg/data — store None then.
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                xsec_pb = float(get_xs(name))
        except Exception:
            xsec_pb = None
        samples.append({
            "name": name,
            "n_files": len(files),
            "xsec_pb": xsec_pb,
            "metadata": _yaml_safe(info.get("metadata", {})),
            "files": files,
        })

    meta = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "coffea_path": str(coffea_path),
        "n_samples": len(samples),
        "selections": selections_used,
        "hist_collections": hist_collections_used,
        "schema": schema,
        "chunksize": chunksize,
        "unweighted_hist": unweighted_hist,
        "sidm_commit": _git_rev(sidm_root),
        "coffea_version": coffea.__version__,
        "samples": samples,
    }
    if extra:
        meta.update(extra)

    sidecar = _sidecar_path(coffea_path)
    with open(sidecar, "w") as f:
        yaml.safe_dump(meta, f, sort_keys=False, default_flow_style=False, width=200)
    return sidecar


def load_run_metadata(coffea_path):
    """Load and return the ``.meta.yaml`` sidecar next to ``coffea_path``."""
    sidecar = _sidecar_path(coffea_path)
    with open(sidecar) as f:
        return yaml.safe_load(f)


def _sidecar_path(coffea_path):
    p = str(coffea_path)
    if p.endswith(".coffea"):
        p = p[: -len(".coffea")]
    return p + ".meta.yaml"


def _find_sidm_root():
    return Path(__file__).resolve().parent.parent.parent


def _git_rev(repo_root):
    try:
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            timeout=5,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return None


def _yaml_safe(obj):
    """Convert a PyYAML-loaded structure into something ``safe_dump`` will
    re-serialize cleanly (no Python-specific tags)."""
    if isinstance(obj, dict):
        return {k: _yaml_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_yaml_safe(x) for x in obj]
    return obj
