"""Browse already-run SIDM coffea outputs on EOS.

Every merged ``.coffea`` in the shared SIDM area carries a ``.meta.yaml``
sidecar (see ``sidm.tools.metadata``) recording the selections, histogram
collections, input samples, and provenance that produced it. This module
crawls those sidecars under
``/store/group/lpcmetx/SIDM/coffea_outputs/<producer>/<study>/`` and builds a
single pandas index, so you can find outputs -- your own or a collaborator's --
by selection, histogram, sample, input file, producer, study, commit, schema,
or date, and recover the exact ``root://`` path (and a ready-to-run ``xrdcp``
command) for each.

The crawl derives the EOS ``.coffea`` location from each sidecar's own EOS
path (``foo.meta.yaml`` -> ``foo.coffea``); the ``coffea_path`` field stored
*inside* the sidecar records the producer's local working path at write time
and is not the shared-EOS location.

Requires a Kerberos ticket (read access to ``/store/group/lpcmetx``). Intended
to run in a JupyterLab kernel on LPC or EAF; the interactive widget front-end
lives in ``sidm/test_notebooks/output_browser.ipynb``.
"""

from __future__ import annotations

import concurrent.futures
import json
import re
import subprocess
from pathlib import PurePosixPath

import pandas as pd
import yaml

REDIRECTOR = "root://cmseos.fnal.gov"
EOS_BASE = "/store/group/lpcmetx/SIDM/coffea_outputs"

# Preferred column order for the index DataFrame. Any column not listed here
# (e.g. future sidecar fields) is appended after these, in discovery order.
INDEX_COLUMNS = [
    "producer", "study", "name", "created_utc", "created",
    "n_samples", "samples", "selections", "n_selections",
    "hist_collections", "hists", "n_hists",
    "n_input_files", "input_files",
    "schema", "chunksize", "unweighted_hist", "n_chunks",
    "sidm_commit", "commit_short", "coffea_version",
    "years", "is_data",
    "coffea_size", "size_human", "coffea_mtime", "coffea_exists",
    "coffea_url", "meta_url", "coffea_path", "meta_path",
    "sample_details",
]

# One `xrdfs ... ls -l -R` line, e.g.
#   -rw-r--r-- murtazas us_cms 70190 2026-06-03 02:06:03 /store/.../x.meta.yaml
# The owner/group columns are not emitted by every xrootd server, so the size
# is anchored as the integer immediately preceding the "date time" stamp, and
# the path is everything after it (verbatim, so a path may contain spaces).
_LS_LINE = re.compile(
    r"^(?P<perms>\S+)\s+.*?(?P<size>\d+)\s+"
    r"(?P<date>\d{4}-\d{2}-\d{2})\s+(?P<time>\d{2}:\d{2}:\d{2})\s+(?P<path>.+)$"
)


# --------------------------------------------------------------------------- #
# Data layer: crawl EOS sidecars -> records -> pandas index -> query.
# --------------------------------------------------------------------------- #
def _run(cmd, timeout=180):
    """Run ``cmd`` and return its stdout (text). Raise on non-zero exit."""
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr.strip()}"
        )
    return proc.stdout


def _parse_ls_long(text):
    """Parse ``xrdfs ... ls -l -R`` output into a list of entry dicts.

    Lines that do not match the expected ``perms [owner group] size date time
    path`` shape are skipped. Directory lines (perms starting with ``d``) are
    kept but flagged via ``is_dir`` so callers can ignore them.
    """
    entries = []
    for line in text.splitlines():
        m = _LS_LINE.match(line.rstrip())
        if not m:
            continue
        try:
            size = int(m["size"])
        except ValueError:
            size = None
        entries.append({
            "path": m["path"],
            "perms": m["perms"],
            "is_dir": m["perms"].startswith("d"),
            "size": size,
            "mtime": f"{m['date']} {m['time']}",
        })
    return entries


def _record(meta_path, meta, by_path, base, redirector):
    """Build one index record from a parsed sidecar and the file listing."""
    rel = PurePosixPath(meta_path).relative_to(base)
    parts = rel.parts  # (producer, study..., name.meta.yaml)
    producer = parts[0] if parts else ""
    study = "/".join(parts[1:-1])  # everything between producer and filename
    name = PurePosixPath(meta_path).name[: -len(".meta.yaml")]
    coffea_path = meta_path[: -len(".meta.yaml")] + ".coffea"
    coffea_entry = by_path.get(coffea_path)
    meta = meta or {}

    sel_meta = meta.get("selections") or []
    hc_meta = meta.get("hist_collections") or []
    samples_meta = meta.get("samples") or []

    selections = [s.get("name") for s in sel_meta if isinstance(s, dict)]
    hist_collections = [h.get("name") for h in hc_meta if isinstance(h, dict)]
    hists = sorted({
        h for hc in hc_meta if isinstance(hc, dict) for h in (hc.get("hists") or [])
    })
    samples = [s.get("name") for s in samples_meta if isinstance(s, dict)]
    input_files = [
        f for s in samples_meta if isinstance(s, dict) for f in (s.get("files") or [])
    ]
    years = sorted({
        str((s.get("metadata") or {}).get("year"))
        for s in samples_meta
        if isinstance(s, dict) and (s.get("metadata") or {}).get("year") is not None
    })
    is_data = any(
        bool((s.get("metadata") or {}).get("is_data"))
        for s in samples_meta if isinstance(s, dict)
    )

    return {
        "producer": producer,
        "study": study,
        "name": name,
        "created_utc": meta.get("created_utc"),
        "n_samples": meta.get("n_samples", len(samples)),
        "samples": samples,
        "sample_details": samples_meta,
        "selections": selections,
        "hist_collections": hist_collections,
        "hists": hists,
        "input_files": input_files,
        "n_input_files": len(input_files),
        "schema": meta.get("schema"),
        "chunksize": meta.get("chunksize"),
        "unweighted_hist": meta.get("unweighted_hist"),
        "n_chunks": meta.get("n_chunks"),
        "sidm_commit": meta.get("sidm_commit"),
        "coffea_version": meta.get("coffea_version"),
        "years": years,
        "is_data": is_data,
        "coffea_path": coffea_path,
        "meta_path": meta_path,
        "coffea_url": f"{redirector}/{coffea_path}",
        "meta_url": f"{redirector}/{meta_path}",
        "coffea_exists": coffea_entry is not None,
        "coffea_size": coffea_entry["size"] if coffea_entry else None,
        "coffea_mtime": coffea_entry["mtime"] if coffea_entry else None,
    }


def crawl(base=EOS_BASE, redirector=REDIRECTOR, producers=None, max_workers=8, verbose=False):
    """Crawl ``.meta.yaml`` sidecars under ``base`` and return a list of records.

    One ``xrdfs ls`` lists every file; then one short-lived ``xrdfs cat`` per
    sidecar reads it (parallelised over ``max_workers`` threads). The crawl
    therefore scales linearly in subprocess spawns -- for a large area pass
    ``producers=`` to scope it, or use ``build_index(cache=...)`` to avoid
    re-crawling.

    Args:
        base: EOS directory holding ``<producer>/<study>/`` subtrees.
        redirector: xrootd redirector URL (FNAL EOS by default).
        producers: optional iterable of producer names to restrict the crawl to
            (e.g. ``["murtazas"]``); ``None`` crawls every producer.
        max_workers: thread count for the per-sidecar ``xrdfs cat`` fetches.
        verbose: print skipped paths and a final count.

    Returns:
        List of record dicts, one per sidecar successfully read.

    Raises:
        RuntimeError: if *every* listing fails (e.g. an expired Kerberos ticket
            or unreachable redirector) -- so an environment problem is not
            mistaken for an empty area.
    """
    roots = [f"{base}/{p}" for p in producers] if producers else [base]
    entries = []
    n_listed = 0
    for root in roots:
        try:
            entries.extend(_parse_ls_long(_run(["xrdfs", redirector, "ls", "-l", "-R", root])))
            n_listed += 1
        except RuntimeError as exc:
            if verbose:
                print(f"skip listing {root}: {exc}")
    if n_listed == 0:
        raise RuntimeError(
            f"could not list any of {roots} on {redirector}. Check your Kerberos "
            "ticket (klist) and EOS access, then retry."
        )

    by_path = {e["path"]: e for e in entries}
    meta_paths = sorted(
        p for p, e in by_path.items() if p.endswith(".meta.yaml") and not e["is_dir"]
    )

    def _fetch(meta_path):
        try:
            meta = yaml.safe_load(_run(["xrdfs", redirector, "cat", meta_path]))
            return _record(meta_path, meta, by_path, base, redirector)
        except Exception as exc:  # noqa: BLE001 - one bad sidecar must not abort the crawl
            if verbose:
                print(f"skip {meta_path}: {exc}")
            return None

    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, max_workers)) as pool:
        for rec in pool.map(_fetch, meta_paths):
            if rec is not None:
                records.append(rec)
    if verbose:
        print(f"indexed {len(records)} outputs from {len(meta_paths)} sidecars under {base}")
    return records


def _human_size(n):
    if n is None or (isinstance(n, float) and pd.isna(n)):
        return ""
    n = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def _records_to_df(records):
    df = pd.DataFrame(records)
    if df.empty:
        return pd.DataFrame(columns=INDEX_COLUMNS)
    df["created"] = pd.to_datetime(df["created_utc"], utc=True, errors="coerce")
    df["commit_short"] = df["sidm_commit"].apply(
        lambda c: c[:8] if isinstance(c, str) else None
    )
    df["n_selections"] = df["selections"].apply(len)
    df["n_hists"] = df["hists"].apply(len)
    df["size_human"] = df["coffea_size"].apply(_human_size)
    df = df.sort_values("created", ascending=False, na_position="last").reset_index(drop=True)
    cols = [c for c in INDEX_COLUMNS if c in df.columns]
    cols += [c for c in df.columns if c not in cols]
    return df[cols]


def build_index(base=EOS_BASE, redirector=REDIRECTOR, producers=None,
                cache=None, refresh=False, max_workers=8, verbose=False):
    """Crawl EOS (or load a cached crawl) and return the output index.

    Args:
        cache: optional path to a JSON cache of the crawl. When given and the
            file exists (and ``refresh`` is False) the index is loaded from it
            instead of re-crawling EOS; a fresh, non-empty crawl is written back
            to it. The cache is NOT auto-invalidated -- pass ``refresh=True`` to
            re-read after new outputs have been written.
        refresh: force a re-crawl even if ``cache`` exists.

    Returns:
        ``pandas.DataFrame`` with one row per output (see ``INDEX_COLUMNS``).
    """
    records = None
    if cache and not refresh:
        try:
            with open(cache) as fh:
                records = json.load(fh)
        except (OSError, json.JSONDecodeError):
            records = None
    if records is None:
        records = crawl(
            base=base, redirector=redirector, producers=producers,
            max_workers=max_workers, verbose=verbose,
        )
        if cache and records:  # never persist an empty/failed crawl over a good cache
            try:
                with open(cache, "w") as fh:
                    json.dump(records, fh)
            except OSError:
                pass
    return _records_to_df(records)


def query(df, *, producer=None, study=None, selection=None, hist=None, sample=None,
          input_file=None, hist_collection=None, commit=None, schema=None,
          coffea_version=None, created_after=None, created_before=None,
          unweighted_hist=None, is_data=None):
    """Filter the index ``df`` and return the matching rows.

    All text filters are case-insensitive **substring** matches. Scalar filters
    (``producer``, ``study``, ``commit``, ``schema``, ``coffea_version``) match
    against the corresponding column; list filters (``selection``, ``hist``,
    ``sample``, ``input_file``, ``hist_collection``) keep a row if ANY element
    of that row's list contains the query. Every text filter also accepts a
    list of substrings -- a row then matches if it satisfies ANY of them.

    Filters are combined with **AND**: a row is kept only if it satisfies every
    filter passed. (Within one list-valued filter the substrings are ORed.)

    Empty / whitespace-only substrings are ignored. ``created_after`` /
    ``created_before`` accept anything ``pandas.to_datetime`` understands
    (e.g. ``"2026-06-01"``); an unparseable value is ignored rather than
    raising. ``unweighted_hist`` / ``is_data`` are exact bool matches.
    """
    if df.empty:
        return df
    mask = pd.Series(True, index=df.index)

    def _scalar(col, val):
        nonlocal mask
        if val is None:
            return
        vals = [val] if isinstance(val, str) else list(val)
        vals = [str(v) for v in vals if str(v).strip()]
        if not vals:
            return
        col_l = df[col].fillna("").astype(str).str.lower()
        hit = pd.Series(False, index=df.index)
        for v in vals:
            hit |= col_l.str.contains(v.lower(), regex=False)
        mask &= hit

    def _listcol(col, val):
        nonlocal mask
        if val is None:
            return
        vals = [val] if isinstance(val, str) else list(val)
        vals = [str(v).lower() for v in vals if str(v).strip()]
        if not vals:
            return
        def _hit(lst):
            low = [str(x).lower() for x in (lst or [])]
            return any(any(v in x for x in low) for v in vals)
        mask &= df[col].apply(_hit)

    def _date(val):
        if val is None:
            return None
        ts = pd.to_datetime(val, utc=True, errors="coerce")
        return None if pd.isna(ts) else ts

    _scalar("producer", producer)
    _scalar("study", study)
    _scalar("sidm_commit", commit)
    _scalar("schema", schema)
    _scalar("coffea_version", coffea_version)
    _listcol("selections", selection)
    _listcol("hists", hist)
    _listcol("samples", sample)
    _listcol("input_files", input_file)
    _listcol("hist_collections", hist_collection)

    after = _date(created_after)
    if after is not None:
        mask &= df["created"] >= after
    before = _date(created_before)
    if before is not None:
        mask &= df["created"] <= before
    if unweighted_hist is not None:
        mask &= df["unweighted_hist"] == unweighted_hist
    if is_data is not None:
        mask &= df["is_data"] == is_data

    return df[mask]


def summary(df):
    """Compact, human-readable view of the index (heavy list columns collapsed)."""
    if df.empty:
        return df
    return pd.DataFrame({
        "producer": df["producer"],
        "study": df["study"],
        "name": df["name"],
        "n_samp": df["n_samples"],
        "selections": df["selections"].apply(lambda x: ", ".join(map(str, x or []))),
        "hist_colls": df["hist_collections"].apply(lambda x: ", ".join(map(str, x or []))),
        "n_hists": df["n_hists"],
        "created": df["created_utc"],
        "size": df["size_human"],
        "commit": df["commit_short"],
    }).reset_index(drop=True)


def xrdcp_command(row, dest=".", which="coffea"):
    """Return an ``xrdcp`` command string to fetch ``row``'s coffea (or meta) file."""
    url = row["coffea_url"] if which == "coffea" else row["meta_url"]
    return f"xrdcp -f {url} {dest}"


def load_coffea(row, tmpdir=None):
    """``xrdcp`` ``row``'s ``.coffea`` to a temp dir and ``coffea.util.load`` it.

    ``coffea.util.load`` opens a local path, so the file is staged locally
    first (it cannot read a ``root://`` URL directly). The staged file is left
    in ``tmpdir`` (a fresh temp dir if none is given) and is not cleaned up.
    """
    import os
    import tempfile
    from coffea.util import load

    tmpdir = tmpdir or tempfile.mkdtemp()
    local = os.path.join(tmpdir, os.path.basename(row["coffea_path"]))
    _run(["xrdcp", "-f", row["coffea_url"], local])
    return load(local)


# --------------------------------------------------------------------------- #
# UI layer: interactive ipywidgets browser over the index (lazy import).
# --------------------------------------------------------------------------- #
def browser(df, dest="."):
    """Build and return an interactive ipywidgets browser over the index ``df``.

    Renders filter controls above a live results table that updates as you
    change any control:

    - Producer, Study, Selection, and Hist-collection **multi-selects** keep a
      row if it matches any *exact* value you pick.
    - Histogram, Sample, and Input-file **text boxes** are case-insensitive
      *substring* searches; Commit and Schema and a created-date range refine
      further.

    An *Inspect an output* dropdown lists the matching outputs; picking one
    shows its full metadata and a ready-to-run ``xrdcp`` command (targeting
    ``dest``).

    Requires a live Jupyter kernel with ipywidgets. Call ``display(...)`` on the
    returned widget (Jupyter auto-displays it if it is the last cell expression).
    """
    import ipywidgets as W
    from IPython.display import clear_output, display

    style = {"description_width": "initial"}

    def _opts(col):
        vals = set()
        for lst in df[col]:
            vals.update(lst or [])
        return sorted(vals)

    producers = sorted(df["producer"].dropna().unique())
    studies = sorted(df["study"].dropna().unique())
    selections = _opts("selections")
    hist_colls = _opts("hist_collections")
    schemas = sorted(s for s in df["schema"].dropna().unique())

    def _rows(seq):
        return min(6, max(2, len(seq)))

    w_producer = W.SelectMultiple(options=producers, description="Producer", rows=_rows(producers), style=style)
    w_study = W.SelectMultiple(options=studies, description="Study", rows=_rows(studies), style=style)
    w_selection = W.SelectMultiple(options=selections, description="Selection", rows=_rows(selections), style=style)
    w_histcoll = W.SelectMultiple(options=hist_colls, description="Hist collection", rows=_rows(hist_colls), style=style)
    w_hist = W.Text(description="Histogram", placeholder="substring, e.g. muon_dxy", style=style)
    w_sample = W.Text(description="Sample", placeholder="substring, e.g. 2Mu2E_500GeV", style=style)
    w_file = W.Text(description="Input file", placeholder="substring, e.g. MDp-1p2", style=style)
    w_commit = W.Text(description="Commit", placeholder="sidm_commit substring", style=style)
    w_schema = W.Dropdown(description="Schema", options=["(any)"] + schemas, value="(any)", style=style)
    w_after = W.Text(description="Created after", placeholder="YYYY-MM-DD", style=style)
    w_before = W.Text(description="Created before", placeholder="YYYY-MM-DD", style=style)
    w_reset = W.Button(description="Reset filters", icon="undo")

    results = W.Output()
    detail = W.Output()
    w_pick = W.Dropdown(description="Inspect", options=[], style=style, layout=W.Layout(width="90%"))

    def _apply(_=None):
        fdf = df
        if w_producer.value:
            fdf = fdf[fdf["producer"].isin(w_producer.value)]
        if w_study.value:
            fdf = fdf[fdf["study"].isin(w_study.value)]
        if w_selection.value:
            picked = set(w_selection.value)
            fdf = fdf[fdf["selections"].apply(lambda lst: bool(picked & set(lst or [])))]
        if w_histcoll.value:
            picked = set(w_histcoll.value)
            fdf = fdf[fdf["hist_collections"].apply(lambda lst: bool(picked & set(lst or [])))]
        with results:
            clear_output(wait=True)
            try:
                fdf = query(
                    fdf,
                    hist=[w_hist.value] if w_hist.value.strip() else None,
                    sample=[w_sample.value] if w_sample.value.strip() else None,
                    input_file=[w_file.value] if w_file.value.strip() else None,
                    commit=w_commit.value.strip() or None,
                    schema=None if w_schema.value == "(any)" else w_schema.value,
                    created_after=w_after.value.strip() or None,
                    created_before=w_before.value.strip() or None,
                )
            except Exception as exc:  # keep the browser responsive on a bad filter
                print(f"filter error: {exc}")
                fdf = df.iloc[0:0]
            if fdf.empty:
                print("no matching outputs")
            else:
                print(f"{len(fdf)} matching output(s)")
                display(summary(fdf))
        w_pick.options = [
            (f"{r['producer']}/{r['study']}/{r['name']}", i) for i, r in fdf.iterrows()
        ]

    def _show_detail(_=None):
        idx = w_pick.value
        with detail:
            clear_output(wait=True)
            if idx is None:
                return
            r = df.loc[idx]
            missing = "" if r["coffea_exists"] else "  [coffea MISSING on EOS]"
            print(f"{r['producer']}/{r['study']}/{r['name']}{missing}")
            print(f"created  : {r['created_utc']}    commit: {r['commit_short']}    coffea: {r['coffea_version']}")
            print(f"schema   : {r['schema']}    chunksize: {r['chunksize']}    unweighted_hist: {r['unweighted_hist']}")
            print(f"samples  : ({r['n_samples']}) " + ", ".join(map(str, r["samples"])))
            print("selections: " + ", ".join(map(str, r["selections"])))
            print("hist colls: " + ", ".join(map(str, r["hist_collections"])) + f"    ({r['n_hists']} hists)")
            print(f"inputs   : {r['n_input_files']} ROOT file(s)    size: {r['size_human']}")
            print(f"coffea   : {r['coffea_url']}")
            display(W.Textarea(
                value=xrdcp_command(r, dest=dest),
                description="fetch", style=style,
                layout=W.Layout(width="95%", height="54px"),
            ))

    def _reset(_):
        for w in (w_producer, w_study, w_selection, w_histcoll):
            w.value = ()
        for w in (w_hist, w_sample, w_file, w_commit, w_after, w_before):
            w.value = ""
        w_schema.value = "(any)"
        _apply()

    for w in (w_producer, w_study, w_selection, w_histcoll, w_hist, w_sample,
              w_file, w_commit, w_schema, w_after, w_before):
        w.observe(_apply, "value")
    w_pick.observe(_show_detail, "value")
    w_reset.on_click(_reset)

    controls = W.VBox([
        W.HBox([w_producer, w_study]),
        W.HBox([w_selection, w_histcoll]),
        W.HBox([w_hist, w_sample, w_file]),
        W.HBox([w_commit, w_schema]),
        W.HBox([w_after, w_before, w_reset]),
    ])
    _apply()
    return W.VBox([
        controls,
        W.HTML("<b>Results</b>"),
        results,
        W.HTML("<b>Inspect an output</b>"),
        w_pick,
        detail,
    ])


if __name__ == "__main__":
    # Headless smoke check: crawl EOS and print the compact index.
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 220)
    _df = build_index(verbose=True)
    print(summary(_df).to_string(index=False))
