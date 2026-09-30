"""Module to define classes and methods that are helpful for scaleout"""

import io
import os
import subprocess
import warnings
import zipfile
from pathlib import Path

import dask
from dask.distributed import Client, PipInstall
from distributed.diagnostics.plugin import UploadDirectory


def make_dask_client(address):
    """Create dask client that includes dependency installer (coffea-casa)."""
    dependencies = [
        "git+https://github.com/btcardwell/SIDM.git",
    ]
    client = Client(address)
    client.register_plugin(PipInstall(packages=dependencies, pip_options=["--upgrade", "--no-cache-dir"]))
    return client


_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_DEFAULT_LPC_CONFIG = _REPO_ROOT / "condor" / "lpc_condor_config"
_DEFAULT_LPC_IMAGE = (
    "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/"
    "coffea-dask-almalinux9:2025.5.0.rc2-py3.11"
)
_DEFAULT_SIDM_LOCAL_DIR = _REPO_ROOT / "sidm"
_PROXY_RENEW_CMD = "voms-proxy-init --valid 192:00 -voms cms"

# What NOT to ship to workers via UploadDirectory. Workers import sidm.tools,
# sidm.definitions and sidm.scripts (see sidm/scripts/census_dask.py) and read
# sidm/configs + sidm/data; they never open a notebook.
#
# UploadDirectory zips the tree into memory, and the scheduler -- which for
# LPCCondorCluster runs inside this very notebook process -- keeps both a
# deserialized and a pickled copy for the life of the cluster and serializes one
# more per nanny that registers, including every worker an adaptive cluster
# starts later. Peak cost in the notebook process is therefore roughly
# (2 + n_workers) x payload. sidm/studies alone is several hundred MB of
# committed notebook outputs, enough to OOM a shared cmslpc interactive node.
#
# NOTE: distributed matches skip_words against components of the ABSOLUTE path
# (distributed/diagnostics/plugin.py, os.walk root joined before the test), not
# the archive-relative one, so a checkout living under a directory named e.g.
# "studies" or "docs" would match every file and ship nothing. _make_upload_plugin
# checks the built payload rather than trusting the skip list.
_UPLOAD_SKIP_WORDS = (
    # UploadDirectory's own defaults, which are replaced rather than extended
    # when skip_words is passed explicitly.
    ".git", ".github", ".pytest_cache", "tests", "docs",
    # SIDM: notebooks and their committed outputs.
    "studies", "test_notebooks", "__pycache__", ".ipynb_checkpoints",
)
# Only compiled artifacts and stray notebooks. Deliberately NOT .root/.coffea:
# the directory skips above already remove every committed output, and sidm/data
# is where a ROOT-format correction payload would land if one is ever added.
_UPLOAD_SKIP_EXTS = (".pyc", ".ipynb")
# Warn when the projected peak, (2 + max_workers) x payload, exceeds this.
_UPLOAD_PEAK_BUDGET_BYTES = 512 * 1024**2
# Must survive the skip list, or the upload is silently useless to a worker.
_UPLOAD_SENTINEL = "tools/sidm_processor.py"


def check_voms_proxy(min_seconds_left=3600):
    """Verify a CMS VOMS proxy exists and is not about to expire.

    Inspects $X509_USER_PROXY (or /tmp/x509up_u<UID> if unset) and asks
    voms-proxy-info how much time is left. Raises RuntimeError with the
    renewal command if the proxy is missing or expiring within
    ``min_seconds_left`` seconds. Sets X509_USER_PROXY in os.environ so that
    downstream tooling (LPCCondorCluster, condor_submit) picks up the same file.

    Returns the proxy file path on success.
    """
    proxy = os.environ.get("X509_USER_PROXY") or f"/tmp/x509up_u{os.getuid()}"
    if not os.path.isfile(proxy):
        raise RuntimeError(
            f"No VOMS proxy found at {proxy}. Renew it on cmslpc with:\n"
            f"  {_PROXY_RENEW_CMD}"
        )

    result = subprocess.run(
        ["voms-proxy-info", "-file", proxy, "-timeleft"],
        capture_output=True, text=True,
    )
    try:
        remaining = int(result.stdout.strip())
    except ValueError:
        raise RuntimeError(
            f"voms-proxy-info on {proxy} did not return a parseable lifetime "
            f"(stdout={result.stdout!r}, stderr={result.stderr!r}). Renew with:\n"
            f"  {_PROXY_RENEW_CMD}"
        )

    if remaining < min_seconds_left:
        hrs = remaining / 3600.0
        raise RuntimeError(
            f"VOMS proxy at {proxy} has only {remaining}s left ({hrs:.1f}h, "
            f"threshold {min_seconds_left}s). Renew with:\n"
            f"  {_PROXY_RENEW_CMD}"
        )

    os.environ["X509_USER_PROXY"] = proxy
    return proxy


def _make_upload_plugin(sidm_local_dir, max_workers):
    """Build the UploadDirectory plugin that ships the local sidm/ tree to workers.

    Excludes notebooks and their committed outputs (see _UPLOAD_SKIP_WORDS), which
    workers never read and which dominate the tree's size. Raises if the payload
    came out empty of SIDM code, and warns when the projected peak memory in this
    process is large enough to threaten the interactive node.
    """
    plugin = UploadDirectory(
        str(sidm_local_dir),
        restart_workers=False,
        skip_words=_UPLOAD_SKIP_WORDS,
        skip=(lambda fn: os.path.splitext(fn)[1] in _UPLOAD_SKIP_EXTS,),
    )

    # A skip list applied to absolute paths can silently match everything, and a
    # missing directory zips to nothing at all. Both yield an empty payload, and
    # the only symptom would be ModuleNotFoundError on a worker much later.
    root = os.path.basename(str(sidm_local_dir).rstrip(os.sep))
    names = zipfile.ZipFile(io.BytesIO(plugin.data)).namelist()
    if f"{root}/{_UPLOAD_SENTINEL}" not in names:
        raise RuntimeError(
            f"UploadDirectory payload for {sidm_local_dir} does not contain "
            f"{root}/{_UPLOAD_SENTINEL} ({len(names)} files). Workers would get no "
            f"SIDM code and fail with ModuleNotFoundError. Either the path is wrong "
            f"or one of its absolute-path components matches _UPLOAD_SKIP_WORDS "
            f"{_UPLOAD_SKIP_WORDS} -- distributed applies skip_words to the whole "
            f"absolute path, so a checkout under a directory named e.g. 'studies' or "
            f"'docs' matches every file."
        )

    projected = len(plugin.data) * (max_workers + 2)
    if projected > _UPLOAD_PEAK_BUDGET_BYTES:
        warnings.warn(
            f"UploadDirectory payload for {sidm_local_dir} is "
            f"{len(plugin.data) / 1024**2:.0f} MB, and the scheduler runs in this "
            f"notebook process: it holds two copies and sends one more per worker, "
            f"so max_workers={max_workers} projects to about "
            f"{projected / 1024**2:.0f} MB. That can exhaust a shared cmslpc "
            f"interactive node and get this kernel OOM-killed. Remove large files "
            f"from the sidm/ tree, extend scaleout._UPLOAD_SKIP_WORDS, or lower "
            f"max_workers.",
            RuntimeWarning,
            stacklevel=3,
        )
    return plugin


def make_lpc_client(
    min_workers=1,
    max_workers=10,
    memory="4GB",
    disk="4GB",
    cores=1,
    death_timeout=600,
    image=_DEFAULT_LPC_IMAGE,
    sidm_local_dir=_DEFAULT_SIDM_LOCAL_DIR,
    condor_config=_DEFAULT_LPC_CONFIG,
    dashboard_port=None,
    **cluster_kwargs,
):
    """Create an LPCCondorCluster + Client for scaling SIDM jobs from a notebook on cmslpc.

    Workers run as Condor jobs inside the coffea-dask apptainer image; the
    client (the notebook process) stays in the sidm_venv outside the apptainer.
    The local sidm/ tree is shipped to each worker via UploadDirectory, so
    uncommitted edits are visible to workers without a git push.

    Args:
        min_workers, max_workers: passed to cluster.adapt(). Default min=1
            keeps a worker around for immediate response. Pass min_workers=0
            for fully adaptive behavior (workers only when work is queued).
        memory, disk: per-worker resource request.
        cores: per-worker CPU cores (also dask threads).
        death_timeout: seconds a worker waits to contact the scheduler before
            self-terminating. Raised from the dask default of 60 because LPC
            condor queues can leave workers idle for several minutes.
        image: apptainer image for workers. Default matches the coffea version
            in sidm_venv (2025.5.0rc2, py3.11). For other coffea versions, look
            under /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/.
        sidm_local_dir: path to a local sidm/ source tree to upload to each
            worker. The local checkout (including uncommitted changes) is what
            workers see, but only sidm/{tools,definitions,scripts,configs,data}
            reaches them: notebooks and their committed outputs are excluded
            (_UPLOAD_SKIP_WORDS), because shipping sidm/studies adds hundreds of
            MB that the in-process scheduler holds and re-sends to every worker.
            Code that runs inside a dask task must therefore live in one of those
            directories, not next to a notebook in sidm/studies/. Set to None to
            skip the upload entirely, but note the apptainer image does not
            contain sidm, so workers would then have no SIDM code at all.
        condor_config: path to a CONDOR_CONFIG file. Defaults to
            condor/lpc_condor_config in this repo, a minimal LPC interactive
            config that omits the cmslpc-local-conf.py include directive
            (which references a per-user file that does not exist on every node).
        dashboard_port: pin the scheduler's Bokeh dashboard to this port on
            the LPC side. None (default) uses dask's default of 8787. Pass an
            int, e.g. 8790, to move the dashboard off 8787 when another one is
            already there (coffea-casa also serves its Dask dashboard on
            8787), then forward that same port in your
            `ssh -L <port>:localhost:<port>` tunnel; or pass 0 to let dask
            pick any free port. Any requested port (8787 or a pinned one) is
            only a request: if it is already in use on the node, distributed
            warns and rebinds to a random free port, so forward whatever port
            `cluster.dashboard_link` actually prints.
        **cluster_kwargs: forwarded to LPCCondorCluster.

    Returns:
        (cluster, client) tuple. Caller is responsible for cluster.close().
    """
    check_voms_proxy()

    # Set CONDOR_CONFIG before importing lpcjobqueue: htcondor caches its
    # config at module-import time, so the env var must be in place first.
    os.environ["CONDOR_CONFIG"] = str(condor_config)

    from lpcjobqueue import LPCCondorCluster

    # Point the dashboard link at localhost so it matches the documented SSH
    # tunnel (`ssh -L 8787:localhost:8787`, README step 7). Importing lpcjobqueue
    # above sets distributed.dashboard.link to the relative "/proxy/{port}/status"
    # -- a jupyter-server-proxy route that does not exist on a plain SSH-forwarded
    # LPC scheduler (the dashboard serves /status at the root), so that link 404s
    # over the tunnel. Rewrite it unless we are genuinely under a JupyterHub, where
    # the proxy route is real. This must run AFTER the import, which would
    # otherwise clobber the setting.
    if "JUPYTERHUB_SERVICE_PREFIX" not in os.environ:
        dask.config.set({"distributed.dashboard.link": "{scheme}://localhost:{port}/status"})

    # Pin the scheduler's Bokeh dashboard to a specific port when asked, e.g. to
    # avoid a clash with another dashboard already on 8787 (coffea-casa serves its
    # Dask dashboard there too). The same port must be forwarded in the
    # `ssh -L <port>:localhost:<port>` tunnel, because cluster.dashboard_link
    # reports the scheduler's actual port. Merge into any caller-supplied
    # scheduler_options without clobbering an explicit dashboard_address; left at
    # dask's default when dashboard_port is None.
    if dashboard_port is not None:
        cluster_kwargs["scheduler_options"] = {
            "dashboard_address": f":{dashboard_port}",
            **(cluster_kwargs.pop("scheduler_options", None) or {}),
        }

    cluster = LPCCondorCluster(
        memory=memory,
        disk=disk,
        cores=cores,
        death_timeout=death_timeout,
        image=image,
        ship_env=False,
        **cluster_kwargs,
    )
    cluster.adapt(minimum=min_workers, maximum=max_workers)
    client = Client(cluster)

    if sidm_local_dir is not None:
        client.register_plugin(_make_upload_plugin(sidm_local_dir, max_workers))

    return cluster, client
