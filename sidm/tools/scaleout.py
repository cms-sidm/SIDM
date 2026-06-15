"""Module to define classes and methods that are helpful for scaleout"""

import os
import subprocess
from pathlib import Path

import dask
from dask.distributed import Client, PipInstall


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
            workers see. Set to None to skip; the worker then sees only what is
            already in the apptainer image.
        condor_config: path to a CONDOR_CONFIG file. Defaults to
            condor/lpc_condor_config in this repo, a minimal LPC interactive
            config that omits the cmslpc-local-conf.py include directive
            (which references a per-user file that does not exist on every node).
        **cluster_kwargs: forwarded to LPCCondorCluster.

    Returns:
        (cluster, client) tuple. Caller is responsible for cluster.close().
    """
    check_voms_proxy()

    # Set CONDOR_CONFIG before importing lpcjobqueue: htcondor caches its
    # config at module-import time, so the env var must be in place first.
    os.environ["CONDOR_CONFIG"] = str(condor_config)

    from lpcjobqueue import LPCCondorCluster
    from distributed.diagnostics.plugin import UploadDirectory

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
        client.register_plugin(
            UploadDirectory(str(sidm_local_dir), restart_workers=False)
        )

    return cluster, client
