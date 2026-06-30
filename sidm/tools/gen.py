"""Cycle-safe gen-particle ancestry for LLPNanoAOD.

Some skimmed ntuples contain self-referential GenPart entries (``genPartIdxMother`` points to
the particle's own index). coffea computes ``distinctParent`` by climbing the same-pdgId mother
chain in an unguarded ``while`` loop (``coffea.nanoevents.transforms._distinctParent_kernel``),
so that self-loop never terminates and the whole job hangs -- a single such event is enough.

Importing this module -- which the LLPNanoAOD schema does -- replaces that kernel with the
IDENTICAL algorithm plus a step cap (a cycle guard). On clean data the result is bit-for-bit
coffea's: validated at 0 mismatches over ~240k gen particles in signal and background, where the
deepest real same-pdgId chain is 5, far below the cap of 1000. A corrupt cyclic chain simply
gets no distinct parent (-1) -- the kernel guarantees termination for ANY cycle, regardless of pdgId.
There is no physics skew: on acyclic (real) ancestry the result is bit-for-bit coffea's (above), and
the only consumer of ``.distinctParent`` in the analysis is the fromGenA cut (sidm/definitions/cuts.py),
which no active channel currently applies -- so production reads it nowhere. (The observed corrupt
entries are also pdgId/status-0 fillers, not real leptons, but the safety rests on the cycle guarantee
plus non-consumption, not on assuming the garbage is always pdgId 0.) The fix is process-wide (it
patches coffea's kernel), so every ``.distinctParent`` access -- a gen source-tracing study, the deep
census under a gen-source channel, any future fromGenA selection -- becomes safe with no caller change.
"""
import numpy
import numba
from coffea.nanoevents import transforms

# Real same-pdgId ancestry chains are a handful of radiative copies (<=5 measured). The cap only
# ever triggers on a corrupt cycle, where the kernel returns -1 (no distinct parent), exactly as
# coffea returns -1 for a particle with no mother.
DISTINCT_PARENT_MAX_STEPS = 1000


@numba.njit
def _capped_distinctParent_kernel(allpart_parent, allpart_pdg):
    """coffea's _distinctParent_kernel + a step cap (cycle guard) and an in-range bound that
    replaces coffea's in-loop RuntimeError with a safe terminate. Same signature, same result on
    any acyclic chain."""
    n = len(allpart_pdg)
    out = numpy.empty(n, dtype=numpy.int64)
    for i in range(n):
        parent = allpart_parent[i]
        if parent < 0:
            out[i] = -1
            continue
        thispdg = allpart_pdg[i]
        steps = 0
        while parent >= 0 and parent < n and allpart_pdg[parent] == thispdg:
            parent = allpart_parent[parent]
            steps += 1
            if steps >= DISTINCT_PARENT_MAX_STEPS:
                parent = -1  # corrupt cycle -> no distinct parent (terminate, do not hang)
                break
        # first different-pdgId ancestor, or -1 (no mother / cycle); an out-of-range mother index
        # (should not survive local2global, but be safe) is also treated as no distinct parent.
        out[i] = parent if parent < n else -1
    return out


def install():
    """Replace coffea's distinctParent kernel with the cycle-safe one. Idempotent. Raises if the
    coffea internal moved, so a coffea upgrade fails loudly here instead of silently re-hanging."""
    if not hasattr(transforms, "_distinctParent_kernel"):
        raise RuntimeError(
            "coffea.nanoevents.transforms._distinctParent_kernel not found; the cycle-safe "
            "distinctParent patch in sidm.tools.gen needs updating for this coffea version.")
    transforms._distinctParent_kernel = _capped_distinctParent_kernel


install()


def distinct_ancestor(particles, generations=1):
    """Walk ``.distinctParent`` up ``generations`` times -- for tracing a lepton back several
    steps to its source (e.g. mu <- its distinct parent <- ...). Uses the cycle-safe kernel
    installed above, so it never hangs on the corrupt skim files. Returns a GenParticle array
    (None where the chain ends)."""
    cur = particles
    for _ in range(generations):
        cur = cur.distinctParent
    return cur
