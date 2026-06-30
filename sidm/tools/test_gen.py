#!/usr/bin/env python3
"""Standalone test for the cycle-safe distinctParent kernel (no network).

    python sidm/tools/test_gen.py

Checks the capped kernel reproduces coffea's distinctParent on acyclic chains AND terminates
(returns -1, no hang) on a self-referential or multi-particle cycle.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen  # noqa: E402 -- importing also installs the kernel patch

_P = _F = 0


def check(name, cond):
    global _P, _F
    if cond:
        _P += 1
    else:
        _F += 1
        print("FAIL:", name)


def test_kernel():
    k = gen._capped_distinctParent_kernel
    # chain mu(13) <- mu(13) <- gamma(22): particles 0,1 have distinct parent = index 2 (gamma)
    out = k(np.array([1, 2, -1], dtype=np.int64), np.array([13, 13, 22], dtype=np.int64))
    check("distinct parent of same-pdg chain is the first different-pdg ancestor",
          out[0] == 2 and out[1] == 2)
    check("particle with no mother -> -1", out[2] == -1)
    # immediate different-pdg mother IS the distinct parent
    out = k(np.array([1, -1], dtype=np.int64), np.array([13, 22], dtype=np.int64))
    check("immediate distinct mother", out[0] == 1)
    # self-referential cycle (mother == own index, same pdg) must TERMINATE -> -1
    out = k(np.array([0], dtype=np.int64), np.array([0], dtype=np.int64))
    check("self-referential cycle terminates -> -1 (no hang)", out[0] == -1)
    # a 2-particle same-pdg cycle (A<->B) also terminates
    out = k(np.array([1, 0], dtype=np.int64), np.array([5, 5], dtype=np.int64))
    check("2-particle cycle terminates -> -1", out[0] == -1 and out[1] == -1)
    # an out-of-range mother index (>= length) is treated as no distinct parent, never an oob store
    out = k(np.array([5], dtype=np.int64), np.array([13], dtype=np.int64))
    check("out-of-range mother -> -1 (no oob)", out[0] == -1)


if __name__ == "__main__":
    test_kernel()
    print(f"\n{_P} passed, {_F} failed")
    sys.exit(1 if _F else 0)
