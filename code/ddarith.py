"""Vectorised double-double arithmetic (error-free transformations, no FMA).

Used for residuals in iterative refinement of the pair-replica system, so
that the refined solution is accurate well beyond double precision even when
the M-matrix is badly conditioned (supercritical, large N).
"""
from fractions import Fraction as Fr
import numpy as np

SPLIT = 134217729.0  # 2^27 + 1


def two_sum(a, b):
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def _split(a):
    t = SPLIT * a
    hi = t - (t - a)
    return hi, a - hi


def two_prod(a, b):
    p = a * b
    ah, al = _split(a)
    bh, bl = _split(b)
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def add(xh, xl, yh, yl):
    s, e = two_sum(xh, yh)
    return two_sum(s, e + xl + yl)


def mul_d(xh, xl, d):
    p, e = two_prod(xh, d)
    return two_sum(p, e + xl * d)


def mul_dd(xh, xl, yh, yl):
    p, e = two_prod(xh, yh)
    return two_sum(p, e + (xh * yl + xl * yh))


def from_fraction(F):
    h = float(F)
    return h, float(Fr(F) - Fr(h))


def matvec(M, xh, xl):
    """Double-double M x for a CSR matrix whose entries are exact doubles (integers here)."""
    n = M.shape[0]
    ip = M.indptr
    cnt = np.diff(ip)
    ph, pl = mul_d(xh[M.indices], xl[M.indices], M.data)
    yh, yl = np.zeros(n), np.zeros(n)
    for t in range(int(cnt.max()) if n else 0):
        rows = np.nonzero(cnt > t)[0]
        pos = ip[rows] + t
        yh[rows], yl[rows] = add(yh[rows], yl[rows], ph[pos], pl[pos])
    return yh, yl
