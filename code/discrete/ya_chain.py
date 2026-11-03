"""Ancestor-count chain Y of Kobayashi et al. (2018) and the Version 2 / Version 3 laws.

The number Y of cultural ancestors of a sample is a Markov chain (Wakano & Aoki 2024, Sec. 1.3.1).
From y ancestors, each draws K potential models uniformly without replacement among N and keeps
each with probability beta; by symmetry the kept set is a uniform j-subset with j ~ Bin(K, beta).
Adding newborns one at a time, the union size grows by a hypergeometric number.

Version 2: C ~ Poisson(m) with m = mu * E L, L = sum_{a>=0} Y_a (W&A eq. 2).
Version 3: C ~ mixed Poisson(mu L), Var C = m + mu^2 Var L; P(C=0) = E exp(-mu L).
Version 1 lies between them in convex order (Theorem D2); its variance is not computed here.

Usage: ya_chain.py N K beta mu n
"""
import sys
from math import comb

import numpy as np

N, K = int(sys.argv[1]), int(sys.argv[2])
BETA, MU = float(eval(sys.argv[3], {})), float(eval(sys.argv[4], {}))
n = int(sys.argv[5])

binK = np.array([comb(K, j) * BETA ** j * (1 - BETA) ** (K - j) for j in range(K + 1)])


def hyper(u, j):
    """Law of the number of elements of a uniform j-subset of [N] outside a fixed u-set."""
    out = np.zeros(j + 1)
    tot = comb(N, j)
    for r in range(j + 1):
        out[r] = comb(N - u, r) * comb(u, j - r) / tot
    return out


# one-newborn kernel on union size u -> u'
step = np.zeros((N + 1, N + 1))
for u in range(N + 1):
    for j in range(K + 1):
        h = hyper(u, j)
        for r, q in enumerate(h):
            if q > 0:
                step[u, u + r] += binK[j] * q

P = np.zeros((N + 1, N + 1))
P[0, 0] = 1.0
row = np.zeros(N + 1)
row[0] = 1.0
for y in range(1, N + 1):
    row = row @ step
    P[y] = row
if not (np.allclose(P.sum(axis=1), 1)): raise SystemExit("CHECK FAILED: " + 'np.allclose(P.sum(axis=1), 1)')
ys = np.arange(N + 1, dtype=float)
Q = P[1:, 1:]
I = np.eye(N)
ell = np.zeros(N + 1)
ell[1:] = np.linalg.solve(I - Q, ys[1:])                      # E L from y
rhs = ys[1:] ** 2 + 2 * ys[1:] * (P[1:] @ ell)
s2 = np.zeros(N + 1)
s2[1:] = np.linalg.solve(I - Q, rhs)                          # E L^2 from y
lap = np.ones(N + 1)
lap[1:] = np.linalg.solve(I - np.exp(-MU * ys[1:])[:, None] * Q, np.exp(-MU * ys[1:]) * P[1:, 0])
m = MU * ell[n]
varL = s2[n] - ell[n] ** 2
print(f"N={N} K={K} beta={BETA} mu={MU} sample n={n} (stationary)")
print(f"E C (all versions)      = {m:.6f}")
print(f"Var C, Version 2        = {m:.6f}   (SD {np.sqrt(m):.3f})")
print(f"Var C, Version 3        = {m + MU ** 2 * varL:.6f}   (SD {np.sqrt(m + MU ** 2 * varL):.3f})")
print(f"P(C=0): V2 = {np.exp(-m):.6e}, V3 = {lap[n]:.6e}")
