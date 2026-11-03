"""Forward Monte Carlo of the three discrete-generation versions with explicit trait identities.

Independent of the theorem: traits are created and copied one by one, forward in time.
Generations g = 0..T; generation 0 learns nothing (the previous generation is empty); in every
generation each newborn first learns from its K potential models (uniform without replacement
from the previous generation), then innovates. Observation at generation T, sample A.

Version 1: G shared by all traits, success marks per trait and link.
Version 2: G and success marks per trait.
Version 3: G and success marks shared by all traits.
Innovation: 'poisson' (Poisson(mu) per newborn) or 'bernoulli' (one trait w.p. mu), the latter
to probe the role of the Poisson assumption.

Usage: mc_versions.py N K beta mu T reps innovation seed
"""
import itertools
import sys

import numpy as np

N, K = int(sys.argv[1]), int(sys.argv[2])
BETA, MU = float(eval(sys.argv[3], {})), float(eval(sys.argv[4], {}))
T, REPS = int(sys.argv[5]), int(sys.argv[6])
INNOV = sys.argv[7]
SEED = int(sys.argv[8])
CHUNK = 20_000

SUBS = np.array([sum(1 << k for k in c) for c in itertools.combinations(range(N), K)], dtype=np.int64)
BITS = 1 << np.arange(N, dtype=np.int64)
SAMPLES = {"A={0}": 0b1, "A={0,1}": 0b11, "A=all": (1 << N) - 1}


def popcount(x):
    c = np.zeros_like(x)
    for k in range(N):
        c += (x >> k) & 1
    return c


def draw_innov(rng, shape):
    if INNOV == "poisson":
        return rng.poisson(MU, size=shape)
    return (rng.random(shape) < MU).astype(np.int64)


def run_chunk(rng, version, R):
    X = N * (T + 1)                       # innovators (generation g, individual j)
    Pi = draw_innov(rng, (R, X))
    S = max(int(Pi.max()), 1)
    alive = np.arange(S)[None, None, :] < Pi[:, :, None]          # (R, X, S) trait exists
    gen_of = np.repeat(np.arange(T + 1), N)                         # innovator generation
    who = np.tile(np.arange(N), T + 1)                              # innovator individual
    pres = np.zeros((R, X, S), dtype=np.int64)                      # carrier bitmask per trait
    for g in range(T + 1):
        if g > 0:
            new = np.zeros_like(pres)
            if version == 3:
                G = SUBS[rng.integers(len(SUBS), size=(R, N))]      # (R, newborn)
                for i in range(N):
                    succ = np.zeros(R, dtype=np.int64)
                    for k in range(N):
                        ok = ((G[:, i] >> k) & 1).astype(bool) & (rng.random(R) < BETA)
                        succ |= np.where(ok, 1 << k, 0)
                    hit = (pres & succ[:, None, None]) != 0
                    new |= np.where(hit, 1 << i, 0)
            else:
                if version == 1:
                    G = SUBS[rng.integers(len(SUBS), size=(R, N))][:, None, None, :]
                    G = np.broadcast_to(G, (R, X, S, N))
                else:
                    G = SUBS[rng.integers(len(SUBS), size=(R, X, S, N))]
                for i in range(N):
                    succ = np.zeros((R, X, S), dtype=np.int64)
                    for k in range(N):
                        ok = ((G[..., i] >> k) & 1).astype(bool) & (rng.random((R, X, S)) < BETA)
                        succ |= np.where(ok, 1 << k, 0)
                    hit = (pres & succ) != 0
                    new |= np.where(hit, 1 << i, 0)
            pres = new
        born = (gen_of == g)
        pres = pres | np.where(born[None, :, None] & alive, (1 << who)[None, :, None], 0)
    out = {}
    for name, A in SAMPLES.items():
        cnt = popcount(pres & A)
        cnt = np.where(alive, cnt, 0)
        n = bin(A).count("1")
        C = (cnt > 0).sum(axis=(1, 2))
        Pk = [(cnt == k).sum(axis=(1, 2)) for k in range(1, n + 1)]
        out[name] = (C, Pk)
    return out


def summarise(x):
    x = x.astype(np.float64)
    m = x.mean()
    v = x.var(ddof=1)
    m4 = ((x - m) ** 4).mean()
    se_m = np.sqrt(v / len(x))
    se_v = np.sqrt(max(m4 - v * v, 0) / len(x))
    return m, se_m, v, se_v


def main():
    print(f"N={N} K={K} beta={BETA} mu={MU} T={T} reps={REPS} innovation={INNOV} seed={SEED}")
    rng = np.random.default_rng(SEED)
    res = {}
    for version in (2, 1, 3):
        acc = {}
        done = 0
        while done < REPS:
            R = min(CHUNK, REPS - done)
            o = run_chunk(rng, version, R)
            for name, (C, Pk) in o.items():
                a = acc.setdefault(name, [[], [[] for _ in Pk]])
                a[0].append(C)
                for j, p in enumerate(Pk):
                    a[1][j].append(p)
            done += R
        res[version] = {nm: (np.concatenate(a[0]), [np.concatenate(p) for p in a[1]])
                        for nm, a in acc.items()}
    for name in SAMPLES:
        print(f"\n{name}")
        for label, getter in [("C", lambda r: r[0])] + [
                (f"P_{k + 1}", (lambda kk: (lambda r: r[1][kk]))(k)) for k in range(bin(SAMPLES[name]).count('1'))]:
            line = f"  {label:>4}:"
            for v in (2, 1, 3):
                m, sm, var, sv = summarise(getter(res[v][name]))
                line += f"  V{v}: mean {m:.4f}({sm:.4f}) var {var:.4f}({sv:.4f})"
            print(line)
        p0 = [np.mean(res[v][name][0] == 0) for v in (2, 1, 3)]
        se0 = [np.sqrt(p * (1 - p) / REPS) for p in p0]
        print("  P(C=0): " + "  ".join(f"V{v}: {p:.5f}({s:.5f})" for v, p, s in zip((2, 1, 3), p0, se0)))


if __name__ == "__main__":
    main()
