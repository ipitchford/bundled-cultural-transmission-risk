"""Exact moments of the repertoire count in the three discrete-generation versions.

Model (Kobayashi, Wakano & Ohtsuki 2018, as described by Wakano & Aoki 2024, Sec. 1.3):
N individuals per non-overlapping generation; each newborn has K potential models drawn
uniformly without replacement from the previous generation; each link is successful with
probability beta; each newborn innovates Poisson(mu) novel traits.  A fixed sample A of the
observation generation 0 is examined; innovations occur in generations -T..0 (generation
-T-1 is empty), or T = infinity (stationary).

Version 1: potential models shared across traits, success independent per trait.
Version 2: potential models and successes drawn independently per trait.
Version 3: potential models and successful links shared by all traits.

Two independent computations are made.
 (i) A replica ("pair") dynamic programme: with L the ancestral volume sum_s |A_s|,
     E C = mu E L and Var C = mu E L + mu^2 (E[L L'] - (E L)^2), where L, L' are volumes of two
     conditionally independent per-trait genealogies given the shared randomness of the version.
     Exact rational arithmetic (fractions).  Also the stationary case by solving linear systems.
 (ii) Brute-force enumeration of every potential-model assignment G over T generations, giving
     the full law of the Version-1 conditional intensity M1 = mu E[L | G], hence Var(M1) and
     P(C = 0) = E exp(-M1) (floating point); Version 3 by enumerating H.
"""
import itertools
import math
import sys
from fractions import Fraction as Fr

N, K = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (3, 2)
BETA = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 2)
MU = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(1, 2)

SUBSETS_K = [sum(1 << k for k in c) for c in itertools.combinations(range(N), K)]


def popcount(x):
    return bin(x).count("1")


def thin_dist(g, beta):
    """Law of the successful subset of potential-model set g (bitmask)."""
    out = {}
    elems = [k for k in range(N) if g >> k & 1]
    for marks in itertools.product([0, 1], repeat=len(elems)):
        p = Fr(1)
        s = 0
        for k, m in zip(elems, marks):
            p *= beta if m else 1 - beta
            if m:
                s |= 1 << k
        out[s] = out.get(s, 0) + p
    return out


THIN = {g: thin_dist(g, BETA) for g in SUBSETS_K}
PG = Fr(1, len(SUBSETS_K))


def newborn_pair_law(in1, in2, version):
    """Joint law of (S1, S2): contributions of one newborn to the two replicas' parent sets."""
    out = {}
    if version == 2:
        d1 = {0: Fr(1)}
        d2 = {0: Fr(1)}
        if in1:
            d1 = {}
            for g in SUBSETS_K:
                for s, p in THIN[g].items():
                    d1[s] = d1.get(s, 0) + PG * p
        if in2:
            d2 = {}
            for g in SUBSETS_K:
                for s, p in THIN[g].items():
                    d2[s] = d2.get(s, 0) + PG * p
        for s1, p1 in d1.items():
            for s2, p2 in d2.items():
                out[(s1, s2)] = out.get((s1, s2), 0) + p1 * p2
        return out
    for g in SUBSETS_K:
        t = THIN[g]
        if version == 3:
            for s, p in t.items():
                key = (s if in1 else 0, s if in2 else 0)
                out[key] = out.get(key, 0) + PG * p
        else:  # version 1
            for s1, p1 in (t.items() if in1 else [(0, Fr(1))]):
                for s2, p2 in (t.items() if in2 else [(0, Fr(1))]):
                    out[(s1, s2)] = out.get((s1, s2), 0) + PG * p1 * p2
    return out


def pair_transition(b1, b2, version):
    dist = {(0, 0): Fr(1)}
    for i in range(N):
        in1, in2 = b1 >> i & 1, b2 >> i & 1
        if not (in1 or in2):
            continue
        law = newborn_pair_law(in1, in2, version)
        new = {}
        for (x1, x2), p in dist.items():
            for (s1, s2), q in law.items():
                key = (x1 | s1, x2 | s2)
                new[key] = new.get(key, 0) + p * q
        dist = new
    return dist


def single_transition(b):
    d = pair_transition(b, 0, 1)
    out = {}
    for (x1, _), p in d.items():
        out[x1] = out.get(x1, 0) + p
    return out


ALL = range(1 << N)
SINGLE = {b: single_transition(b) for b in ALL}
PAIR = {v: {(b1, b2): pair_transition(b1, b2, v) for b1 in ALL for b2 in ALL} for v in (1, 2, 3)}


def finite_horizon(A, T, version):
    """Return (E L, E[L L']) for depth T (innovations in generations -T..0)."""
    ell = {b: Fr(popcount(b)) for b in ALL}  # depth 0 remaining
    F = {(b1, b2): Fr(popcount(b1) * popcount(b2)) for b1 in ALL for b2 in ALL}
    for _ in range(T):
        ell_new = {}
        for b in ALL:
            ell_new[b] = popcount(b) + sum(p * ell[x] for x, p in SINGLE[b].items())
        F_new = {}
        for (b1, b2), tr in PAIR[version].items():
            acc = Fr(popcount(b1) * popcount(b2))
            for (x1, x2), p in tr.items():
                acc += p * (popcount(b1) * ell[x2] + popcount(b2) * ell[x1] + F[(x1, x2)])
            F_new[(b1, b2)] = acc
        ell, F = ell_new, F_new
    return ell[A], F[(A, A)]


def solve_linear(states, coef, rhs):
    """Solve x_s = rhs_s + sum_t coef[s][t] x_t exactly (Gaussian elimination)."""
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    M = [[Fr(0)] * (n + 1) for _ in range(n)]
    for s in states:
        i = idx[s]
        M[i][i] += 1
        for t, c in coef[s].items():
            if t in idx:
                M[i][idx[t]] -= c
        M[i][n] = rhs[s]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return {s: M[idx[s]][n] for s in states}


def stationary(A, version):
    states = [b for b in ALL if b]
    ell = solve_linear(states, {b: SINGLE[b] for b in states}, {b: Fr(popcount(b)) for b in states})
    ell[0] = Fr(0)
    pstates = [(b1, b2) for b1 in ALL for b2 in ALL if b1 and b2]
    coef, rhs = {}, {}
    for s in pstates:
        b1, b2 = s
        tr = PAIR[version][s]
        coef[s] = tr
        rhs[s] = popcount(b1) * popcount(b2) + sum(
            p * (popcount(b1) * ell[x2] + popcount(b2) * ell[x1]) for (x1, x2), p in tr.items())
    F = solve_linear(pstates, coef, rhs)
    return ell[A], F[(A, A)]


def moments(EL, ELL):
    mean = MU * EL
    var = MU * EL + MU * MU * (ELL - EL * EL)
    return mean, var


# ---------------- brute-force enumeration (independent route) ----------------

def v1_enumeration(A, T):
    """Enumerate all potential-model assignments of all N newborns in generations 0..-(T-1)."""
    thinf = {g: {s: float(p) for s, p in THIN[g].items()} for g in SUBSETS_K}
    mu = float(MU)
    vals = []
    nconf = len(SUBSETS_K) ** (N * T)
    for conf in itertools.product(SUBSETS_K, repeat=N * T):
        dist = {A: 1.0}
        EL = float(popcount(A))
        for t in range(T):
            gs = conf[t * N:(t + 1) * N]
            new = {}
            for b, p in dist.items():
                d = {0: p}
                for i in range(N):
                    if b >> i & 1:
                        nd = {}
                        for x, q in d.items():
                            for s, r in thinf[gs[i]].items():
                                nd[x | s] = nd.get(x | s, 0.0) + q * r
                        d = nd
                for x, q in d.items():
                    new[x] = new.get(x, 0.0) + q
            dist = new
            EL += sum(q * popcount(x) for x, q in dist.items())
        vals.append(mu * EL)
    w = 1.0 / nconf
    m = sum(vals) * w
    var_m = sum((v - m) ** 2 for v in vals) * w
    p0 = sum(math.exp(-v) for v in vals) * w
    return m, m + var_m, p0


def v3_law_L(A, T):
    """Exact law of the shared ancestral volume L (Version 3) by DP over (set, volume)."""
    hlaw = {}
    for g in SUBSETS_K:
        for s, p in THIN[g].items():
            hlaw[s] = hlaw.get(s, 0) + PG * p
    dist = {(A, popcount(A)): Fr(1)}
    for _ in range(T):
        new = {}
        for (b, vol), p in dist.items():
            d = {0: p}
            for i in range(N):
                if b >> i & 1:
                    nd = {}
                    for x, q in d.items():
                        for s, r in hlaw.items():
                            nd[x | s] = nd.get(x | s, 0) + q * r
                    d = nd
            for x, q in d.items():
                key = (x, vol + popcount(x))
                new[key] = new.get(key, 0) + q
        dist = new
    law = {}
    for (_, vol), p in dist.items():
        law[vol] = law.get(vol, 0) + p
    return law


def main():
    print(f"N={N}, K={K}, beta={BETA}, mu={MU}; innovation Poisson(mu) per newborn")
    print("Exact (rational) mean and variance of the sample repertoire count C, by version\n")
    samples = {"A={0} (n=1)": 0b001, "A={0,1} (n=2)": 0b011, "A=all (n=3)": 0b111}
    for name, A in samples.items():
        print(name)
        for T in [0, 1, 2, 3, 4, "stationary"]:
            row = []
            for v in (2, 1, 3):
                EL, ELL = stationary(A, v) if T == "stationary" else finite_horizon(A, T, v)
                row.append(moments(EL, ELL))
            means = {r[0] for r in row}
            if not (len(means) == 1): raise SystemExit("CHECK FAILED: " + "means differ")
            v2, v1, v3 = (r[1] for r in row)
            if not (v2 <= v1 <= v3): raise SystemExit("CHECK FAILED: " + 'v2 <= v1 <= v3')
            if not (v2 == row[0][0]): raise SystemExit("CHECK FAILED: " + "Version 2 variance must equal its mean (Poisson)")
            print(f"  T={T!s:>10}: mean={float(row[0][0]):.6f}  Var V2={float(v2):.6f}  "
                  f"V1={float(v1):.6f}  V3={float(v3):.6f}  "
                  f"[V1-V2={v1 - v2}, V3-V1={v3 - v1}]")
        print()
    print("Brute-force enumeration check (floating point) of Version 1 and Version 3")
    for name, A in samples.items():
        for T in [1, 2, 3]:
            m1, var1, p01 = v1_enumeration(A, T)
            EL, ELL = finite_horizon(A, T, 1)
            mean_dp, var_dp = moments(EL, ELL)
            law = v3_law_L(A, T)
            muf = float(MU)
            EL3 = sum(float(p) * v for v, p in law.items())
            EL3sq = sum(float(p) * v * v for v, p in law.items())
            var3 = muf * EL3 + muf ** 2 * (EL3sq - EL3 ** 2)
            p03 = sum(float(p) * math.exp(-muf * v) for v, p in law.items())
            p02 = math.exp(-float(mean_dp))
            _, var3_dp = moments(*finite_horizon(A, T, 3))
            print(f"  {name} T={T}: V1 var enum={var1:.10f} pair-DP={float(var_dp):.10f} | "
                  f"V3 var enum={var3:.10f} pair-DP={float(var3_dp):.10f} | "
                  f"P(C=0): V2={p02:.6f} V1={p01:.6f} V3={p03:.6f}")
            if not (abs(var1 - float(var_dp)) < 1e-9 and abs(var3 - float(var3_dp)) < 1e-9): raise SystemExit("CHECK FAILED: " + 'abs(var1 - float(var_dp)) < 1e-9 and abs(var3 - float(var3_dp)) < 1e-9')
            if not (abs(m1 - float(mean_dp)) < 1e-12): raise SystemExit("CHECK FAILED: " + 'abs(m1 - float(mean_dp)) < 1e-12')
            if not (p02 <= p01 + 1e-15 <= p03 + 2e-15): raise SystemExit("CHECK FAILED: " + 'p02 <= p01 + 1e-15 <= p03 + 2e-15')
    print("\nAll checks passed.")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
