"""Small exact extras for the revision.

Usage: python3 -I -B code/exact_extras.py [--out DIR] [--cmax 45]

(a) N=2, u=1, beta=2 (rho=1), mu=1/2 (eta=1/2): the stationary law of the
    whole-population count C = X+Y+Z from the forward chain on (X,Y,Z),
    truncated at X+Y+Z <= cmax (innovations blocked at the boundary), at
    b in {1/4, 1/2, 1}. Reports P(C <= q), q = 0..6, and a truncation bound
    (regenerative coupling; see bound_derivation in the JSON).
(b) b = 1: exact stationary P(C_V = 0) = E exp(-mu L_N) for the homogeneous
    model (beta=13/10, mu=1/10, u=1) at N in {2,10,20,50,100}, by
    Feynman-Kac on the birth-death ancestral chain, in rational arithmetic,
    next to the no-learning upper bound prod_i u_i/(u_i+mu_i), and the
    implied complete-loss frequency and mean nonempty episode (Corollary 2).
(c) N=2, u=1, beta=2, mu=1/2: exact Var P_1 and Var P_2 at b in {0,1/4,1/2,1}
    from the forward closure.
Writes results/exact_extras.json.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import math
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from scipy.optimize import minimize_scalar
from scipy.stats import binom

from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from forward_closure import ForwardMoments, homogeneous


# ---------------------------------------------------------------------------
# (a) truncated forward chain for N = 2

def n2_truncated(b, rho=1.0, eta=0.5, cmax=40):
    states = [(x, y, z) for x in range(cmax + 1) for y in range(cmax + 1 - x)
              for z in range(cmax + 1 - x - y)]
    idx = {s: p for p, s in enumerate(states)}
    n = len(states)
    rows, cols, vals = [], [], []
    blocked = np.zeros(n)
    pmf = [binom.pmf(np.arange(x + 1), x, b) for x in range(cmax + 1)]

    def add(p, s, rate):
        if rate > 0 and s != states[p]:
            rows.append(p); cols.append(idx[s]); vals.append(rate)

    for p, (x, y, z) in enumerate(states):
        c = x + y + z
        if c < cmax:
            add(p, (x + 1, y, z), eta)
            add(p, (x, y + 1, z), eta)
        else:
            blocked[p] = 2 * eta
        add(p, (0, y + z, 0), 1.0)          # person 1 resets
        add(p, (x + z, 0, 0), 1.0)          # person 2 resets
        enc = rho / b
        if x:
            pm = pmf[x]
            for m in range(1, x + 1):
                add(p, (x - m, y, z + m), enc * pm[m])
        if y:
            pm = pmf[y]
            for m in range(1, y + 1):
                add(p, (x, y - m, z + m), enc * pm[m])
    Q = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
    out = np.asarray(Q.sum(axis=1)).ravel()
    Q = Q - sps.diags(out)
    # stationary law: pi Q = 0, sum pi = 1. Fix pi(empty)=1 and solve the
    # remaining equations (no dense normalisation row, so no LU fill-in).
    QT = Q.T.tocsc()
    keep = np.arange(1, n)
    pi_rest = spla.spsolve(QT[keep][:, keep].tocsc(), -QT[keep][:, [0]].toarray().ravel())
    pi = np.r_[1.0, pi_rest]
    pi /= pi.sum()
    C = np.array([sum(s) for s in states])
    law = np.bincount(C, weights=pi, minlength=cmax + 1)
    # truncation bound: auxiliary chain on nonempty states, absorbed at empty or
    # at the first blocked innovation (event F).
    ne = [p for p, s in enumerate(states) if sum(s) > 0]
    sub = Q[ne][:, ne].tocsc()
    A = -(sub - sps.diags(blocked[ne]))
    hF = spla.spsolve(A, blocked[ne])                  # P_s(F)
    gF = spla.spsolve(A, hF)                           # E_s[T 1_F]
    pos = {p: q for q, p in enumerate(ne)}
    start = [pos[idx[(1, 0, 0)]], pos[idx[(0, 1, 0)]]]
    PF = float(np.mean(hF[start]))
    ETF = float(np.mean(gF[start]))
    # uniform bound K on expected time to empty from any state (block argument)
    rate_quiet = 2 * eta + 2 * rho / b
    f = lambda t: t * math.exp(rate_quiet * t) / (1 - math.exp(-t)) ** 2
    opt = minimize_scalar(f, bounds=(1e-3, 5), method="bounded")
    K = float(opt.fun)
    alpha = 2 * eta
    delta = ETF + K * PF
    bound = 2 * alpha * delta
    mean = float(np.dot(np.arange(cmax + 1), law))
    var = float(np.dot(np.arange(cmax + 1) ** 2, law) - mean ** 2)
    return dict(b=b, cmax=cmax, states=n, law=law.tolist(), cdf=np.cumsum(law).tolist(),
                mean=mean, variance=var, P_F=PF, E_T_1F=ETF, K_uniform_time_to_empty=K,
                block_length=float(opt.x), truncation_bound=bound,
                boundary_mass=float(law[cmax]))


# ---------------------------------------------------------------------------
# (b) Feynman-Kac P(C_V = 0) at b = 1

def p_empty_b1(N, beta, mu, u=Fr(1)):
    """phi(n) = E_n exp(-mu int |B_t| dt); returns phi(N) exactly (Thomas algorithm)."""
    beta, mu, u = Fr(beta), Fr(mu), Fr(u)
    lam = [beta * n * (N - n) / N for n in range(N + 1)]
    dth = [u * n for n in range(N + 1)]
    # rows n=1..N: (lam+dth+mu n) phi_n - lam phi_{n+1} - dth phi_{n-1} = 0, phi_0 = 1
    diag = [lam[n] + dth[n] + mu * n for n in range(1, N + 1)]
    upper = [-lam[n] for n in range(1, N + 1)]
    lower = [-dth[n] for n in range(1, N + 1)]
    rhs = [Fr(0)] * N
    rhs[0] = dth[1]           # moves phi_0 = 1 to the right-hand side
    c, d = [Fr(0)] * N, [Fr(0)] * N
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for k in range(1, N):
        den = diag[k] - lower[k] * c[k - 1]
        c[k] = upper[k] / den if k < N - 1 else Fr(0)
        d[k] = (rhs[k] - lower[k] * d[k - 1]) / den
    phi = [Fr(0)] * N
    phi[N - 1] = d[N - 1]
    for k in range(N - 2, -1, -1):
        phi[k] = d[k] - c[k] * phi[k + 1]
    return phi[N - 1], phi


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="output directory (default: results/)")
    parser.add_argument("--cmax", type=int, default=45)
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    target = check_writable(out / "exact_extras.json", DEFAULT_RESULTS)
    ck = Checker()
    t0 = time.time()
    report = {"description": __doc__}

    # (a)
    exact_var = {0.25: 11 / 4, 0.5: 23 / 8, 1.0: 13 / 4}
    part_a = []
    for b in [0.25, 0.5, 1.0]:
        r = n2_truncated(b, cmax=args.cmax)
        r_small = n2_truncated(b, cmax=args.cmax - 10)
        r["cdf_change_vs_cmax_minus_10"] = float(np.max(np.abs(np.array(r["cdf"][:7]) - np.array(r_small["cdf"][:7]))))
        ck.check(r["truncation_bound"] < 1e-7, f"(a) truncation bound small at b={b}", bound=r["truncation_bound"])
        ck.close(f"(a) truncated mean = 3/2 at b={b}", r["mean"], 1.5, 1e-7)
        ck.close(f"(a) truncated variance = exact at b={b}", r["variance"], exact_var[b], 1e-6)
        part_a.append(r)
        print(f"(a) b={b}: P(C<=q), q=0..6:", [f"{x:.6f}" for x in r["cdf"][:7]],
              f"bound={r['truncation_bound']:.2e}", flush=True)
    ck.close("(a) P(C=0) at b=1 equals 4/11", part_a[2]["law"][0], 4 / 11, 1e-7)
    cdf = {r["b"]: r["cdf"][:7] for r in part_a}
    pattern = []
    for q in range(7):
        seq = [cdf[b][q] for b in [0.25, 0.5, 1.0]]
        if seq[0] < seq[1] < seq[2]:
            trend = "increasing in b"
        elif seq[0] > seq[1] > seq[2]:
            trend = "decreasing in b"
        else:
            trend = "not monotone in b"
        pattern.append(dict(q=q, cdf=seq, trend=trend))
    ck.check(any(p["trend"] != pattern[0]["trend"] for p in pattern),
             "(a) P(C<=q) is not monotone in b uniformly over thresholds q")
    report["a_n2_stationary_law"] = dict(
        parameters="N=2, u=1, beta=2 (rho=1), mu=1/2 (eta=1/2)", runs=part_a, cdf_pattern=pattern,
        bound_derivation=(
            "Both the true and the truncated chain regenerate at the empty state and coincide "
            "until the first innovation that would exceed cmax (event F, time T*). With delta = "
            "E[T* 1_F] + K P(F), where K bounds the expected time to empty from any state (a block "
            "of length t with no innovation or encounter and a reset of both sites, probability "
            "exp(-(2 eta + 2 rho/b) t)(1-exp(-t))^2), every stationary probability differs by at "
            "most 2 delta / E[cycle] <= 2 alpha delta, alpha = 2 eta. P(F), E[T* 1_F] are solved "
            "on the truncated chain in floating point."))

    # (b)
    part_b = []
    phi2, _ = p_empty_b1(2, 2, Fr(1, 2))
    ck.check(phi2 == Fr(4, 11), "(b) Feynman-Kac reproduces P(C_V=0)=4/11 for N=2, beta=2, mu=1/2")
    for N in [2, 10, 20, 50, 100]:
        p0, _ = p_empty_b1(N, Fr(13, 10), Fr(1, 10))
        bound = Fr(10, 11) ** N
        alpha = N * Fr(1, 10)
        ck.check(0 < p0 <= bound, f"(b) N={N}: 0 < P(C_V=0) <= prod u/(u+mu)")
        rec = dict(N=N, p_empty=float(p0), p_empty_exact=str(p0) if N <= 20 else None,
                   no_learning_bound=float(bound), ratio_to_bound=float(p0 / bound),
                   loss_episode_frequency=float(alpha * p0),
                   mean_nonempty_episode=float((1 - p0) / (alpha * p0)),
                   log10_p_empty=math.log10(p0.numerator) - math.log10(p0.denominator))
        part_b.append(rec)
        print(f"(b) N={N}: P(C_V=0)={float(p0):.6e} bound={float(bound):.6e} "
              f"ratio={float(p0 / bound):.4g} mean episode={rec['mean_nonempty_episode']:.4g}", flush=True)
    report["b_p_empty_at_b1"] = dict(parameters="homogeneous, u=1, beta=13/10, mu=1/10, b=1",
                                     method="phi(n)=E_n exp(-mu int n_t dt) on the single ancestral birth-death chain, exact rationals",
                                     runs=part_b,
                                     check_N2_beta2_mu_half=str(phi2))

    # (c)
    u, m, a = homogeneous(2, Fr(2), 1, Fr(1, 2))
    part_c = []
    for b in [Fr(0), Fr(1, 4), Fr(1, 2), Fr(1)]:
        fm = ForwardMoments(2, u, m, a, b)
        (m1, v1), (m2, v2) = fm.exactly(1), fm.exactly(2)
        cm, cv = fm.union([0, 1])
        part_c.append(dict(b=str(b), E_P1=str(m1), Var_P1=str(v1), E_P2=str(m2), Var_P2=str(v2),
                           Var_C=str(cv), Var_P1_float=float(v1), Var_P2_float=float(v2)))
        print(f"(c) b={b}: E P1={m1} Var P1={v1} ({float(v1):.6f}); E P2={m2} Var P2={v2} ({float(v2):.6f}); Var C={cv}",
              flush=True)
    ck.check([Fr(x["Var_C"]) for x in part_c] == [Fr(53, 20), Fr(11, 4), Fr(23, 8), Fr(13, 4)],
             "(c) forward closure reproduces the paper's N=2 variance table")
    ck.check(all(Fr(x["Var_P1"]) < Fr(y["Var_P1"]) for x, y in zip(part_c, part_c[1:])), "(c) Var P1 increasing")
    ck.check(all(Fr(x["Var_P2"]) < Fr(y["Var_P2"]) for x, y in zip(part_c, part_c[1:])), "(c) Var P2 increasing")
    report["c_n2_popularity_variances"] = dict(parameters="N=2, u=1, beta=2, mu=1/2", runs=part_c)

    report.update(ck.summary())
    report["seconds"] = time.time() - t0
    report["optimize_flag"] = sys.flags.optimize
    target.write_text(json.dumps(report, indent=1))
    print(target)
    ck.finish()


if __name__ == "__main__":
    main()
