"""Scaling study of the whole-population variance V_N(b) = Var(C_V,b).

Usage:
    python3 -I -B code/scaling_study.py [--out DIR] [--N 10,20,...] [--beta 0.7,...]
                                        [--no-direct]

Solves the exact lumped pair-replica system of Section 6.1 (homogeneous
complete graph, u=1, mu=0.1) for every (beta, N, b) on the grid below, in its
cancellation-free covariance form (see PairSystem): Jacobi-preconditioned
BiCGSTAB inner solves, warm-started along b, inside iterative refinement whose
residuals are computed in double-double arithmetic from exact rational rates.
For N <= 50 every solve is repeated with a direct sparse LU (scipy spsolve);
for N <= 30 also with the original h-formulation of Section 6.1. At b=1 the result is compared with the
independent single-genealogy second-moment recursion. For each solve the
script reports the relative residual and an a-posteriori bound

    |c - c_true|_inf <= T_max * |A c - Gamma|_inf,

where A = -Q is the nonsingular M-matrix of the pair chain and T_max, the
largest expected one-replica extinction time, bounds |A^{-1}|_inf (expected
exit time of the pair chain). The bound is evaluated in floating point.
The exact derivative V'(b) is obtained from one further solve,
dc/db = A^{-1} (dGamma/db - (dA/db) c), giving the sensitivity V'(b)/V(b).

Writes results/scaling.json. figures are drawn by code/plot_scaling.py.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from fractions import Fraction as Fr

import ddarith as dd
from transmission import single_moments

B_GRID = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0]
B_EXTRA_13 = [13 / 14]          # Wakano et al. (2026) Fig. 16 setting, added for beta = 1.3
BETAS = [0.7, 1.0, 1.3, 2.0]
NS = [10, 20, 30, 50, 75, 100, 150]
# (N, beta) pairs skipped by default: at beta/u = 2, N = 150 the pair chain's
# expected exit time is ~1.6e12, the double-precision Krylov inner solves do
# not converge (an exploratory run diverged at b = 0.75), so no value is reported.
SKIP = {(150, 2.0)}
U, MU = 1.0, 0.1


class LumpedPair:
    """Structure matrices of the lumped pair generator: A(b) = a (E0 + b E1) + u Dm."""

    def __init__(self, N):
        self.N = N
        st = [(i, j, k) for k in range(N + 1) for i in range(N - k + 1)
              for j in range(N - k - i + 1) if i + k > 0 and j + k > 0]
        S = np.array(st, dtype=np.int64)
        I, J, K = S[:, 0], S[:, 1], S[:, 2]
        n = len(st)
        lookup = -np.ones((N + 2, N + 2, N + 2), dtype=np.int64)
        lookup[I, J, K] = np.arange(n)
        self.states, self.lookup = S, lookup
        D = N - I - J - K
        rows = np.arange(n)

        def offdiag(di, dj, dk, rate):
            ii, jj, kk = I + di, J + dj, K + dk
            ok = (rate != 0) & (ii >= 0) & (jj >= 0) & (kk >= 0)
            dest = np.full(n, -1)
            dest[ok] = lookup[ii[ok], jj[ok], kk[ok]]
            keep = dest >= 0
            return rows[keep], dest[keep], -rate[keep].astype(float)

        def build(parts, diag):
            r = np.concatenate([p[0] for p in parts] + [rows])
            c = np.concatenate([p[1] for p in parts] + [rows])
            v = np.concatenate([p[2] for p in parts] + [diag.astype(float)])
            return sps.csr_matrix((v, (r, c)), shape=(n, n))

        # encounter rates per unit a, at b = 0 and the coefficient of b
        r_i0, r_j0 = D * (I + K), D * (J + K)
        r_ii, r_jj = I * (J + K), J * (I + K)
        self.E0 = build([offdiag(1, 0, 0, r_i0), offdiag(0, 1, 0, r_j0),
                         offdiag(-1, 0, 1, r_ii), offdiag(0, -1, 1, r_jj)],
                        r_i0 + r_j0 + r_ii + r_jj)
        dk = D * K
        self.E1 = build([offdiag(1, 0, 0, -dk), offdiag(0, 1, 0, -dk), offdiag(0, 0, 1, dk)],
                        -dk)   # net out-rate change: -dk - dk + dk
        self.Dm = build([offdiag(-1, 0, 0, I), offdiag(0, -1, 0, J), offdiag(0, 0, -1, K)], I + J + K)
        self.target = int(lookup[0, 0, N])

    def matrix(self, beta, b, u=U):
        a = beta / self.N
        return (a * (self.E0 + b * self.E1) + u * self.Dm).tocsr()

    def rhs(self, ell):
        I, J, K = self.states.T
        return (I + K) * ell[J + K] + (J + K) * ell[I + K]


def max_extinction_time(N, beta, u=U):
    """max_n E_n[extinction time] of the single ancestral birth-death chain (= start at N)."""
    inc = np.zeros(N + 2)
    for k in range(N, 0, -1):
        inc[k] = 1 / (u * k) + beta * (N - k) / (N * u) * inc[k + 1]
    # E_k[time to reach k-1] = 1/(u k) + (lambda_k/(u k)) E_{k+1}[time to reach k]
    return float(inc[1:N + 1].sum())


def bicg(A, rhs, x0=None, rtol=1e-12, maxiter=20000):
    dinv = 1.0 / A.diagonal()
    M = spla.LinearOperator(A.shape, lambda x: dinv * x)
    it = [0]

    def cb(_):
        it[0] += 1
    scale = float(np.max(np.abs(rhs))) or 1.0
    x, info = spla.bicgstab(A, rhs / scale, x0=None if x0 is None else x0 / scale, M=M,
                            rtol=rtol, atol=0, maxiter=maxiter, callback=cb)
    return x * scale, int(info), it[0]


class PairSystem:
    """Exact data of the cancellation-free covariance system for one (N, beta).

    c(i,j,k) = Cov(R_1, R_2)/mu^2 for the two replicas solves A(b) c = Gamma(b),
    where Gamma collects the joint jumps of the single-set mean ell:
        Gamma = u k D_{n1} D_{n2} + (beta/N) b (N-i-j-k) k D_{n1+1} D_{n2+1},
    n1 = i+k, n2 = j+k, D_n = ell_n - ell_{n-1} (D_{N+1} = 0). Then
        Var(C_V) = mu ell_N + mu^2 c(0,0,N).
    (Derivation: apply the pair generator to ell(n1) ell(n2); the cross term
    is the carre-du-champ Gamma, so F - ell x ell solves the same system with
    right-hand side Gamma >= 0. No subtraction of large numbers is needed.)
    """

    def __init__(self, N, beta):
        self.N = N
        self.betaF = Fr(beta).limit_denominator(1000)
        self.beta = float(self.betaF)
        self.lp = LumpedPair(N)
        D = [Fr(0)] * (N + 2)
        for n in range(N, 0, -1):
            D[n] = 1 + self.betaF * (N - n) / N * D[n + 1]     # u = 1
        self.ellN = sum(D[1:N + 1])
        I, J, K = self.lp.states.T
        n1, n2, Dd = I + K, J + K, N - I - J - K
        Ph, Pl = np.zeros((N + 2, N + 2)), np.zeros((N + 2, N + 2))
        for p in range(N + 2):
            for q in range(p, N + 2):
                h, l = dd.from_fraction(D[p] * D[q])
                Ph[p, q] = Ph[q, p] = h
                Pl[p, q] = Pl[q, p] = l
        self.G0 = dd.mul_d(Ph[n1, n2], Pl[n1, n2], K.astype(float))
        m1, m2 = np.minimum(n1 + 1, N + 1), np.minimum(n2 + 1, N + 1)
        self.G1 = dd.mul_d(Ph[m1, m2], Pl[m1, m2], (Dd * K).astype(float))   # times a*b
        self.a = dd.from_fraction(self.betaF / N)

    def Ax(self, b, xh, xl):
        bb = dd.from_fraction(Fr(b).limit_denominator(10 ** 6))
        e0 = dd.matvec(self.lp.E0, xh, xl)
        e1 = dd.matvec(self.lp.E1, xh, xl)
        dm = dd.matvec(self.lp.Dm, xh, xl)
        t = dd.add(*e0, *dd.mul_dd(*e1, *bb))
        return dd.add(*dd.mul_dd(*t, *self.a), *dm)

    def gamma(self, b):
        bb = dd.from_fraction(Fr(b).limit_denominator(10 ** 6))
        return dd.add(*self.G0, *dd.mul_dd(*dd.mul_dd(*self.G1, *self.a), *bb))

    def dgamma(self):
        return dd.mul_dd(*self.G1, *self.a)

    def E1x(self, xh, xl):
        return dd.mul_dd(*dd.matvec(self.lp.E1, xh, xl), *self.a)


def refine(sysm, b, A, rhs_dd, x0=None, max_steps=8, tol=1e-17):
    """Iterative refinement with double-double residuals; Krylov inner solves."""
    t = time.time()
    n = A.shape[0]
    xh, xl = (np.zeros(n), np.zeros(n)) if x0 is None else x0
    hist, its_total, infos = [], 0, []
    scale = float(np.max(np.abs(rhs_dd[0])))
    r = None
    for step in range(max_steps):
        yh, yl = sysm.Ax(b, xh, xl)
        rh, rl = dd.add(*rhs_dd, -yh, -yl)
        r = rh + rl
        rel = float(np.max(np.abs(r)) / scale)
        hist.append(rel)
        if not np.isfinite(rel) or (step and rel < 1e-28):
            break
        e, info, its = bicg(A, r)
        its_total += its
        infos.append(info)
        if not np.all(np.isfinite(e)):
            hist.append(float("nan"))
            break
        xh, xl = dd.add(xh, xl, e, np.zeros(n))
        corr = float(np.max(np.abs(e)) / max(np.max(np.abs(xh)), 1e-300))
        if corr < tol:
            yh, yl = sysm.Ax(b, xh, xl)
            rh, rl = dd.add(*rhs_dd, -yh, -yl)
            r = rh + rl
            hist.append(float(np.max(np.abs(r)) / scale))
            break
    return (xh, xl), dict(method="Jacobi-BiCGSTAB inner solves + double-double residual refinement",
                          refinement_residual_history=hist, inner_iterations=its_total,
                          inner_info=infos, final_relative_residual=hist[-1],
                          residual_inf=float(np.max(np.abs(r))) if np.isfinite(hist[-1]) else float("nan"),
                          seconds=time.time() - t)


def run(N, beta, bs, ck, direct, h_check=False):
    t0 = time.time()
    sysm = PairSystem(N, beta)
    lp = sysm.lp
    beta = sysm.beta
    ell, ell2 = single_moments(N, beta, U)
    Tmax = max_extinction_time(N, beta)
    mean = float(Fr(1, 10) * sysm.ellN)
    v_full = MU * ell[N] + MU ** 2 * (ell2[N] - ell[N] ** 2)   # independent b=1 value
    rows = []
    x0 = None
    for b in sorted(bs):
        A = lp.matrix(beta, b)
        (ch, cl), info = refine(sysm, b, A, sysm.gamma(b), x0)
        ok = np.isfinite(info["final_relative_residual"]) and info["final_relative_residual"] < 1e-12
        if not ok:   # recorded as not computed (precision-limited), not as a reported value
            print(f"  N={N} beta={beta} b={b}: refinement FAILED {info['refinement_residual_history']}", flush=True)
            return dict(N=N, beta=beta, failed=True, states=int(lp.states.shape[0]), by_b=rows,
                        failed_at_b=b, failed_history=info["refinement_residual_history"],
                        note="double-precision Krylov inner solves did not converge; not reported")
        x0 = (ch, cl)
        rhs_d = dd.add(*sysm.dgamma(), *[-v for v in sysm.E1x(ch, cl)])
        (dh, dl), info_d = refine(sysm, b, A, rhs_d, None, max_steps=4, tol=1e-12)
        c_target = ch[lp.target] + cl[lp.target]
        V = mean + MU ** 2 * c_target
        dV = MU ** 2 * (dh[lp.target] + dl[lp.target])
        err_c = Tmax * info["residual_inf"]
        rec = dict(b=b, variance=float(V), dV_db=float(dV), sensitivity=float(dV / V),
                   variance_error_bound=float(MU ** 2 * err_c),
                   variance_relative_error_bound=float(MU ** 2 * err_c / V), solve=info,
                   derivative_solve=dict(final_relative_residual=info_d["final_relative_residual"],
                                         seconds=info_d["seconds"]))
        if direct:
            g = sysm.gamma(b)[0]
            cd = spla.spsolve(A.tocsc(), g)
            Vd = mean + MU ** 2 * cd[lp.target]
            rec["direct_variance"] = float(Vd)
            rec["direct_relative_difference"] = float(abs(Vd - V) / abs(V))
            rec["direct_relative_residual"] = float(np.max(np.abs(A @ cd - g)) / np.max(g))
            ck.check(rec["direct_relative_difference"] < 1e-10,
                     f"N={N} beta={beta} b={b}: refined iterative vs direct sparse LU",
                     diff=rec["direct_relative_difference"])
            if h_check:
                hd = spla.spsolve(A.tocsc(), lp.rhs(ell))      # Section 6.1 formulation
                Vh = mean + MU ** 2 * (hd[lp.target] - ell[N] ** 2)
                rec["h_formulation_variance"] = float(Vh)
                rec["h_formulation_relative_difference"] = float(abs(Vh - V) / V)
                ck.check(rec["h_formulation_relative_difference"] < 1e-8,
                         f"N={N} beta={beta} b={b}: covariance vs h formulation",
                         diff=rec["h_formulation_relative_difference"])
        ck.check(rec["variance_relative_error_bound"] < 1e-10, f"N={N} beta={beta} b={b}: a-posteriori bound",
                 bound=rec["variance_relative_error_bound"])
        rows.append(rec)
        print(f"  N={N:3d} beta={beta} b={b:<8.6g} V={V:.12g} V'={dV:.6g} "
              f"relres={info['final_relative_residual']:.1e} bound={rec['variance_relative_error_bound']:.1e} "
              f"{info['seconds']:.1f}s", flush=True)
    byb = {r["b"]: r for r in rows}
    V1 = byb[1.0]["variance"]
    V0 = byb[0.0]["variance"]
    ck.close(f"N={N} beta={beta}: b=1 pair system vs single-genealogy recursion", V1, v_full, 1e-9,
             relative=True)
    vs = [r["variance"] for r in rows]
    ck.check(all(x < y for x, y in zip(vs, vs[1:])), f"N={N} beta={beta}: V strictly increasing in b")
    ck.check(all(r["dV_db"] > 0 for r in rows), f"N={N} beta={beta}: V'(b) > 0")
    for r in rows:
        r["ratio_to_b1"] = r["variance"] / V1
    return dict(N=N, beta=beta, u=U, mu=MU, states=int(lp.states.shape[0]), mean=float(mean),
                genealogy_mean=float(sysm.ellN), variance_b1_single_genealogy=float(v_full),
                b1_relative_difference=abs(V1 - v_full) / v_full,
                T_max=Tmax, seconds=time.time() - t0,
                V0=V0, V1=V1, V1_over_V0=V1 / V0,
                endpoint_slope_fd=(V1 - byb[0.99]["variance"]) / 0.01,
                endpoint_slope_exact=byb[1.0]["dV_db"],
                fraction_of_gain_above_0_9=(V1 - byb[0.9]["variance"]) / (V1 - V0),
                fraction_of_gain_above_0_99=(V1 - byb[0.99]["variance"]) / (V1 - V0),
                fraction_of_V1_lost_at_0_99=(V1 - byb[0.99]["variance"]) / V1,
                by_b=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="output directory (default: results/)")
    parser.add_argument("--N", default=",".join(map(str, NS)))
    parser.add_argument("--beta", default=",".join(map(str, BETAS)))
    parser.add_argument("--no-direct", action="store_true", help="skip direct cross-checks")
    parser.add_argument("--direct-max-N", type=int, default=50)
    parser.add_argument("--name", default="scaling.json")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    target = check_writable(out / args.name, DEFAULT_RESULTS)
    ck = Checker()
    t0 = time.time()
    results = []
    for beta in [float(x) for x in args.beta.split(",")]:
        for N in [int(x) for x in args.N.split(",")]:
            if (N, beta) in SKIP:
                print(f"beta={beta} N={N}: skipped (precision-infeasible in double-precision inner solves)", flush=True)
                continue
            bs = B_GRID + (B_EXTRA_13 if abs(beta - 1.3) < 1e-12 else [])
            res = run(N, beta, bs, ck, direct=(not args.no_direct and N <= args.direct_max_N),
                      h_check=(N <= 30))
            results.append(res)
            if res.get("failed"):
                continue
            print(f"beta={beta} N={N}: states={res['states']} mean={res['mean']:.6f} V0={res['V0']:.6g} "
                  f"V1={res['V1']:.6g} V1/V0={res['V1_over_V0']:.4g} slope={res['endpoint_slope_fd']:.4g} "
                  f"frac>0.9={res['fraction_of_gain_above_0_9']:.3f} ({res['seconds']:.1f}s)", flush=True)
            report = dict(description=__doc__, parameters=dict(u=U, mu=MU, b_grid=B_GRID,
                          extra_b_for_beta_1_3=B_EXTRA_13, skipped=sorted(SKIP)), results=results, status="in progress")
            target.write_text(json.dumps(report, indent=1))
    report = dict(description=__doc__, parameters=dict(u=U, mu=MU, b_grid=B_GRID,
                  extra_b_for_beta_1_3=B_EXTRA_13, skipped=sorted(SKIP)), results=results,
                  not_computed=[dict(N=r["N"], beta=r["beta"], failed_at_b=r["failed_at_b"])
                                for r in results if r.get("failed")])
    report.update(ck.summary())
    report["seconds"] = time.time() - t0
    report["optimize_flag"] = sys.flags.optimize
    report["software"] = {"python": sys.version, "numpy": np.__version__,
                          "scipy": __import__("scipy").__version__}
    target.write_text(json.dumps(report, indent=1))
    print(target)
    ck.finish()


if __name__ == "__main__":
    main()
