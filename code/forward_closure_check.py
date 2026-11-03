"""Independent forward moment-closure check of Theorems 1-2 and the pair system.

Usage: python3 -I -B code/forward_closure_check.py [--out DIR]

For each case the exact forward closure on carrier-set counts n_S
(code/forward_closure.py) gives the mean and variance of
  * C_A, the union count of the sample A;
  * P_k, the number of traits held by exactly k sites (whole population);
  * G_k, the number of traits held by at least k sites (G_1 = C_V).
C_A is compared exactly (rational equality) with the backward pair system:
the bundle's lumped generator (transmission.pair_transitions, evaluated in
Fractions) for homogeneous cases, and an independent labelled backward pair
generator for all cases. Checks: means constant in b; variances nondecreasing
in b for every functional; strict increase where Proposition 1 or the chain
example predicts it. Writes results/forward_closure_check.json.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import time
from fractions import Fraction as Fr

from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from forward_closure import ForwardMoments, backward_union, homogeneous, solve_rational
from transmission import pair_states, pair_transitions


def lumped_backward_exact(N, beta, u, mu, b, n):
    """Bundle's lumped homogeneous pair generator, solved in exact arithmetic."""
    beta, u, mu, b = Fr(beta), Fr(u), Fr(mu), Fr(b)
    a = beta / N
    # single-set expected lengths: tridiagonal birth-death, solved exactly
    rows = []
    for k in range(1, N + 1):
        birth, death = a * k * (N - k), u * k
        row = {k - 1: birth + death}
        if k < N:
            row[k] = -birth
        if k > 1:
            row[k - 2] = -death
        rows.append(row)
    ell = [Fr(0)] + solve_rational(rows, [Fr(k) for k in range(1, N + 1)])
    states = pair_states(N)
    index = {s: p for p, s in enumerate(states)}
    rows2, rhs2 = [], []
    for state in states:
        out = Fr(0)
        row = {}
        for dest, rate in pair_transitions(state, N, beta, u, b):
            out += rate
            if dest in index:
                row[index[dest]] = row.get(index[dest], Fr(0)) - rate
        row[index[state]] = row.get(index[state], Fr(0)) + out
        i, j, k = state
        rows2.append(row)
        rhs2.append((i + k) * ell[j + k] + (j + k) * ell[i + k])
    h = solve_rational(rows2, rhs2)
    hn = h[index[(0, 0, n)]]
    mean = mu * ell[n]
    return mean, mean + mu * mu * (hn - ell[n] ** 2)


def frac(x):
    return {"exact": str(x), "float": float(x)}


def run_case(ck, name, N, u, mu, a, bs, samples, homog=None, strict_union=True):
    rows = []
    for b in bs:
        fm = ForwardMoments(N, u, mu, a, b)
        entry = {"b": str(b), "union": {}, "exactly": {}, "at_least": {}}
        for A in samples:
            fmean, fvar = fm.union(A)
            bmean, bvar = backward_union(N, u, mu, a, b, A)
            ck.check(fmean == bmean and fvar == bvar,
                     f"{name}: forward == labelled backward for C_A, A={A}, b={b}",
                     forward=str(fvar), backward=str(bvar))
            rec = {"mean": frac(fmean), "variance": frac(fvar),
                   "labelled_backward_equal": fmean == bmean and fvar == bvar}
            if homog is not None:
                beta, uu, mm = homog
                lmean, lvar = lumped_backward_exact(N, beta, uu, mm, b, len(A))
                ck.check(fmean == lmean and fvar == lvar,
                         f"{name}: forward == bundle lumped backward, n={len(A)}, b={b}",
                         forward=str(fvar), lumped=str(lvar))
                rec["lumped_backward_equal"] = fmean == lmean and fvar == lvar
            entry["union"][str(A)] = rec
        for k in range(1, N + 1):
            mk, vk = fm.exactly(k)
            gk, wk = fm.at_least(k)
            entry["exactly"][k] = {"mean": frac(mk), "variance": frac(vk)}
            entry["at_least"][k] = {"mean": frac(gk), "variance": frac(wk)}
        ck.check(fm.at_least(1) == fm.union(list(range(N))), f"{name}: G_1 equals C_V, b={b}")
        rows.append(entry)

    def series(group, key):
        return [Fr(r[group][key]["mean"]["exact"]) for r in rows], \
               [Fr(r[group][key]["variance"]["exact"]) for r in rows]

    summary = {}
    groups = [("union", str(A)) for A in samples] + \
             [("exactly", k) for k in range(1, N + 1)] + [("at_least", k) for k in range(1, N + 1)]
    for group, key in groups:
        means, vars_ = series(group, key)
        const = all(x == means[0] for x in means)
        nondec = all(x <= y for x, y in zip(vars_, vars_[1:]))
        strict = all(x < y for x, y in zip(vars_, vars_[1:]))
        label = {"union": "C_A", "exactly": "P_k", "at_least": "G_k(>=k)"}[group]
        ck.check(const, f"{name}: mean of {label} [{key}] constant in b")
        ck.check(nondec, f"{name}: variance of {label} [{key}] nondecreasing in b",
                 variances=[float(v) for v in vars_])
        if group == "union" and strict_union:
            ck.check(strict, f"{name}: variance of C_A [{key}] strictly increasing in b",
                     variances=[float(v) for v in vars_])
        summary[f"{label}[{key}]"] = {
            "mean": str(means[0]), "mean_constant": const,
            "variances": [float(v) for v in vars_], "nondecreasing": nondec,
            "strictly_increasing": strict,
            "trivially_zero": all(v == 0 for v in vars_)}
    return {"name": name, "N": N, "b_grid": [str(b) for b in bs],
            "u": [str(x) for x in u], "mu": [str(x) for x in mu],
            "a_donor_to_learner": [[str(x) for x in row] for row in a],
            "samples": samples, "by_b": rows, "summary": summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="output directory (default: results/)")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    target = check_writable(out / "forward_closure_check.json", DEFAULT_RESULTS)
    ck = Checker()
    t0 = time.time()
    cases = []
    quarter = [Fr(0), Fr(1, 4), Fr(1, 2), Fr(1)]

    for N, beta, mu in [(3, Fr(13, 10), Fr(1, 5)), (4, Fr(2), Fr(1, 2))]:
        u, m, a = homogeneous(N, beta, 1, mu)
        cases.append(run_case(ck, f"homogeneous N={N} beta={beta} mu={mu} u=1", N, u, m, a, quarter,
                              [[0], list(range(N))], homog=(beta, 1, mu)))
        print(cases[-1]["name"], "done", round(time.time() - t0, 1), "s", flush=True)

    # Reviewer's spot values (methodology report W2), as an external cross-check.
    c3 = cases[0]["by_b"]
    ck.check(c3[2]["union"]["[0, 1, 2]"]["variance"]["exact"] == "401845979721590906/304364967034359375",
             "N=3 b=1/2 variance equals the methodology reviewer's value")
    ck.check(c3[3]["union"]["[0, 1, 2]"]["variance"]["exact"] == "1854829/1265625",
             "N=3 b=1 variance equals the methodology reviewer's value")
    ck.check(cases[1]["by_b"][3]["union"]["[0, 1, 2, 3]"]["variance"]["exact"] == "1361/64",
             "N=4 b=1 variance equals the methodology reviewer's value")

    # Heterogeneous directed N=3 network (methodology reviewer's example).
    u = [Fr(1), Fr(2), Fr(1, 2)]
    mu = [Fr(1), Fr(0), Fr(1, 3)]
    a = [[Fr(0), Fr(3, 2), Fr(0)], [Fr(1, 2), Fr(0), Fr(2)], [Fr(1), Fr(0), Fr(0)]]
    cases.append(run_case(ck, "heterogeneous directed N=3", 3, u, mu, a,
                          [Fr(0), Fr(1, 3), Fr(2, 3), Fr(1)], [[1], [0, 1, 2]]))
    print(cases[-1]["name"], "done", round(time.time() - t0, 1), "s", flush=True)

    # Directed chain 0 -> 1 -> 2, innovation only at 0, sample {2}. Proposition 1's
    # sufficient condition fails for A={2} (its only in-neighbour has mu=0).
    u = [Fr(1)] * 3
    mu = [Fr(1), Fr(0), Fr(0)]
    a = [[Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(1)], [Fr(0), Fr(0), Fr(0)]]
    cases.append(run_case(ck, "directed chain 0->1->2, mu only at 0, sample {2}", 3, u, mu, a,
                          [Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1)], [[2], [0, 1, 2]]))
    print(cases[-1]["name"], "done", round(time.time() - t0, 1), "s", flush=True)

    report = {"description": __doc__, "cases": cases}
    report.update(ck.summary())
    report["seconds"] = time.time() - t0
    report["optimize_flag"] = sys.flags.optimize
    target.write_text(json.dumps(report, indent=1))
    for c in cases:
        print("\n" + c["name"])
        for key, s in c["summary"].items():
            print(f"  {key:16s} mean={s['mean']:>14s} var={['%.6g' % v for v in s['variances']]}"
                  f" nondec={s['nondecreasing']} strict={s['strictly_increasing']}")
    print(target)
    ck.finish()


if __name__ == "__main__":
    main()
