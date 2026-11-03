"""Partial forgetting extension: exact forward-closure checks.

Usage: python3 -I -B code/partial_forgetting_check.py [--out DIR]

Generalised model: loss events at site i occur at rate u_i/d and delete EACH
trait held there independently with probability d (d = 1 is the paper's
whole-site reset); encounters j -> i occur at rate a_ji/b and copy each donor
trait independently with probability b. Every trait therefore sees the same
per-trait loss rate u_i and learning rate a_ji for all (b, d). The limits
b -> 0 and d -> 0 are taken analytically in the closure.

Exact rational checks, homogeneous N=2 (u=1, beta=2, mu=1/2) and N=3
(u=1, beta=13/10, mu=1/5), over (b, d) in {0 (limit), 1/4, 1/2, 1}^2:
  * means of C_V, of every P_k and every G_k (held by >= k sites) invariant;
  * variances nondecreasing in b at fixed d and in d at fixed b;
  * at b = d = 0 (both limits) every count is Poisson: Var = mean;
  * at d = 1 the closure reproduces the backward pair system.
Also derives the N=2 closed form of Var C_V as a function of (b, d, rho, eta)
symbolically with sympy and checks that it reduces to the paper's formula at
d = 1. Writes results/partial_forgetting.json.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import time
from fractions import Fraction as Fr

import sympy as sp

from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from forward_closure import ForwardMoments, backward_union, homogeneous


def symbolic_n2():
    """Independent symbolic forward closure for N=2 (u=1, edge rate rho, innovation eta)."""
    b, d, rho, eta = sp.symbols("b d rho eta", positive=True)
    # carrier sets: x = {1} only, y = {2} only, z = {1,2}
    names = ["x", "y", "z"]
    m = sp.symbols("mx my mz")
    S = {}
    for p in range(3):
        for q in range(p, 3):
            S[(p, q)] = S[(q, p)] = sp.Symbol(f"s_{names[p]}{names[q]}")
    Sig = sp.Matrix(3, 3, lambda p, q: S[(p, q)])
    mv = sp.Matrix(m)
    X, Y, Z = 0, 1, 2

    def vec(**kw):
        v = sp.zeros(3, 1)
        for k, val in kw.items():
            v[names.index(k)] = val
        return v

    # families: (per-trait rate, success prob, {source index: move vector})
    fam = [
        (1, d, {X: vec(x=-1), Z: vec(z=-1, y=1)}),       # loss at site 1
        (1, d, {Y: vec(y=-1), Z: vec(z=-1, x=1)}),       # loss at site 2
        (rho, b, {X: vec(x=-1, z=1)}),                   # encounter 1 -> 2
        (rho, b, {Y: vec(y=-1, z=1)}),                   # encounter 2 -> 1
    ]
    innov = [vec(x=1), vec(y=1)]
    dm = sp.zeros(3, 1)
    dS = sp.zeros(3, 3)
    for rate, p, moves in fam:
        K = sp.zeros(3, 3)
        for s, v in moves.items():
            K[:, s] = v
        dm += rate * K * mv
    for e in innov:
        dm += eta * e
    msol = sp.solve(list(dm), list(m), dict=True)[0]
    mval = mv.subs(msol)
    for rate, p, moves in fam:
        K = sp.zeros(3, 3)
        noise = sp.zeros(3, 3)
        for s, v in moves.items():
            K[:, s] = v
            noise += mval[s] * v * v.T
        dS += rate * (K * Sig + Sig * K.T + p * K * Sig * K.T + (1 - p) * noise)
    for e in innov:
        dS += eta * (mval * e.T + e * mval.T + e * e.T)
    unknowns = sorted({S[k] for k in S}, key=str)
    eqs = [dS[p, q] for p in range(3) for q in range(p, 3)]
    sol = sp.solve(eqs, unknowns, dict=True)[0]
    Sv = Sig.subs(sol)
    one = sp.Matrix([1, 1, 1])
    mean = sp.factor((one.T * mval)[0])
    var = sp.factor((one.T * Sv * one)[0] - mean ** 2)
    excess = sp.factor(sp.simplify((var - mean) / eta ** 2))
    paper = eta * (2 + rho) + eta ** 2 * (rho ** 3 + (8 - 2 * b) * rho ** 2 + 10 * rho + 4) / (2 + (3 - 2 * b) * rho)
    at_d1 = sp.simplify(var.subs(d, 1) - paper)
    num, den = sp.fraction(sp.together(excess))
    return dict(symbols=(b, d, rho, eta), mean=mean, var=var, excess=excess,
                numerator=sp.expand(num), denominator=sp.expand(den),
                reduces_to_paper_at_d1=(at_d1 == 0),
                dV_db=sp.factor(sp.diff(var, b)), dV_dd=sp.factor(sp.diff(var, d)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="output directory (default: results/)")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    target = check_writable(out / "partial_forgetting.json", DEFAULT_RESULTS)
    ck = Checker()
    t0 = time.time()
    grid = [Fr(0), Fr(1, 4), Fr(1, 2), Fr(1)]
    report = {"description": __doc__, "cases": []}

    for N, beta, mu in [(2, Fr(2), Fr(1, 2)), (3, Fr(13, 10), Fr(1, 5))]:
        u, m, a = homogeneous(N, beta, 1, mu)
        V = list(range(N))
        names = ["C_V"] + [f"P_{k}" for k in range(1, N + 1)] + [f"G_{k}" for k in range(2, N + 1)]
        table = {}  # (b,d) -> name -> (mean, var)
        for b in grid:
            for d in grid:
                fm = ForwardMoments(N, u, m, a, b, d)
                vals = {"C_V": fm.union(V)}
                for k in range(1, N + 1):
                    vals[f"P_{k}"] = fm.exactly(k)
                for k in range(2, N + 1):
                    vals[f"G_{k}"] = fm.at_least(k)
                table[(b, d)] = vals
                if d == 1:
                    bm, bv = backward_union(N, u, m, a, b, V)
                    ck.check(vals["C_V"] == (bm, bv), f"N={N} d=1 b={b}: forward == backward pair system")
        case = {"N": N, "beta": str(beta), "mu": str(mu), "u": "1", "grid": [str(g) for g in grid],
                "values": {}, "checks": {}}
        for name in names:
            means = {key: table[key][name][0] for key in table}
            ck.check(len(set(means.values())) == 1, f"N={N}: mean of {name} invariant over (b,d) grid")
            nondec_b = all(table[(grid[p], d)][name][1] <= table[(grid[p + 1], d)][name][1]
                           for d in grid for p in range(3))
            nondec_d = all(table[(b, grid[p])][name][1] <= table[(b, grid[p + 1])][name][1]
                           for b in grid for p in range(3))
            strict_b = all(table[(grid[p], d)][name][1] < table[(grid[p + 1], d)][name][1]
                           for d in grid for p in range(3))
            strict_d = all(table[(b, grid[p])][name][1] < table[(b, grid[p + 1])][name][1]
                           for b in grid for p in range(3))
            ck.check(nondec_b, f"N={N}: Var {name} nondecreasing in b at every d")
            ck.check(nondec_d, f"N={N}: Var {name} nondecreasing in d at every b")
            mean0, var0 = table[(Fr(0), Fr(0))][name]
            ck.check(var0 == mean0, f"N={N}: {name} Poisson at b=d=0 (Var = mean)",
                     mean=str(mean0), var=str(var0))
            case["checks"][name] = {"mean": str(next(iter(means.values()))),
                                   "nondecreasing_in_b": nondec_b, "nondecreasing_in_d": nondec_d,
                                   "strict_in_b": strict_b, "strict_in_d": strict_d,
                                   "poisson_at_b0_d0": var0 == mean0}
            case["values"][name] = {f"b={b},d={d}": {"variance": str(table[(b, d)][name][1]),
                                                     "float": float(table[(b, d)][name][1])}
                                    for (b, d) in table}
        report["cases"].append(case)
        print(f"N={N} done ({time.time() - t0:.1f}s)", flush=True)
        print("  Var C_V rows b, columns d =", [str(g) for g in grid])
        for b in grid:
            print(f"   b={str(b):>4s}", [f"{float(table[(b, d)]['C_V'][1]):.6f}" for d in grid])

    sym = symbolic_n2()
    b, d, rho, eta = sym["symbols"]
    ck.check(sym["reduces_to_paper_at_d1"], "N=2 symbolic Var C_V reduces to the paper's formula at d=1")
    ck.check(sp.simplify(sym["mean"] - eta * (2 + rho)) == 0, "N=2 symbolic mean is eta(2+rho)")
    ck.check(sp.simplify(sym["var"].subs({b: 0, d: 0}) - sym["mean"]) == 0,
             "N=2 symbolic Var C_V Poisson at b=d=0")
    # A structured form of the same closed form (hand-collected; checked here):
    A_ = (rho**2 + 2*rho + 2)*d**2 + (rho**4 + rho**3 - 7*rho**2 - 16*rho - 10)*d \
        + 4*(rho + 1)*(2*rho**2 + 4*rho + 3)
    B_ = d**2 + (rho**2 - 3*rho - 5)*d + 2*(rho + 1)*(rho + 3)
    structured = (d*A_ - 2*b*rho*(rho + 1)*(d*(rho + 1) - 1)) / ((2 - d)*B_ - 2*b*rho*(rho + 2 - d))
    ck.check(sp.simplify(structured - sym["excess"]) == 0, "N=2 structured closed form equals sympy result")
    # numeric agreement of symbolic formula with the rational grid
    n2 = report["cases"][0]["values"]["C_V"]
    for key, rec in n2.items():
        bv, dv = (Fr(x.split("=")[1]) for x in key.split(","))
        val = sym["var"].subs({b: sp.Rational(bv.numerator, bv.denominator),
                               d: sp.Rational(dv.numerator, dv.denominator), rho: 1, eta: sp.Rational(1, 2)})
        ck.check(sp.nsimplify(val) == sp.Rational(rec["variance"]), f"N=2 symbolic formula matches grid at {key}")
    report["n2_closed_form"] = {
        "parameters": "u=1 (time unit), rho = beta/(2u) directed-edge rate, eta = mu/u",
        "mean": str(sym["mean"]),
        "variance": str(sym["var"]),
        "excess_over_mean_divided_by_eta2": str(sym["excess"]),
        "excess_numerator": str(sym["numerator"]),
        "excess_denominator": str(sym["denominator"]),
        "dVar_db": str(sym["dV_db"]),
        "dVar_dd": str(sym["dV_dd"]),
        "latex_excess": sp.latex(sym["excess"]),
        "structured_excess": "[d*A - 2*b*rho*(rho+1)*(d*(rho+1)-1)] / [(2-d)*B - 2*b*rho*(rho+2-d)], "
                             "A = (rho^2+2rho+2)d^2 + (rho^4+rho^3-7rho^2-16rho-10)d + 4(rho+1)(2rho^2+4rho+3), "
                             "B = d^2 + (rho^2-3rho-5)d + 2(rho+1)(rho+3)",
        "reduces_to_paper_formula_at_d1": bool(sym["reduces_to_paper_at_d1"]),
    }
    print("N=2 (Var C_V - E C_V)/eta^2 =", sym["excess"])
    print("dVar/db =", sym["dV_db"])
    print("dVar/dd =", sym["dV_dd"])
    report.update(ck.summary())
    report["seconds"] = time.time() - t0
    report["optimize_flag"] = sys.flags.optimize
    target.write_text(json.dumps(report, indent=1))
    print(target)
    ck.finish()


if __name__ == "__main__":
    main()
