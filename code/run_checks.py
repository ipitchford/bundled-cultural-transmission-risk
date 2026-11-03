"""Reproduce exact, labelled-state and independent forward checks.

Usage:
    python3 -I -B code/run_checks.py [--out DIR] [--reuse-forward FILE]

All checks are explicit (no `assert`), so they also fail under `python -O`.
The status string is set from the recorded results and the script exits with
status 1 if any check fails. --reuse-forward re-verifies the saved forward
simulations from an earlier verification.json instead of rerunning the
400,000-cycle simulation (the deterministic exact checks are always rerun).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import time
from fractions import Fraction

import numpy as np

from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from transmission import (exact_pair_second_N2, repertoire_moments, n2_population_variance,
                          n2_identify_b, labelled_pair_moments, stationary_forward_cycles)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="output directory (default: results/)")
    parser.add_argument("--reuse-forward", metavar="FILE",
                        help="re-check saved forward simulations instead of rerunning them")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    target = check_writable(out / "verification.json", DEFAULT_RESULTS)
    ck = Checker()
    t0 = time.time()
    report = {"convention": "beta is effective per-learner rate, directed edge beta/N",
              "exact": [], "labelled": [], "degenerate": [], "b1": [], "forward": []}

    for b in [Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1)]:
        h = exact_pair_second_N2(b=b)
        got = repertoire_moments(2, 2, mu=.5, b=float(b))
        mean, var = n2_population_variance(2, 1, .5, float(b))
        ck.close(f"N=2 rational vs sparse pair second moment, b={b}", got["weighted_ancestry_second"], float(h), 1e-11)
        ck.close(f"N=2 sparse variance vs closed form, b={b}", got["variance"], var, 1e-11)
        recovered = n2_identify_b(1, .5, var)
        ck.close(f"N=2 inverse recovers b={b}", recovered, float(b), 1e-12)
        report["exact"].append(dict(b=str(b), E_W2=str(h), mean=mean, variance=var, identified_b=recovered))

    for N in [2, 3, 4]:
        for n in [1, N]:
            for b in [0, .2, .65, 1]:
                got = repertoire_moments(N, 1.3, mu=.2, b=b, n=n)
                labelled = labelled_pair_moments(N, 1.3, 1, b, n)
                error = abs(labelled - got["weighted_ancestry_second"])
                ck.check(error < 1e-10, f"labelled vs lumped N={N} n={n} b={b}", error=error)
                report["labelled"].append(dict(N=N, n=n, b=b, error=error))

    for N, n in [(1, 1), (5, 1), (5, 3), (10, 10)]:
        for b in [0, .3, 1]:
            got = repertoire_moments(N, 0, mu=.3, b=b, n=n)
            ck.close(f"no-learning mean N={N} n={n} b={b}", got["mean"], .3 * n, 1e-12)
            ck.close(f"no-learning variance N={N} n={n} b={b}", got["variance"], (.3 + .3 ** 2) * n, 1e-12)
            report["degenerate"].append(dict(N=N, n=n, b=b, mean=got["mean"], variance=got["variance"]))

    for N in [2, 5, 10, 20, 30]:
        for beta in [.6, 1, 1.3]:
            got = repertoire_moments(N, beta, mu=.1, b=1)
            error = abs(got["weighted_ancestry_second"] - got["full_copy_genealogy_second"]) / got["full_copy_genealogy_second"]
            ck.check(error < 1e-10, f"b=1 pair system vs single-genealogy recursion N={N} beta={beta}", relative_error=error)
            report["b1"].append(dict(N=N, beta=beta, relative_error=error))

    if args.reuse_forward:
        saved = json.loads(Path(args.reuse_forward).read_text())["forward"]
        report["forward_source"] = f"re-checked saved simulations from {Path(args.reuse_forward).name}; not rerun"
        sims = saved
    else:
        report["forward_source"] = "simulated in this run"
        sims = None
    for j, b in enumerate([.1, .25, .5, 1]):
        if sims is None:
            sim = stationary_forward_cycles(2, 2, 1, .5, b, 100000, seed=20261102 + 100 * j)
        else:
            sim = dict(sims[j])
            ck.check(sim["b"] == b and sim["cycles"] == 100000 and sim["seed"] == 20261102 + 100 * j,
                     f"saved forward run {j} has the prespecified b, cycles and seed", b=sim.get("b"))
        mean, var = n2_population_variance(2, 1, .5, b)
        sim["exact_mean"] = mean
        sim["exact_variance"] = var
        sim["mean_z"] = (sim["mean"] - mean) / sim["mean_se"]
        sim["variance_z"] = (sim["variance"] - var) / sim["variance_se"]
        ck.check(abs(sim["mean_z"]) < 6, f"forward simulation mean within 6 SE, b={b}", z=sim["mean_z"])
        ck.check(abs(sim["variance_z"]) < 6, f"forward simulation variance within 6 SE, b={b}", z=sim["variance_z"])
        report["forward"].append(sim)
        if sims is None:
            report["status"] = "in progress"
            target.write_text(json.dumps(report, indent=2))
        print("forward", b, "mean", sim["mean"], "z", round(sim["mean_z"], 3), "var", sim["variance"],
              "z", round(sim["variance_z"], 3), "P0", sim["p_zero"], flush=True)

    report["max_abs_mean_z"] = max(abs(s["mean_z"]) for s in report["forward"])
    report["max_abs_variance_z"] = max(abs(s["variance_z"]) for s in report["forward"])
    report.update(ck.summary())
    report["software"] = {"python": sys.version, "numpy": np.__version__,
                          "optimize_flag": sys.flags.optimize}
    report["seconds"] = time.time() - t0
    target.write_text(json.dumps(report, indent=2))
    print(target, flush=True)
    ck.finish()


if __name__ == "__main__":
    main()
