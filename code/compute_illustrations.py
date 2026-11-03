"""Reproduce the larger-population illustrations.

Usage: python3 -I -B code/compute_illustrations.py [--large] [--out DIR]

Default: six N=50 values, using a direct sparse solve. Optional --large adds
the N=100, b=1/2 source-scale comparison with a Jacobi-preconditioned
BiCGSTAB solve (repeated at N=50 as a check against the direct solve). Both solve
the same exact finite Markov reward system. The full-copy endpoint is also
computed by a positive one-dimensional recursion.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab
from certify import Checker, DEFAULT_RESULTS, output_dir, check_writable
from transmission import pair_generator, single_moments, repertoire_moments


def iterative(N):
    start = time.time()
    Q, states, index = pair_generator(N, 1.3, 1.0, 0.5)
    A = (-Q).tocsc()
    m, m2 = single_moments(N, 1.3, 1.0)
    rhs = np.array([(i + k) * m[j + k] + (j + k) * m[i + k]
                    for i, j, k in states])
    # Jacobi-preconditioned BiCGSTAB (the earlier ILU/GMRES variant was slow at N=100).
    dinv = 1.0 / A.diagonal()
    iterations = [0]
    def callback(_):
        iterations[0] += 1
    h, info = bicgstab(A.tocsr(), rhs, M=LinearOperator(A.shape, lambda x: dinv * x),
                       rtol=1e-13, atol=0, maxiter=200000, callback=callback)
    if info:
        raise RuntimeError(f"BiCGSTAB did not converge: {info}")
    absolute_residual = float(np.max(np.abs(A @ h - rhs)))
    relative_residual = absolute_residual / float(np.max(rhs))
    increments = np.zeros(N + 2)
    for k in range(N, 0, -1):
        increments[k] = 1 / k + 1.3 * (N - k) / N * increments[k + 1]
    maximum_absorption_mean = float(increments.sum())
    variance = 0.1 * m[N] + 0.01 * (h[index[(0, 0, N)]] - m[N] ** 2)
    return dict(N=N, b=0.5, mean=float(0.1 * m[N]), variance=float(variance),
                residual=relative_residual, residual_infinity=absolute_residual,
                max_single_absorption_mean=maximum_absorption_mean,
                variance_residual_error_bound=0.01 * maximum_absorption_mean * absolute_residual,
                bound_note="M-matrix residual bound evaluated in ordinary floating point, not interval arithmetic.",
                full_copy_variance=float(0.1 * m[N] + 0.01 * (m2[N] - m[N] ** 2)),
                info=int(info), iterations=iterations[0], seconds=time.time() - start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--large", action="store_true", help="also run the 100-person iterative calculation")
    parser.add_argument("--out", help="output directory (default: results/)")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_RESULTS)
    targets = [check_writable(out / "large_N_moments.json", DEFAULT_RESULTS)]
    if args.large:
        targets += [check_writable(out / f"iterative_N{N}.json", DEFAULT_RESULTS) for N in [50, 100]]
    ck = Checker()
    values = []
    for b in [0, 0.1, 0.5, 0.9, 0.99, 1]:
        result = repertoire_moments(50, 1.3, u=1, mu=0.1, b=b)
        values.append(result)
        ck.check(result["relative_residual"] < 1e-10, f"N=50 direct solve residual b={b}",
                 residual=result["relative_residual"])
        print(f"N=50 b={b:g} mean={result['mean']:.8f} variance={result['variance']:.8f}", flush=True)
    ck.check(np.ptp([v["mean"] for v in values]) < 1e-9, "N=50 mean constant in b")
    ck.check(all(x["variance"] < y["variance"] for x, y in zip(values, values[1:])),
             "N=50 variance strictly increasing in b")
    ck.close("N=50 b=1 pair system vs single-genealogy recursion", values[-1]["weighted_ancestry_second"],
             values[-1]["full_copy_genealogy_second"], 1e-9, relative=True)
    targets[0].write_text(json.dumps(values, indent=2))
    if args.large:
        direct50 = values[2]["variance"]
        for N, target in zip([50, 100], targets[1:]):
            result = iterative(N)
            ck.check(result["residual"] < 1e-9, f"iterative N={N} relative residual", residual=result["residual"])
            if N == 50:
                ck.close("iterative vs direct N=50 b=1/2 variance", result["variance"], direct50, 1e-8)
            target.write_text(json.dumps(result, indent=2))
            print(result, flush=True)
    ck.finish()
