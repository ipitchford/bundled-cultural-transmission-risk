"""Independent backward-ancestry validation of partial-transmission moments.

This uses ancestral-set rewards, not the forward count-moment implementation.
All rates are per unit time. beta_total is the homogeneous effective rate per
learner including the conventional self-choice probability; beta_total/N is
the effective rate of each directed non-self edge. b=0 denotes the limiting
independent-per-trait contact process.
"""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def reward_solution(states, events, reward):
    index = {s: j for j, s in enumerate(states)}
    rows, cols, vals = [], [], []
    for state in states:
        row = index[state]
        total = 0.0
        for target, rate in events(state):
            if rate <= 0 or target == state:
                continue
            total += rate
            if target in index:
                rows.append(row); cols.append(index[target]); vals.append(-rate)
        rows.append(row); cols.append(row); vals.append(total)
    minus_q = coo_matrix((vals, (rows, cols)), shape=(len(states), len(states))).tocsr()
    rhs = np.asarray([reward(s) for s in states], dtype=float)
    sol = spsolve(minus_q, rhs)
    residual = float(np.max(np.abs(minus_q @ sol - rhs)))
    return {s: float(sol[index[s]]) for s in states}, residual


def lumped(N, beta_total, reset, innovation, b, sample_size=None):
    """Backward single and pair reward equations after exchangeable lumping."""
    if sample_size is None:
        sample_size = N
    edge = beta_total / N
    def single_events(n):
        return [(n+1, edge*n*(N-n)), (n-1, reset*n)]
    first, residual1 = reward_solution(
        list(range(1, N+1)), single_events, lambda n: innovation*n)
    first[0] = 0.0
    states = [(i,j,k) for i in range(N+1) for j in range(N+1-i)
              for k in range(N+1-i-j) if i+k > 0 and j+k > 0]
    def pair_events(state):
        i,j,k = state
        empty = N-i-j-k
        return [
            ((i+1,j,k), edge*empty*(i+(1-b)*k)),
            ((i,j+1,k), edge*empty*(j+(1-b)*k)),
            ((i,j,k+1), edge*b*empty*k),
            ((i-1,j,k+1), edge*i*(j+k)),
            ((i,j-1,k+1), edge*j*(i+k)),
            ((i-1,j,k), reset*i),
            ((i,j-1,k), reset*j),
            ((i,j,k-1), reset*k),
        ]
    def rhs(state):
        i,j,k = state
        return innovation*((i+k)*first[j+k]+(j+k)*first[i+k])
    cross, residual2 = reward_solution(states, pair_events, rhs)
    mean = first[sample_size]
    factorial_second = cross[(0,0,sample_size)]
    return dict(mean=mean, variance=mean+factorial_second-mean**2,
                factorial_second=factorial_second, states=len(states),
                residual=max(residual1,residual2))


def labelled(beta, reset, innovation, b, sample=None):
    """Arbitrary directed-network solver; beta[receiver, donor]."""
    beta, reset, innovation = map(np.asarray, (beta,reset,innovation))
    N = len(reset)
    sample = (1 << N)-1 if sample is None else sample
    masks = list(range(1, 1 << N))
    def weight(mask):
        return sum(innovation[i] for i in range(N) if mask & (1 << i))
    def single_events(mask):
        out = []
        for i in range(N):
            if mask & (1 << i):
                out.append((mask & ~(1 << i), reset[i]))
                for j in range(N):
                    if not mask & (1 << j):
                        out.append((mask | (1 << j), beta[i,j]))
        return out
    first, r1 = reward_solution(masks, single_events, weight)
    first[0] = 0.0
    states = [(A,B) for A in masks for B in masks]
    def pair_events(state):
        A,B = state
        out = []
        for i in range(N):
            bit_i = 1 << i
            out.append(((A & ~bit_i, B & ~bit_i), reset[i]))
            for j in range(N):
                if i == j or beta[i,j] == 0:
                    continue
                bit_j = 1 << j
                newA = A | bit_j if A & bit_i else A
                newB = B | bit_j if B & bit_i else B
                out.append(((newA,B), beta[i,j]*(1-b)))
                out.append(((A,newB), beta[i,j]*(1-b)))
                out.append(((newA,newB), beta[i,j]*b))
        return out
    cross, r2 = reward_solution(states, pair_events,
        lambda state: weight(state[0])*first[state[1]]+weight(state[1])*first[state[0]])
    mean = first[sample]
    f2 = cross[(sample,sample)]
    return dict(mean=mean,variance=mean+f2-mean**2,
                factorial_second=f2,states=len(states),residual=max(r1,r2))


def n2_formula(edge, reset, innovation, b):
    rho, eta = edge/reset, innovation/reset
    mean = eta*(2+rho)
    scaled = (rho**3+(8-2*b)*rho**2+10*rho+4)/(2+(3-2*b)*rho)
    variance = mean+eta**2*scaled
    derivative = 2*eta**2*rho*(rho+1)*(rho+2)**2/(2+(3-2*b)*rho)**2
    inverse = (scaled*(2+3*rho)-(rho**3+8*rho**2+10*rho+4))/(2*rho*(scaled-rho))
    return dict(mean=mean,variance=variance,derivative=derivative,inverse=inverse)


class Checks:
    """Explicit checks (no assert): they also fail under python -O."""
    def __init__(self):
        self.count, self.failures = 0, []
    def check(self, condition, label, **detail):
        self.count += 1
        if not bool(condition):
            self.failures.append(dict(label=label, **detail))
            print("CHECK FAILED:", label, detail, file=sys.stderr, flush=True)
    def status(self):
        if not self.failures:
            return f"all {self.count} checks passed"
        return f"FAILED: {len(self.failures)} of {self.count} checks failed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="output directory (default: agent_math/, next to this script)")
    args = parser.parse_args()
    default = Path(__file__).resolve().parent
    outdir = default if args.out is None else Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / 'independent_pair_validation.json'
    if path.exists() and outdir != default:
        raise SystemExit(f"Refusing to overwrite existing file outside the default directory: {path}")
    ck = Checks()
    n2 = []
    for edge, reset, innovation in [(0.2,1.0,0.1),(1.0,1.0,0.5),(3.0,0.7,0.9)]:
        for b in [0.0,0.1,0.5,0.9,1.0]:
            got = lumped(2,2*edge,reset,innovation,b)
            expected = n2_formula(edge,reset,innovation,b)
            err = abs(got['variance']-expected['variance'])
            ck.check(err <= 1e-9*max(1.0,expected['variance']),
                     f"N=2 lumped vs formula edge={edge} b={b}", error=err)
            ck.check(abs(expected['inverse']-b) < 2e-13,
                     f"N=2 inverse edge={edge} b={b}", inverse=expected['inverse'])
            n2.append(dict(edge=edge,reset=reset,innovation=innovation,b=b,
                           **got,formula_variance=expected['variance'],error=err,
                           inverse=expected['inverse']))
    cross_checks = []
    for N in [2,3,4]:
        beta_total, reset_rate, innovation_rate = 1.3,1.0,0.1
        beta = np.full((N,N),beta_total/N); np.fill_diagonal(beta,0.0)
        for b in [0.0,0.35,1.0]:
            for sample_size in [1,N]:
                a = lumped(N,beta_total,reset_rate,innovation_rate,b,sample_size)
                c = labelled(beta,np.full(N,reset_rate),np.full(N,innovation_rate),
                             b,(1 << sample_size)-1)
                err = abs(a['variance']-c['variance'])
                ck.check(err < 1e-9, f"lumped vs labelled N={N} b={b} n={sample_size}", error=err)
                cross_checks.append(dict(N=N,b=b,sample_size=sample_size,
                                         lumped=a,labelled=c,error=err))
    beta = np.array([[0.,0.9,0.],[0.2,0.,0.4],[0.7,0.,0.]])
    reset = np.array([0.8,1.1,0.6]); innovation = np.array([0.1,0.3,0.2])
    network = [dict(b=b,**labelled(beta,reset,innovation,b))
               for b in [0.0,0.1,0.3,0.6,0.9,1.0]]
    ck.check(max(x['mean'] for x in network)-min(x['mean'] for x in network) < 1e-12,
             "directed network: mean constant in b")
    ck.check(all(x['variance'] < y['variance'] for x,y in zip(network,network[1:])),
             "directed network: variance strictly increasing in b",
             variances=[x['variance'] for x in network])
    out=dict(n2_formula_checks=n2,labelled_lumped_checks=cross_checks,
             directed_heterogeneous_example=network,
             max_n2_absolute_error=max(x['error'] for x in n2),
             max_labelled_lumped_absolute_error=max(x['error'] for x in cross_checks),
             status=ck.status(), checks_run=ck.count, failures=ck.failures,
             optimize_flag=sys.flags.optimize)
    path.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k] for k in ['max_n2_absolute_error',
                      'max_labelled_lumped_absolute_error','directed_heterogeneous_example']},indent=2))
    print(path)
    print(ck.status(), flush=True)
    if ck.failures:
        sys.exit(1)


if __name__ == '__main__':
    main()
