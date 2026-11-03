"""Exact FORWARD moment closure on carrier-set counts, plus an exact labelled
BACKWARD pair-replica solver, both in rational arithmetic.

Forward state: n_S = number of traits whose carrier set is exactly S, for each
nonempty subset S of the N sites. Events are affine in n, so the first and
second moments close exactly:

  innovation at i (rate mu_i):     n_{{i}} += 1
  loss event at i (rate u_i/d):    for each S containing i, Binomial(n_S, d)
                                   traits move S -> S minus i (lost if empty).
                                   d = 1 is the whole-site reset of the paper.
  encounter j -> i (rate a_ji/b):  for each S with j in S and i not in S,
                                   Binomial(n_S, b) traits move S -> S plus i.

The limits b -> 0 and d -> 0 are taken analytically: with v_S the move vector
of one trait, the drift of E[n n^T] from an event family of per-trait rate r
and per-trait success probability p is

  r [K S + S K^T + p K S K^T + (1 - p) sum_S m_S v_S v_S^T],

which is polynomial in p. This representation shares no code with the
backward pair systems in transmission.py and agent_math/pair_ancestry.py.

Conventions: a[j][i] is the effective rate of the directed edge j -> i
(donor j, learner i); a[i][i] is ignored (silent self-choice).
"""
from __future__ import annotations

from fractions import Fraction as Fr

from sympy import QQ
from sympy.polys.matrices import DomainMatrix


# ----------------------------------------------------------------------------
# exact linear algebra

def _q(x):
    x = Fr(x)
    return QQ(x.numerator, x.denominator)


def solve_rational(rows: list[dict[int, Fr]], rhs: list[Fr]) -> list[Fr]:
    """Solve A x = rhs exactly; A given as sparse rows {col: value}."""
    n = len(rhs)
    dense = [[QQ(0)] * n for _ in range(n)]
    for r, row in enumerate(rows):
        for c, v in row.items():
            if v:
                dense[r][c] = _q(v)
    A = DomainMatrix(dense, (n, n), QQ)
    B = DomainMatrix([[_q(v)] for v in rhs], (n, 1), QQ)
    X = A.lu_solve(B)
    out = []
    for i in range(n):
        e = X[i, 0].element
        out.append(Fr(int(e.numerator), int(e.denominator)))
    return out


# ----------------------------------------------------------------------------
# forward closure

def homogeneous(N: int, beta, u=1, mu=Fr(1, 10)):
    """Homogeneous complete graph: a_ji = beta/N for i != j."""
    beta, u, mu = Fr(beta), Fr(u), Fr(mu)
    a = [[beta / N if i != j else Fr(0) for i in range(N)] for j in range(N)]
    return [u] * N, [mu] * N, a


class ForwardMoments:
    """Stationary first and second moments of the carrier-set counts n_S."""

    def __init__(self, N, u, mu, a, b, d=1):
        self.N = N
        self.b, self.d = Fr(b), Fr(d)
        if not (0 <= self.b <= 1 and 0 <= self.d <= 1):
            raise ValueError("b and d must lie in [0, 1]")
        u = [Fr(x) for x in u]
        mu = [Fr(x) for x in mu]
        a = [[Fr(x) for x in row] for row in a]
        self.subsets = list(range(1, 1 << N))
        K = len(self.subsets)
        idx = {S: p for p, S in enumerate(self.subsets)}
        self.idx = idx

        # Event families: (per-trait rate, success probability, {source col: {row: coef}})
        families = []
        for i in range(N):
            if u[i] == 0:
                continue
            moves = {}
            for S in self.subsets:
                if (S >> i) & 1:
                    T = S & ~(1 << i)
                    v = {idx[S]: Fr(-1)}
                    if T:
                        v[idx[T]] = Fr(1)
                    moves[idx[S]] = v
            families.append((u[i], self.d, moves))
        for j in range(N):
            for i in range(N):
                if i == j or a[j][i] == 0:
                    continue
                moves = {}
                for S in self.subsets:
                    if (S >> j) & 1 and not (S >> i) & 1:
                        moves[idx[S]] = {idx[S]: Fr(-1), idx[S | (1 << i)]: Fr(1)}
                families.append((a[j][i], self.b, moves))

        # First moments: 0 = sum_f r_f K_f m + sum_i mu_i e_i
        rows = [dict() for _ in range(K)]
        rhs = [Fr(0)] * K
        for rate, _, moves in families:
            for s, v in moves.items():
                for x, cf in v.items():
                    rows[x][s] = rows[x].get(s, Fr(0)) + rate * cf
        for i in range(N):
            rhs[idx[1 << i]] -= mu[i]
        m = solve_rational(rows, rhs)
        self.m = m

        # Second moments E[n_x n_y], x <= y.
        pairs = [(x, y) for x in range(K) for y in range(x, K)]
        pid = {p: q for q, p in enumerate(pairs)}

        def key(x, y):
            return pid[(x, y) if x <= y else (y, x)]

        P = len(pairs)
        rows2 = [dict() for _ in range(P)]
        rhs2 = [Fr(0)] * P

        def add(r, col, val):
            if val:
                rows2[r][col] = rows2[r].get(col, Fr(0)) + val

        for rate, prob, moves in families:
            krow = {}  # x -> [(s, coef)] : (K)_{x s}
            for s, v in moves.items():
                for x, cf in v.items():
                    krow.setdefault(x, []).append((s, cf))
            for (x, y) in pairs:
                r = pid[(x, y)]
                for s, cf in krow.get(x, []):
                    add(r, key(s, y), rate * cf)
                for s, cf in krow.get(y, []):
                    add(r, key(x, s), rate * cf)
                if prob:
                    for s, cf in krow.get(x, []):
                        for t, cg in krow.get(y, []):
                            if s == t:
                                # E[X_s^2] part: p^2 E n_s^2 + p(1-p) E n_s;
                                # after dividing by p: p E n_s^2 + (1-p) m_s
                                add(r, key(s, s), rate * prob * cf * cg)
                            else:
                                add(r, key(s, t), rate * prob * cf * cg)
                if prob != 1:
                    for s, v in moves.items():
                        if x in v and y in v:
                            rhs2[r] -= rate * (1 - prob) * m[s] * v[x] * v[y]
        for i in range(N):
            ei = idx[1 << i]
            for (x, y) in pairs:
                r = pid[(x, y)]
                val = Fr(0)
                if y == ei:
                    val += m[x]
                if x == ei:
                    val += m[y]
                if x == ei and y == ei:
                    val += 1
                rhs2[r] -= mu[i] * val
        sig = solve_rational(rows2, rhs2)
        self._sig, self._key = sig, key

    def second(self, x, y):
        return self._sig[self._key(x, y)]

    def functional(self, coeff):
        """Mean and variance of sum_S coeff(S) n_S (coeff integer-valued)."""
        c = [Fr(coeff(S)) for S in self.subsets]
        nz = [p for p, v in enumerate(c) if v]
        mean = sum(c[p] * self.m[p] for p in nz)
        second = sum(c[p] * c[q] * self.second(p, q) for p in nz for q in nz)
        return mean, second - mean * mean

    # Common functionals -------------------------------------------------
    def union(self, sample):
        mask = sum(1 << x for x in sample)
        return self.functional(lambda S: 1 if S & mask else 0)

    def exactly(self, k):
        return self.functional(lambda S: 1 if bin(S).count("1") == k else 0)

    def at_least(self, k):
        return self.functional(lambda S: 1 if bin(S).count("1") >= k else 0)


# ----------------------------------------------------------------------------
# labelled backward pair-replica system (exact)

def backward_union(N, u, mu, a, b, sample):
    """Mean and variance of C_A from the labelled backward pair system."""
    b = Fr(b)
    u = [Fr(x) for x in u]
    mu = [Fr(x) for x in mu]
    a = [[Fr(x) for x in row] for row in a]
    masks = list(range(1, 1 << N))
    i1 = {B: p for p, B in enumerate(masks)}

    def g(B):
        return sum((mu[i] for i in range(N) if (B >> i) & 1), Fr(0))

    rows = []
    for B in masks:
        tr = {}
        for i in range(N):
            if (B >> i) & 1:
                D = B & ~(1 << i)
                tr[D] = tr.get(D, Fr(0)) + u[i]
        for j in range(N):
            for i in range(N):
                if i != j and a[j][i] and (B >> i) & 1 and not (B >> j) & 1:
                    D = B | (1 << j)
                    tr[D] = tr.get(D, Fr(0)) + a[j][i]
        row = {i1[B]: sum(tr.values(), Fr(0))}
        for D, rate in tr.items():
            if D:
                row[i1[D]] = row.get(i1[D], Fr(0)) - rate
        rows.append(row)
    m = solve_rational(rows, [g(B) for B in masks])

    def mm(B):
        return m[i1[B]] if B else Fr(0)

    states = [(x, y) for x in masks for y in masks]
    i2 = {s: p for p, s in enumerate(states)}
    rows2, rhs2 = [], []
    for (B1, B2) in states:
        tr = {}

        def add(D, rate):
            if D != (B1, B2) and rate:
                tr[D] = tr.get(D, Fr(0)) + rate

        for i in range(N):
            add((B1 & ~(1 << i), B2 & ~(1 << i)), u[i])
        for j in range(N):
            for i in range(N):
                if i == j or not a[j][i]:
                    continue
                x1 = B1 | (1 << j) if (B1 >> i) & 1 else B1
                x2 = B2 | (1 << j) if (B2 >> i) & 1 else B2
                add((x1, x2), a[j][i] * b)
                add((x1, B2), a[j][i] * (1 - b))
                add((B1, x2), a[j][i] * (1 - b))
        row = {i2[(B1, B2)]: sum(tr.values(), Fr(0))}
        for D, rate in tr.items():
            if D[0] and D[1]:
                row[i2[D]] = row.get(i2[D], Fr(0)) - rate
        rows2.append(row)
        rhs2.append(g(B1) * mm(B2) + g(B2) * mm(B1))
    F = solve_rational(rows2, rhs2)
    A = sum(1 << x for x in sample)
    mean = m[i1[A]]
    return mean, mean + F[i2[(A, A)]] - mean * mean
