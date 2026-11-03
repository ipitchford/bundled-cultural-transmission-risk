"""Exact moments for neutral cultural transmission with common encounters.

Convention: N sites, death rate u at each site, innovation rate mu at each
site, effective per-learner transmission rate beta. A distinct ordered pair
has encounter rate beta/(N*b); every parent trait is copied independently
with probability b. Self encounters have no effect.

The b=0 argument in moment routines is the limit of infinitely frequent,
infinitesimal copying, not a literal finite-rate encounter model.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def single_moments(N: int, beta: float, u: float = 1.0):
    """E L_n and E L_n^2 for all n, using positive birth-death recursions."""
    if N < 1 or beta < 0 or u <= 0:
        raise ValueError("Require N>=1, beta>=0, u>0")
    delta = np.zeros(N + 2)
    for j in range(N, 0, -1):
        delta[j] = 1 / u + beta * (N - j) / (N * u) * delta[j + 1]
    mean = np.r_[0.0, np.cumsum(delta[1:N + 1])]
    delta2 = np.zeros(N + 2)
    for j in range(N, 0, -1):
        delta2[j] = 2 * mean[j] / u + beta * (N - j) / (N * u) * delta2[j + 1]
    second = np.r_[0.0, np.cumsum(delta2[1:N + 1])]
    return mean, second


def pair_states(N):
    """Counts (A-only,B-only,both); stop when either ancestral set is empty."""
    return [(i, j, k) for k in range(N + 1)
            for i in range(N - k + 1) for j in range(N - k - i + 1)
            if i + k > 0 and j + k > 0]


def pair_transitions(state, N, beta, u, b):
    """Off-diagonal rates, derived directly from the eight event types."""
    i, j, k = state
    d = N - i - j - k
    events = [
        ((i - 1, j, k), u * i),
        ((i, j - 1, k), u * j),
        ((i, j, k - 1), u * k),
        ((i + 1, j, k), beta / N * d * (i + (1 - b) * k)),
        ((i, j + 1, k), beta / N * d * (j + (1 - b) * k)),
        ((i, j, k + 1), beta / N * b * d * k),
        ((i - 1, j, k + 1), beta / N * i * (j + k)),
        ((i, j - 1, k + 1), beta / N * j * (i + k)),
    ]
    return [(s, rate) for s, rate in events if rate > 0]


def pair_generator(N, beta, u, b):
    if not 0 <= b <= 1:
        raise ValueError("b must lie in [0,1]")
    states = pair_states(N)
    index = {s: p for p, s in enumerate(states)}
    rows, cols, vals = [], [], []
    for row, state in enumerate(states):
        outgoing = 0.0
        for dest, rate in pair_transitions(state, N, beta, u, b):
            outgoing += rate
            if dest in index:
                rows.append(row); cols.append(index[dest]); vals.append(rate)
        rows.append(row); cols.append(row); vals.append(-outgoing)
    Q = coo_matrix((vals, (rows, cols)), shape=(len(states), len(states))).tocsr()
    return Q, states, index


def repertoire_moments(N, beta, u=1.0, mu=0.1, b=1.0, n=None):
    """Exact finite linear-system result, evaluated in floating point."""
    n = N if n is None else n
    if not 1 <= n <= N or mu < 0:
        raise ValueError("Invalid sample size or innovation rate")
    m, m2 = single_moments(N, beta, u)
    Q, states, index = pair_generator(N, beta, u, b)
    rhs = np.array([(i + k) * m[j + k] + (j + k) * m[i + k]
                    for i, j, k in states])
    cross = spsolve(-Q, rhs)
    residual = np.max(np.abs((-Q) @ cross - rhs)) / max(1.0, np.max(rhs))
    W2 = float(cross[index[(0, 0, n)]])
    Wvar = W2 - m[n] ** 2
    mean = mu * m[n]
    return dict(N=N, n=n, beta=beta, u=u, mu=mu, b=b,
                mean=float(mean), variance=float(mean + mu ** 2 * Wvar),
                environmental_variance=float(mu ** 2 * Wvar),
                genealogy_mean=float(m[n]), weighted_ancestry_second=W2,
                full_copy_genealogy_second=float(m2[n]),
                states=len(states), relative_residual=float(residual))


def n2_population_variance(beta, u, mu, b):
    """Closed expression; rho uses directed-edge, not per-learner, rate."""
    rho, eta = beta / (2 * u), mu / u
    mean = eta * (2 + rho)
    V = (rho ** 3 + (8 - 2 * b) * rho ** 2 + 10 * rho + 4) / (2 + (3 - 2 * b) * rho)
    return mean, mean + eta ** 2 * V


def n2_identify_b(mean_singletons, mean_doubletons, variance):
    """Population-law identification; no finite-sample uncertainty claimed."""
    eta = mean_singletons / 2
    rho = mean_doubletons / eta
    mean = mean_singletons + mean_doubletons
    V = (variance - mean) / eta ** 2
    A = rho ** 3 + 8 * rho ** 2 + 10 * rho + 4
    return (V * (2 + 3 * rho) - A) / (2 * rho * (V - rho))


def exact_pair_second_N2(beta=Fraction(2), u=Fraction(1), b=Fraction(1, 2), n=2):
    """Five-state rational Gaussian elimination, independent of sparse solver."""
    N = 2
    states = pair_states(N)
    idx = {s: j for j, s in enumerate(states)}
    m = [Fraction(0), (1 + beta / (2 * u)) / u, (2 + beta / (2 * u)) / u]
    A = [[Fraction(0) for _ in range(len(states) + 1)] for _ in states]
    for row, state in enumerate(states):
        for dest, rate in pair_transitions(state, N, beta, u, b):
            A[row][row] += rate
            if dest in idx:
                A[row][idx[dest]] -= rate
        i, j, k = state
        A[row][-1] = (i + k) * m[j + k] + (j + k) * m[i + k]
    for p in range(len(states)):
        pivot = next(q for q in range(p, len(states)) if A[q][p])
        A[p], A[pivot] = A[pivot], A[p]
        divisor = A[p][p]
        A[p] = [x / divisor for x in A[p]]
        for q in range(len(states)):
            if q != p:
                multiplier = A[q][p]
                A[q] = [x - multiplier * y for x, y in zip(A[q], A[p])]
    return A[idx[(0, 0, n)]][-1]


def labelled_pair_moments(N, beta, u, b, n):
    """Independent generator over labelled bitmask pairs, for small N only.

    Enumerates every parent, learner and two Bernoulli outcomes. This avoids
    assuming the count-level transition formulas and checks their lumping.
    """
    masks = range(1, 1 << N)
    states = list(product(masks, repeat=2))
    index = {s: j for j, s in enumerate(states)}
    M = len(states)
    Q = np.zeros((M, M))
    m, _ = single_moments(N, beta, u)
    rhs = np.empty(M)
    for row, (A, B) in enumerate(states):
        a, c = A.bit_count(), B.bit_count()
        rhs[row] = a * m[c] + c * m[a]
        transitions = {}
        for site in range(N):
            bit = 1 << site
            dest = (A & ~bit, B & ~bit)
            if dest != (A, B):
                transitions[dest] = transitions.get(dest, 0) + u
        if b > 0:
            for parent in range(N):
                for learner in range(N):
                    if parent == learner:
                        continue
                    for x, y in product((0, 1), repeat=2):
                        prob = (b if x else 1 - b) * (b if y else 1 - b)
                        dest = (A | ((1 << parent) if x and (A >> learner) & 1 else 0),
                                B | ((1 << parent) if y and (B >> learner) & 1 else 0))
                        if dest != (A, B) and prob:
                            transitions[dest] = transitions.get(dest, 0) + beta / (N * b) * prob
        else:
            # Limit b down to zero: successful arrows for the two traits
            # are independent, while the death clocks remain common.
            for parent in range(N):
                for learner in range(N):
                    if parent == learner:
                        continue
                    dest_A = (A | ((1 << parent) if (A >> learner) & 1 else 0), B)
                    dest_B = (A, B | ((1 << parent) if (B >> learner) & 1 else 0))
                    for dest in (dest_A, dest_B):
                        if dest != (A, B):
                            transitions[dest] = transitions.get(dest, 0) + beta / N
        for dest, rate in transitions.items():
            Q[row, row] -= rate
            if dest in index:
                Q[row, index[dest]] += rate
    cross = np.linalg.solve(-Q, rhs)
    return float(cross[index[((1 << n) - 1, (1 << n) - 1)]])


def stationary_forward_cycles(N, beta, u, mu, b, cycles, seed=20261102, shortfall=3):
    """Independent exact typed forward simulation, regenerative at empty culture.

    Every complete cycle comprises an empty wait and a nonempty excursion.
    Ratio-of-iid-cycle estimators have no burn-in parameter. Reported standard
    errors use the regenerative delta method; they are Monte Carlo errors.
    """
    rng = np.random.default_rng(seed)
    death_rate = N * u
    innovation_rate = N * mu
    contact_rate = beta * (N - 1) / b
    total = death_rate + innovation_rate + contact_rate
    rewards = np.zeros((cycles, 5))  # time, C, C^2, zero, (k-C)+
    next_trait = 0
    for cyc in range(cycles):
        traits = [set() for _ in range(N)]
        empty_wait = rng.exponential(1 / innovation_rate)
        rewards[cyc] = [empty_wait, 0, 0, empty_wait, shortfall * empty_wait]
        traits[int(rng.integers(N))].add(next_trait)
        next_trait += 1
        C = 1
        while C:
            dt = rng.exponential(1 / total)
            rewards[cyc] += dt * np.array([1, C, C * C, 0, max(shortfall - C, 0)])
            event = rng.random() * total
            if event < death_rate:
                traits[int(rng.integers(N))].clear()
            elif event < death_rate + innovation_rate:
                traits[int(rng.integers(N))].add(next_trait)
                next_trait += 1
            else:
                parent = int(rng.integers(N))
                learner = int(rng.integers(N - 1))
                learner += learner >= parent
                missing = traits[parent] - traits[learner]
                if b == 1:
                    traits[learner].update(missing)
                else:
                    traits[learner].update(t for t in missing if rng.random() < b)
            C = len(set().union(*traits))
    durations = rewards[:, 0]
    means = rewards.sum(axis=0) / durations.sum()
    influences = rewards - durations[:, None] * means
    se = influences.std(axis=0, ddof=1) / np.sqrt(cycles) / durations.mean()
    variance = means[2] - means[1] ** 2
    var_influence = influences[:, 2] - 2 * means[1] * influences[:, 1]
    var_se = var_influence.std(ddof=1) / np.sqrt(cycles) / durations.mean()
    return dict(N=N, beta=beta, u=u, mu=mu, b=b, cycles=cycles, seed=seed,
                total_time=float(durations.sum()), mean=float(means[1]),
                mean_se=float(se[1]), variance=float(variance), variance_se=float(var_se),
                p_zero=float(means[3]), p_zero_se=float(se[3]),
                shortfall_k=shortfall, expected_shortfall=float(means[4]),
                expected_shortfall_se=float(se[4]))
