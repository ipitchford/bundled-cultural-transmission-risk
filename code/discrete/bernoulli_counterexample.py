"""Exact variances when each newborn innovates at most one trait (Bernoulli(mu) innovation).

Case N=3, K=1, beta=9/10, mu=1/2, sample of one individual, innovations in generations -T..0.
With K=1 and n=1 the ancestry of the sample is a single potential-model path.
 Version 1: the path is fixed by G; the trait of the path's innovator at depth a reaches the
            sample with probability beta^a, independently across depths (own success marks).
 Version 2: each of the N innovators at depth a >= 1 lies on the trait's own random path with
            probability 1/N, so the depth-a contribution is a sum of N independent
            Bernoulli(mu beta^a / N); depth 0 is the sampled individual itself.
 Version 3: the shared path survives to depth a with probability beta^a (nested events).
Compare with mc_N3K1_b9o10_mu1o2_T3_bernoulli.txt (A={0}).
"""
from fractions import Fraction as Fr

N, BETA, MU, T = 3, Fr(9, 10), Fr(1, 2), 3
p = [MU * BETA ** a for a in range(T + 1)]
mean = sum(p)
v1 = sum(q * (1 - q) for q in p)
v2 = p[0] * (1 - p[0]) + sum(N * (q / N) * (1 - q / N) for q in p[1:])
cov = sum(MU * MU * BETA ** b * (1 - BETA ** a) for a in range(T + 1) for b in range(a + 1, T + 1))
v3 = v1 + 2 * cov
print(f"mean = {mean} = {float(mean):.6f} (all versions)")
print(f"Var V1 = {v1} = {float(v1):.6f}")
print(f"Var V2 = {v2} = {float(v2):.6f}")
print(f"Var V3 = {v3} = {float(v3):.6f}")
print("Ordering V2 <= V1 <= V3 holds:", v2 <= v1 <= v3)
