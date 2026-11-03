# Exact partial transmission: mathematical results and independent validation

Prepared 2 November 2026. The results below concern the explicitly defined stochastic cultural-transmission model. They do not assert empirical realism, first priority or a Research Excellence Framework grade. Poisson marking, conditional-expectation ordering and replica reward equations are established mathematical tools; their role and the closest model precedents must be credited in the paper.

## 1. Model, rates and scope

Let `V={1,...,N}` be a finite set of labelled individuals or population units. Individual `i` experiences independent resetting events at rate `delta_i>0`. A reset removes every trait held by that individual. Innovations occur at individual `i` as an independent Poisson process of rate `nu_i>=0`; every innovation is a previously unseen trait. Traits do not interact, compete for capacity or affect any event rate.

For each directed edge `j -> i`, fix an effective per-trait learning rate `beta_ij>=0`. At completeness parameter `b in (0,1]`, potential contacts on that edge occur at rate `beta_ij/b`. At a contact, each trait held by `j` is independently transmitted to `i` with probability `b`, independently across traits and contacts. The learner keeps its existing traits. Thus every individual trait has exactly the same transmission rate `beta_ij` for every `b`, while correlations between the histories of different traits change with `b`.

The homogeneous convention used by Wakano et al. (2026) is potential contact rate `s` per learner, uniformly choosing among all `N` units, including itself. In our homogeneous fixed-effective-rate comparison,

\[
\beta_{ij}=\beta/N\quad(i\ne j),\qquad \beta=sb.
\]

Here `beta` is a total rate including the conventional self-choice probability; the rate of an actual non-self edge is `beta/N`. Self-contact has no effect. This distinction matters for the two-person formulas below.

Fix a nonempty sampled set `A subset V` and observation time zero. Let `C_{A,b}` be the number of distinct traits held by at least one sampled individual at stationarity. A unique stationary construction exists by running the Poisson environment from the infinite past: for finite positive `b`, the finite ancestral-set process with resets is absorbed almost surely, erasing dependence on the initial trait configuration. All moments used below are finite.

The boundary `b=0` denotes the limiting model in which different traits have independent successful-contact clocks of rate `beta_ij`, while the resetting environment remains shared. It is not a finite contact rate with zero transmission probability. Conditional on the reset history, traits then evolve independently.

## 2. Weighted ancestry and the exact mixed-Poisson law

Let `D` denote the complete resetting history. Let `H` be a graph of successful learning arrows at rates `beta_ij`, independent of `D`. Trace a sample `A` backwards through `H,D`: a lineage at a learner branches to include its donor at a learning arrow, while a reset removes the corresponding ancestral node. Denote the ancestral set at backward time `t` by `B_A(t;H,D)`.

Define the innovation-weighted ancestral volume

\[
R_A(H,D)=\int_0^\infty\sum_{i\in B_A(t;H,D)}\nu_i\,dt.
\]

This is a random conditional mean, not the number of realised innovations. With common innovation rate `mu`, it is `mu` times total ancestral branch length.

Let `G_b` be the potential-contact graph at rates `beta_ij/b`. Conditional on `G_b`, let `H` be its independent `b`-thinning, so that its unconditional rates are `beta_ij`. Put

\[
M_{A,b}(G_b,D)=E[R_A(H,D)\mid G_b,D].
\]

For `b=1`, this is the ordinary weighted ancestry volume itself. For `b=0`, define `M_{A,0}(D)=E[R_A(H,D)|D]`.

**Theorem 1 (partial-transmission representation).**

\[
\boxed{C_{A,b}\mid G_b,D\ \sim\ \operatorname{Poisson}(M_{A,b})}
\qquad (0<b\leq1),
\]

and `C_{A,0}|D ~ Poisson(M_{A,0})` in the limiting independent-trait model. Consequently

\[
E C_{A,b}=E R_A,
\qquad
\operatorname{Var}(C_{A,b})=E R_A+\operatorname{Var}(M_{A,b}),
\]

\[
E[(C_{A,b})_k]=E[M_{A,b}^k]
\]

for each positive integer `k`, where `(x)_k=x(x-1)...(x-k+1)`.

**Proof.** Conditional on potential contacts and resets, innovations form a Poisson random measure on space-time. Give each innovation its own independent Bernoulli marks on the potential arrows. Its presence in the sample is determined by its own marked paths and the fixed reset history. Survival marks of different innovations are therefore independent conditional on `G_b,D`. Independent thinning of the innovation Poisson measure gives a Poisson number of retained innovations. Its mean is the integrated survival probability, which equals `M_{A,b}` by the definition of weighted ancestry and Tonelli's theorem. The moment formulas are those of a mixed Poisson law.

For integrability, the successful-arrow ancestral process on a finite state space is absorbed with an exponential tail. The bounded reward rate implies that `R_A` has all moments, and an exponential moment in a neighbourhood of zero. Conditional Jensen gives `E M_{A,b}^k <= E R_A^k` uniformly in `b`. This also justifies the stationary limiting constructions. The same argument applies conditional on resets to `b=0`. ∎

**Mean-spectrum consequence.** Every additive first-moment statistic of individual traits is invariant in `b` at fixed effective rates. This includes expected counts for each exact present-day carrier subset, the expected popularity spectrum and the expected age-popularity spectrum, whenever their means exist. A single trait's successful-arrow graph always has law `H`; only the dependence between different traits changes. Thus even the full collection of mean pattern counts cannot separate contact frequency from per-contact completeness.

## 3. A complete convex-order comparison

For integrable random variables `X,Y`, write `X <=cx Y` when `E phi(X)<=E phi(Y)` for every convex function for which the expectations are defined. Equal means are part of this ordering. It compares dispersion, and does not imply ordinary stochastic ordering or ordering of all tail probabilities.

**Theorem 2 (more complete contacts increase dispersion).** For fixed innovation, reset and effective learning rates, and any fixed sampled set,

\[
\boxed{0\leq b_1\leq b_2\leq1
\quad\Longrightarrow\quad
M_{A,b_1}\leq_{\rm cx}M_{A,b_2}
\quad\text{and}\quad
C_{A,b_1}\leq_{\rm cx}C_{A,b_2}.}
\]

**Proof.** First suppose `0<b1<b2<=1`. Independently on each edge, generate the successful graph `H` at rate `beta_ij`, add Poisson noise of rate `beta_ij(1/b2-1)` to obtain `G2`, and then add independent noise of rate `beta_ij(1/b1-1/b2)` to obtain `G1`. Retain the same independent reset graph `D` throughout.

Conditional on `G2`, the successful graph is a `b2`-thinning; conditional on `G1`, it is a `b1`-thinning. Moreover `H -> G2 -> G1` is a Markov channel. Therefore

\[
M_{A,b_1}=E[M_{A,b_2}\mid G_1,D].
\]

Conditional Jensen proves the first convex-order assertion. For `b1=0`, use `M_{A,0}=E[M_{A,b2}|D]`.

For a convex function on the nonnegative integers define `psi(x)=E phi(Poisson(x))`. It is convex: for suitable finite approximations,

\[
\psi''(x)=E[\phi(P_x+2)-2\phi(P_x+1)+\phi(P_x)]\geq0,
\]

where `P_x~Poisson(x)`. The general statement follows by standard convex approximation, after subtracting an affine lower bound. Apply the first ordering to `psi` and use Theorem 1. ∎

In particular,

\[
\operatorname{Var}(C_{A,b_1})\leq\operatorname{Var}(C_{A,b_2}),
\]

\[
\boxed{P(C_{A,b_1}=0)\leq P(C_{A,b_2}=0),}
\]

and, for every real threshold `k`,

\[
E(C_{A,b_1}-k)_+\leq E(C_{A,b_2}-k)_+,
\qquad
E(k-C_{A,b_1})_+\leq E(k-C_{A,b_2})_+.
\]

The zero-count inequality follows directly from the strictly convex function `exp(-M)`. The other two quantities are expected excess and expected shortfall. This theorem does not assert that `P(C<=k)` is ordered for every positive `k`, or that temporal collapse rates, durations or recovery times are ordered. The probability of an empty repertoire is a stationary occupancy probability.

### 3.1 Strictness

**Proposition 1.** Suppose there is a non-self edge `j -> i` with `beta_ij>0`, `nu_j>0` and `i in A`, and every reset rate is positive. Then for `0<=b1<b2<=1`,

\[
\operatorname{Var}(C_{A,b_1})<\operatorname{Var}(C_{A,b_2}),
\qquad
P(C_{A,b_1}=0)<P(C_{A,b_2}=0).
\]

**Proof.** In the coupling used above, consider a finite interval before observation. There is positive probability that an earlier subinterval contains a reset of every node and no potential contacts, eliminating older ancestry. Subsequently let the only relevant potential arrow be `j -> i`; let `j` reset after that arrow and `i` survive to observation. These requirements involve finitely many positive-probability restrictions on independent Poisson processes. Innovations located on a positive-length segment of `j` before the arrow can then reach `A` through that arrow, but cannot do so if it is absent.

Conditional on `G1,D`, the intermediate graph `G2` retains or deletes that arrow with probabilities strictly between zero and one. Retaining it increases `M_{A,b2}` by a strictly positive amount, since `nu_j>0` and the remaining `b2`-thinning has positive retention probability. Thus `M_{A,b2}` is not conditionally constant given `G1,D` on an event of positive probability. Strict conditional Jensen, applied to the square and to `exp(-x)`, proves both statements. At `b1=0`, the same event shows nonconstancy conditional on `D`, and the same argument applies. ∎

The condition is sufficient, not claimed necessary. It automatically holds for any nonempty sample in the homogeneous complete graph with `N>=2`, positive learning, positive innovation and positive reset rates. A longer accessible directed path can also give strictness, but this broader graph characterisation is unnecessary for the identification result and is not asserted here without a separate formal statement.

## 4. Exact second moments from two ancestral replicas

The following supplies a computable formula, rather than only an implicit mixing distribution. The two-replica method has antecedents in stochastic gene-transfer models; it is not a new general mathematical method.

Let `B1,B2` be two backward ancestral sets. Both see the same resets. On each edge they see a common successful-arrow process of rate `beta_ij b`, and separate successful-arrow processes of rate `beta_ij(1-b)` for each replica. Their marginal learning rates are therefore `beta_ij`. This is exactly the joint law obtained by independently `b`-thinning the same potential-contact graph twice.

Write

\[
g(B)=\sum_{i\in B}\nu_i,
\qquad m(B)=E_B\int_0^\infty g(B_t)\,dt.
\]

If `Q` is the single-set generator, then `-Qm=g`, with `m(empty)=0`. If `Q_b^(2)` is the pair-set generator, let

\[
F_b(B_1,B_2)=E\left[
\left(\int_0^\infty g(B_{1,t})\,dt\right)
\left(\int_0^\infty g(B_{2,t})\,dt\right)\right].
\]

First-step analysis gives the linear reward equation

\[
\boxed{-Q_b^{(2)}F_b(B_1,B_2)
=g(B_1)m(B_2)+g(B_2)m(B_1),}
\]

with `F_b=0` whenever either set is empty. Conditional independence of the two thinnings gives

\[
E[M_{A,b}^2]=F_b(A,A),
\]

and hence

\[
\boxed{\operatorname{Var}(C_{A,b})
=m(A)+F_b(A,A)-m(A)^2.}
\]

There is no closure approximation. For a general directed network the transient pair state space has `(2^N-1)^2` states, which limits direct computation. It is not appropriate to describe this general calculation as polynomial-time in population size.

### 4.1 Exact homogeneous lumping

Under exchangeability, use `(i,j,k)` for the numbers of nodes exclusive to replica 1, exclusive to replica 2, and in both. Put `n0=N-i-j-k` and `a_edge=beta/N`. The nonzero transitions are:

| Change in `(i,j,k)` | Rate |
|---|---:|
| `(i+1,j,k)` | `a_edge n0 [i+(1-b)k]` |
| `(i,j+1,k)` | `a_edge n0 [j+(1-b)k]` |
| `(i,j,k+1)` | `a_edge b n0 k` |
| `(i-1,j,k+1)` | `a_edge i(j+k)` |
| `(i,j-1,k+1)` | `a_edge j(i+k)` |
| `(i-1,j,k)` | `u i` |
| `(i,j-1,k)` | `u j` |
| `(i,j,k-1)` | `u k` |

The transient states satisfy `i+k>0`, `j+k>0`, `i+j+k<=N`. Their number is

\[
\binom{N+3}{3}-2(N+1)+1,
\]

which is cubic in `N`. A sample of size `n` starts at `(0,0,n)`. The reward right-hand side is

\[
\mu\{(i+k)m_{j+k}+(j+k)m_{i+k}\}.
\]

Here `m_n` solves the ordinary single-ancestry birth-death reward system with birth rate `beta n(N-n)/N`, death rate `u n` and reward `mu n`. The pair generator has only the eight listed transition types. Floating-point conditioning can deteriorate in strongly persistent regimes; a cubic number of states does not itself guarantee numerical precision.

## 5. Mean spectrum and structural identification of completeness

Write `P_k` for the number of traits held by exactly `k` of the `N` units. Set

\[
\eta=\mu/u,\qquad \tau=\beta/u=sb/u.
\]

An individual trait's carrier count is a birth-death process with birth rate `beta k(N-k)/N` and death rate `u k`, regardless of `b`. Stationary innovation flow gives

\[
E P_1=N\eta,
\qquad
\frac{E P_{k+1}}{E P_k}=
\frac{\tau k(N-k)}{N(k+1)}.
\]

Equivalently,

\[
E P_k=\frac{N\eta}{k}
\left(\frac{\tau}{N}\right)^{k-1}
(N-1)^{\underline{k-1}}.
\]

These means are existing single-trait spectrum calculations, presented to show what the second moment adds.

**Corollary 3 (structural identification).** In the homogeneous model with known `N>=2` and positive innovation, reset and learning rates, the following population quantities identify `(eta,tau,b)`:

\[
\big(E P_1,E P_2,\operatorname{Var}(C_{V,b})\big).
\]

Indeed,

\[
\eta=\frac{E P_1}{N},\qquad
\tau=\frac{2N}{N-1}\frac{E P_2}{E P_1}.
\]

For fixed `N,eta,tau`, Proposition 1 makes total-repertoire variance strictly increasing in `b`. It is continuous, as follows either from the finite pair generator and its invertible transient resolvent or from the coupled representation. Consequently the variance determines at most one `b`; a value outside the model's variance range rejects the model. The relative encounter rate follows as `s/u=tau/b` for positive `b`.

No absolute time scale is identified: multiplying `mu,u,s` by one common positive constant leaves the stationary law unchanged. The result requires population expectations and a stationary variance, not merely one realised spectrum. Replicated comparable populations, or an adequate stationary time series with ergodic averaging, would be needed to estimate those quantities. It does not supply a finite-sample estimator, confidence interval or robustness to population heterogeneity. Identification becomes weak in practice when innovation or effective learning is small, because the variance contrast is then small.

## 6. An explicit two-person formula and inverse

For `N=2`, denote the effective rate of each directed edge by `lambda_edge=beta/2`. Set

\[
\rho=\lambda_{\rm edge}/u=\tau/2,
\qquad \eta=\mu/u.
\]

Then

\[
\boxed{E C_{V,b}=\eta(2+\rho),}
\]

\[
\boxed{\operatorname{Var}(C_{V,b})
=\eta(2+\rho)+
\eta^2\frac{\rho^3+(8-2b)\rho^2+10\rho+4}
{2+(3-2b)\rho}.}
\]

The derivative is explicitly positive:

\[
\boxed{\frac{\partial}{\partial b}\operatorname{Var}(C_{V,b})
=\frac{2\eta^2\rho(\rho+1)(\rho+2)^2}
{[2+(3-2b)\rho]^2}>0.}
\]

For complete copying,

\[
\operatorname{Var}(C_{V,1})
=\eta(2+\rho)+\eta^2(\rho^2+4\rho+2).
\]

The difference from that endpoint is

\[
\operatorname{Var}(C_{V,1})-\operatorname{Var}(C_{V,b})
=\frac{2\eta^2(1-b)\rho(\rho+1)(\rho+2)}
{2+(3-2b)\rho}.
\]

Let

\[
V=\frac{\operatorname{Var}(C)-E C}{\eta^2},
\qquad A=\rho^3+8\rho^2+10\rho+4.
\]

The inverse is

\[
\boxed{b=\frac{V(2+3\rho)-A}{2\rho(V-\rho)}.}
\]

Its denominator is positive throughout the model because

\[
V-\rho=\frac{(\rho+1)(\rho+2)^2}{2+(3-2b)\rho}>0.
\]

### Derivation by direct moment equations

Rescale time to `u=1`, and let `X,Y,Z` count traits held only by person 1, only by person 2, and by both. A contact from person 1 to person 2 moves `Binomial(X,b)` traits from `X` to `Z`, at rate `rho/b`. The reverse event acts similarly on `Y`. A reset of person 1 maps `(X,Y,Z)` to `(0,Y+Z,0)`; the other reset is symmetric. Innovations add one to `X` or `Y` at rate `eta` each.

The stationary first moments are `E X=E Y=eta`, `E Z=rho eta`. With

\[
f=E[X(X-1)],\quad v=E[XY],\quad w=E[XZ],\quad g=E[Z(Z-1)],
\]

the generator gives

\[
v=\frac{\eta^2}{1+\rho},
\]

\[
[1+2\rho(1-b)]f-2(1+\rho)w=2\eta^2,
\]

\[
(2+\rho)w=\rho\eta^2+\rho v+\rho(1-b)f,
\qquad g=\rho(2w+bf).
\]

Solving,

\[
f=\frac{2\eta^2(1+\rho)(2+\rho)}{2+(3-2b)\rho},
\qquad
w=\frac{\rho\eta^2}{1+\rho}
+\frac{\rho(1-b)}{2+\rho}f.
\]

Since `E[C(C-1)]=2f+g+2v+4w`, substitution yields the displayed variance. These equations are also valid in the `b=0` limit.

## 7. Boundary limits and what is not established

The `b downarrow 0` stationary law is the conditional-on-resets independent-trait model. One direct proof first fixes a finite observation window. Conditional on its finite Poisson number of innovations, simultaneous transmissions of two specified traits occur at rate proportional to `b`, while each individual rate stays fixed. The finite trait process therefore converges to independent successful-arrow processes sharing resets. The expected number of surviving innovations older than the window is independent of `b` and tends to zero as the window expands, because each trait has the same exponentially decaying finite-population survival law. This permits passage to stationarity. Uniform moment bounds follow from Theorem 2 and the finite exponential moment at `b=1`.

No large-population limit, phase transition or singularity at `b=1` is proved here. Apparent concentration at any fixed `b<1` would require a separate asymptotic argument controlling correlations under the shared reset environment. The finite-population theorems stand without that argument. Numerical results must not be promoted to such a theorem.

## 8. Independent implementation checks

`agent_math/pair_ancestry.py` implements the backward reward equations independently of the root agent's forward count-moment code. It compares the exchangeable triple-state generator with a separate fully labelled two-set generator, and checks the two-person rational formula and its inverse.

The executed checks covered 15 two-person parameter cases, 18 labelled-versus-lumped comparisons at `N=2,3,4` and sample sizes one and the full population, and a directed heterogeneous three-node example. The largest absolute error in the two-person variance formula was `1.14e-13`; the largest discrepancy between labelled and lumped variances was `5.56e-16`.

For the heterogeneous example, the mean was constant at `1.14690787845513`. Its variance increased from `1.6332430179089623` at `b=0` to `1.947440465105384` at `b=1`, with all tested intermediate values increasing. This is a check of the implementation and a finite illustration, not the proof of the ordering.

The machine-readable output is `agent_math/independent_pair_validation.json`. Numerical linear-system checks do not establish priority or empirical applicability.

## 9. A nonnegative power series for the variance

This is an application of established Poisson-space isometry and Mehler thinning identities, not a new general probability theorem. It gives a stronger shape restriction than monotonicity for this model.

Write `V_A(b)=Var(C_{A,b})` and `m_D=E_H R_A(H,D)`. Let `Lambda` be the intensity measure of successful arrows: counting measure on directed edges, weighted by their effective rates, times Lebesgue measure on backward time. For a fixed reset history `D`, define the add-one operator

\[
\Delta_xR_D(H)=R_A(H+\delta_x,D)-R_A(H,D),
\]

and let `Delta^k` denote its iteration. Since `E R_A^2<infinity`, for almost every `D` the functional `R_D` is square integrable under the Poisson arrow law. Define

\[
T_{k,D}(x_1,\ldots,x_k)
=E_H[\Delta_{x_1}\cdots\Delta_{x_k}R_D(H)],
\qquad
a_k=\frac1{k!}E_D\|T_{k,D}\|_{L^2(\Lambda^k)}^2.
\]

The kernels are well defined almost everywhere by the Poisson isometry; conditional expectations can be taken as jointly measurable versions. All `a_k` are nonnegative.

**Corollary 4 (absolute monotonicity).** For `0<=b<=1`,

\[
\boxed{V_A(b)=V_A(0)+\sum_{k=1}^{\infty}a_k b^k,}
\qquad
V_A(0)=E R_A+\operatorname{Var}_D(m_D),
\]

\[
\sum_{k=1}^{\infty}a_k=V_A(1)-V_A(0)<\infty.
\]

Thus `V_A` is analytic for `|b|<1`, all derivatives are nonnegative for `0<=b<1`, and

\[
\boxed{V_A(0)+a_1 b\ \leq\ V_A(b)
\ \leq\ (1-b)V_A(0)+bV_A(1).}
\]

In particular, the variance is convex. Under Proposition 1's positive-edge condition, `a_1>0`, so `V_A'(b)>=a_1>0` for every `0<=b<1`.

**Proof.** Conditional on `D`, generate two independent `b`-thinnings `H1,H2` of the same potential graph. Their common arrows have intensity `b Lambda`, and each has private arrows of intensity `(1-b)Lambda`. Conditional on `H1`, the second graph is a `b`-thinning of `H1` plus an independent Poisson graph of intensity `(1-b)Lambda`. This is the Poisson Mehler operator `P_b`. Therefore

\[
E[M_{A,b}^2\mid D]=E_H[R_D(H)(P_bR_D)(H)].
\]

The Fock-space isometry of Last and Penrose (2011, Theorem 1.1, equation (1.5)) and the Mehler kernel identity of Last, Peccati and Schulte (2016, equation (3.6)) give

\[
T_k(P_bR_D)=b^k T_k(R_D),
\]

\[
E_H[R_D P_bR_D]
=m_D^2+\sum_{k\geq1}\frac{b^k}{k!}\|T_{k,D}\|^2.
\]

Integrating over `D` and applying Theorem 1 yields the series. Tonelli's theorem applies because its terms are nonnegative. At `b=1`, the isometry gives the summability assertion. Derivative and chord inequalities follow from the nonnegative power-series coefficients.

For strictness, the first coefficient is

\[
a_1=E_D\int
\{E_H[\Delta_xR_D(H)]\}^2\,\Lambda(dx).
\]

Adding a learning arrow can only enlarge the ancestral set, so its add-one cost is nonnegative. Under the positive-edge condition, select a nontrivial interval of possible arrow times `x=(j->i,t)`. There is positive probability of a reset history on which `i` survives from observation back through that interval; `j` has a recent reset before that interval and a further reset after it; and all sampled nodes have reset by a later finite time. Conditional on such a reset history, the event of no other successful arrows in this finite window has positive probability. Adding `j->i` then creates a positive-length ancestral segment at `j`, weighted by `nu_j>0`, which is absent without that arrow. The mean add-one cost is therefore positive on a set of positive `P_D times Lambda` measure, proving `a_1>0`. ∎

**Finite-model endpoint refinement.** The power-series argument alone does not bound derivatives at `b=1` for an arbitrary square-integrable Poisson functional. Here, however, the finite pair-set generator in Section 4 has entries affine in `b`, and its transient submatrix is invertible throughout `[0,1]`. Invertibility follows, for example, from a positive-probability block with no arrows and at least one reset at every node, which clears both replicas. Its finite-dimensional resolvent, and hence `V_A(b)`, is rational and analytic in a neighbourhood of each point of `[0,1]`. All derivatives therefore extend finitely to the endpoints and remain nonnegative there. No population-size-uniform derivative bound follows.

The decomposition separates two sources of excess variance. `Var_D(m_D)` remains when different traits have independent successful-contact clocks, because the reset history is shared. The positive series measures the additional variance due to contacts shared across traits. It does not imply independent replicas at `b=0` before conditioning on resets.

### Primary references for the Poisson identities

- Last, G., and Penrose, M. D. (2011). *Poisson process Fock space representation, chaos expansion and covariance inequalities*. Probability Theory and Related Fields **150**(3–4), 663–690. DOI: [10.1007/s00440-010-0288-5](https://doi.org/10.1007/s00440-010-0288-5). [Author manuscript](https://arxiv.org/pdf/0909.3205). Theorem 1.1 and equations (1.3)–(1.5), manuscript page 2, supply the add-one kernels and isometry.
- Last, G., Peccati, G., and Schulte, M. (2016). *Normal approximation on Poisson spaces: Mehler’s formula, second order Poincaré inequalities and stabilization*. Probability Theory and Related Fields **165**(3–4), 667–723. DOI: [10.1007/s00440-015-0643-7](https://doi.org/10.1007/s00440-015-0643-7). [Author manuscript](https://arxiv.org/pdf/1401.7568). Section 3, equation (3.1), Lemma 3.1 and equation (3.6), manuscript pages 14–15, supply the thinning/immigration operator and its action on the kernels.

## 10. A precise temporal consequence for complete population loss

This extension concerns the empty state of the **entire finite population**, not a sampled repertoire or partial loss. Let

\[
\alpha=\sum_i\nu_i>0,
\qquad
\pi_b=P(C_{V,b}=0).
\]

Whenever the population is empty, only a fresh innovation can make it nonempty. The waiting time `J` to this departure is exponential with rate `alpha`. Let `B_b` denote the positive excursion duration after that first innovation until the next completely empty state. Innovations continue during the excursion. The location of the initial innovation has probabilities `nu_i/alpha`. By the strong Markov property, successive empty holding times and positive excursions form independent identically distributed cycles, with `J` independent of `B_b`.

**Corollary 5 (complete-loss recurrence and positive-excursion transforms).** In the stationary finite model,

\[
\boxed{\text{frequency of entries into the empty state}=\alpha\pi_b,}
\]

\[
\boxed{E B_b=\frac{1-\pi_b}{\alpha\pi_b},
\qquad E(J+B_b)=\frac1{\alpha\pi_b}.}
\]

Hence the frequency of complete-empty entries is nondecreasing in `b`, while the mean positive excursion and the mean cycle duration are nonincreasing. Under Proposition 1's condition for `A=V`, these inequalities are strict.

Moreover, for `0<=b1<=b2<=1` and every `z>0`,

\[
\boxed{E e^{-zB_{b_1}}\leq E e^{-zB_{b_2}}.}
\]

This is Laplace-transform order; it is not a claim of ordinary stochastic ordering of excursion durations or of ordered time-dependent hazards.

**Proof of finite regeneration times.** For `b>0`, a fixed finite time block with no innovation, no potential contacts and at least one reset at every node has a positive probability independent of the current finite trait configuration. Such a block forces emptiness. Repeated disjoint blocks give a geometric bound on the time to empty and in particular a finite mean. At `b=0`, do not invoke a no-potential-contact event, because that representation has infinite potential-contact rate. Instead use the countable Markov chain of numbers of traits in each nonempty carrier subset. It is nonexplosive: the number of distinct traits can increase only through the finite-rate innovation process, and the total rate of individual-trait moves is bounded by a constant times that number. The already constructed stationary distribution puts strictly positive mass `E exp(-M_{V,0})` on the empty state. A state having positive mass under an invariant probability distribution is positive recurrent. Restricting to the communicating class of empty gives the finite-mean regenerative cycles. Empty is accessible from every finite configuration by a finite reset sequence with no innovations and no successful transmissions of its finitely many traits.

The renewal-reward fraction of time empty is

\[
\pi_b=\frac{EJ}{EJ+EB_b}
=\frac{1}{1+\alpha EB_b}.
\]

One entry occurs per complete cycle, giving the frequency and mean formulas.

For the transform ordering, start the population empty at time zero and write

\[
p_{00,b}(t)=P(C_{V,b}(t)=0\mid C_{V,b}(0)=0).
\]

Only innovations in `[0,t]` contribute. Apply the same conditional Poisson construction and Markov-channel comparison as in Theorems 1–2 to the ancestry truncated at elapsed time `t`. For every `t>=0`,

\[
p_{00,b_1}(t)\leq p_{00,b_2}(t).
\]

Let `F_b` be the distribution of `B_b` and put `\widehat p_b(z)=integral_0^infinity exp(-zt)p_{00,b}(t)dt`. Either the first empty holding interval has not ended by `t`, or a complete first cycle has ended and the process restarts from empty. Thus

\[
p_{00,b}(t)=e^{-\alpha t}
+\int_{[0,t]}p_{00,b}(t-s)\,d(G_\alpha*F_b)(s),
\]

where `G_alpha` is the exponential waiting-time distribution. Taking Laplace transforms yields

\[
\boxed{\widehat p_b(z)
=\frac{1}{z+\alpha-\alpha E e^{-zB_b}}.}
\]

Pointwise ordering of `p00` gives ordering of its Laplace transform. Since the displayed expression is strictly increasing in `E exp(-zB_b)`, the asserted transform order follows. Under the positive-edge condition, the finite-window single-arrow witness fits into any interval of length `t>0` and gives strict `p00` inequalities; the transform inequality is then strict for each `z>0`. ∎

An entry here means a transition from a nonempty population to a completely empty one. Reset or contact events that occur while the population is already empty are not additional entries. The result assumes the stated stationary or empty-initial experiment; it does not compare collapse from arbitrary initial repertoires, partial-threshold losses, recovery to a nonzero repertoire threshold, or arbitrary sample-empty episodes. Empty holding times themselves have the same exponential law for every `b`.

