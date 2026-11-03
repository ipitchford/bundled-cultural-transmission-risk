# Independent sceptical checks before manuscript review

**Date:** 2 November 2026. **Status:** working mathematical and interpretive audit; not the final eleven-part academic review. Materials inspected include `target.yaml`, `code/transmission.py`, `math_results.md`, the numerical outputs, `novelty_note.md`, and the closest primary cultural papers specified below. This note records independent checks and cautions while the manuscript is being prepared.

## 1. Core representation and coupling: passed

Conditional on potential contacts and deaths, each unique innovation can be assigned its own independent Bernoulli marks on every contact. Whether that innovation reaches the observation sample is then an independent mark of its innovation point. Poisson marking therefore yields the proposed mixed-Poisson count. The representation requires independent, unique innovations and trait-specific copying marks; it would not follow for intrinsically coupled trait acquisition or recurrent innovation identities.

For `0<b1<b2`, a correct construction is `H -> G2 -> G1`: start with the effective graph, add independent Poisson noise to obtain the potential graph at `b2`, then add further independent noise to obtain the graph at `b1`. These form a Markov channel. Conditional independence gives

\[
E\{E[R(H,D)\mid G_2,D]\mid G_1,D\}
=E[R(H,D)\mid G_1,D].
\]

The observed graph sigma-fields should **not** simply be asserted to be nested: knowing the superposed graph does not reveal which arrows belong to its thinned subgraph. The Markov conditional-independence argument is the necessary justification. The revised `math_results.md` now includes it.

Conditional Jensen gives the convex order of mixing intensities; Poissonisation preserves convex order. These are standard probability principles. The potentially original unit is their exact application to this cultural-transmission family and its model-specific consequences, not a new general stochastic-order theorem.

## 2. Strictness: a positive-probability witness

For a homogeneous complete graph with `N>=2`, choose a sampled learner and a different donor. On a finite backward-time window, require all other sampled ancestors to die before a chosen contact time; require a single donor-to-learner contact, followed backwards by the learner's death and then the donor's older death. No other contacts occur, and the window clears all possible ancestors before its older endpoint. These finite Poisson restrictions have positive probability.

Conditioning on the larger potential graph and death history, the intermediate graph retains this sole useful arrow with probability `b1/b2`. If absent, the mixing volume is a baseline `M0`; if present, it is `M0+b2*nu_donor*Delta`, with `Delta>0`. Hence the conditional variance is positive on an event of positive probability. Strict Jensen for the square and for `exp(-x)` proves strict count variance and strict stationary zero-count probability. The same idea works for the sufficient directed-edge condition stated in `math_results.md`.

The claim correctly fails in degenerate cases such as one individual, no effective learning, or no innovations. The zero-copying boundary removes shared-contact dependence but leaves the shared death environment. It does not generally produce an ordinary Poisson count or unconditional independence among traits.

## 3. Structural identification: passed with the stated observational boundary

Let `Yk` count traits held by exactly `k` of all `N` individuals. The carrier birth and death rates are

\[
\beta k(N-k)/N,\qquad uk.
\]

Innovation balance and zero stationary flux between adjacent positive carrier counts give

\[
E Y_1=N\eta,\qquad E Y_2=\eta\tau(N-1)/2,
\quad \eta=\mu/u,\ \tau=\beta/u.
\]

Thus the first two full-population spectrum means identify `eta,tau` when `N>=2` is known and the rates are nondegenerate. Strict variance monotonicity then identifies at most one `b`, and a variance outside the feasible range rejects the model. Absolute time scale is not identified.

This is a population-law statement. One realised spectrum is not its expectation, and variance across individuals within one population is not the stationary variance of the population's union repertoire. Replicated comparable populations or appropriately sampled stationary temporal data would be needed. No finite-sample estimation or confidence procedure has yet been established. The same inversion should not be claimed for an arbitrary subsample spectrum without a separate argument.

## 4. New temporal consequence: passed for the full population

Write `alpha=sum_i nu_i>0` and let `pi0_b` be stationary probability that the **entire population** has no traits. The full empty configuration is a regenerative state. The scalar repertoire count need not itself be Markov for this argument.

Let `J` be the time spent empty and `B_b` the next strictly positive repertoire excursion. Then `J~Exp(alpha)`, the initiating innovation occurs at site `i` with probability `nu_i/alpha`, and future independent Poisson increments make `J` independent of `B_b`. Innovations continue during `B_b`; it is not the extinction time of the initiating trait alone. For positive `b`, finite mean return time follows from a positive-probability block with no contacts or innovations and at least one reset at every site, which clears any initial configuration.

Renewal reward gives

\[
\pi_{0,b}=\frac{1}{1+\alpha E B_b},\quad
\text{empty-entry frequency}=\alpha\pi_{0,b},\quad
E(J+B_b)=\frac{1}{\alpha\pi_{0,b}}.
\]

Thus the stationary zero-count ordering also orders a genuine temporal frequency: more complete contacts, compensated by fewer contacts, give more entries into complete absence and shorter mean nonempty episodes. Entries mean transitions into absence separated by nonempty excursions, not resets occurring while already empty.

Starting empty at time zero, the finite-horizon count has the same conditional-Poisson/garbling representation with a truncated ancestral integral. Consequently `p00_b(t)=P(C_b(t)=0 | empty at 0)` increases in `b`. The initial empty condition matters; this proof does not automatically cover arbitrary deterministic nonempty initial repertoires.

The alternating-renewal equation gives, for `z>0`,

\[
\widehat p_{00,b}(z)
=\frac{1}{z+\alpha-\alpha E e^{-zB_b}}.
\]

Since the left side increases in `b`, so does `E exp(-zB_b)`. Therefore the positive episode duration decreases in **Laplace-transform order**. This is not ordinary stochastic order, an ordered first-passage law from arbitrary states, or an ordered instantaneous collapse hazard. The regenerative identities themselves are standard; their application is an additional consequence of the model-specific comparison. A sample's empty state is generally not regenerative when unsampled carriers can reintroduce traits.

## 5. Code and exact checks

The eight lumped pair-generator transitions are consistent with common deaths, independent trait-specific successes and the probability `b^2` of a common success. The reward right-hand side correctly accounts for the remaining reward of one replica after the other becomes empty.

Independent recalculation for `N=2,beta=2,u=1` gives:

| Sample size | `b=1/4` | `b=1/2` | `b=1` |
|---:|---:|---:|---:|
| 1: `E[W^2]` | 8 | 17/2 | 10 |
| 2: `E[W^2]` | 14 | 29/2 | 16 |

These rational calculations agree with both the sparse lumped and independently enumerated labelled generators. The direct two-person forward moment equations in `math_results.md` also check out, including the variance derivative and inverse formula.

An initially observed helper bug at `b=0` omitted contacts in the labelled generator. It was reported and corrected during the review. The corrected helper now implements independent successful arrows for the two replicas while retaining common deaths. For `N=2,beta=2,u=1,n=2,b=0`, both generators give `E[W^2]=13.6`. This resolved issue should not be reported as an outstanding defect.

Large-population floating-point solutions require appropriately rounded presentation. A small linear-system residual alone is not an accuracy bound in an ill-conditioned persistent regime. The available endpoint and small exact checks are useful; no asymptotic theorem follows from the large-N table.

## 6. Primary literature independently checked and early significance assessment

I independently opened Wakano et al. (2026), especially §4 and equations (36)–(37), and Wakano and Aoki (2024), §§1.6.1–1.6.2. The former specifies the compensating `s/b` comparison and explores partial copying numerically. The latter explicitly identifies distributions beyond means and correlated genealogies as work remaining in its discrete-generation setting. The present continuous-time theorem answers a precise related model question; it does not solve that complete discrete-generation problem or give a closed-form probability mass function in general.

I also independently opened Baumdicker and Pfaffelhuber (2014), Remark 2.4 and §6.3. Remark 2.4 explicitly says that holding the products of event rates and per-gene transfer/loss probabilities fixed preserves first moments while second moments differ. Thus “same means, different variance” is established prior art. It does not state the directional comparison checked here. The exact order, strict identification and regenerative consequences are the additional units requiring assessment. The replica and overlap-count method is also established in that paper.

The current novelty note documents the further generic probability predecessors. No contrary exact cultural comparison was located in the inspected primary material, but this is bounded negative evidence, not proof of priority. The latest Scheie et al. full text and the 2018 founding paper remain access-limited checks documented by the literature auditor.

**Early reviewer judgement:** the result is worth developing if the model-specific originality survives external scrutiny and the final exposition makes the general network, strictness, identifiability and temporal consequences clear. The theoretical significance can be substantial without adding an empirical study, but model assumptions and observation requirements must remain explicit. Stationary repertoire size does not measure cultural value, and complete absence in this idealised model is not a demonstrated historical collapse mechanism.

### Primary sources consulted for this note

Baumdicker, F., & Pfaffelhuber, P. (2014). The infinitely many genes model with horizontal gene transfer. *Electronic Journal of Probability, 19*, Article 115, 1–27. https://doi.org/10.1214/EJP.v19-2642 — full text https://arxiv.org/html/1301.6547v2

Wakano, J. Y., & Aoki, K. (2024). Forward and backward modeling of cultural evolutionary processes. *Proceedings of the National Academy of Sciences, 121*(48), e2322880121. https://doi.org/10.1073/pnas.2322880121 — full text https://pmc.ncbi.nlm.nih.gov/articles/PMC11621620/

Wakano, J. Y., Ohtsuki, H., Kobayashi, Y., & Baake, E. (2026). Evolving genealogies in cultural evolution, the descendant process, and the number of cultural traits. *Theoretical Population Biology, 168*, 1–18. https://doi.org/10.1016/j.tpb.2025.12.003 — full text https://arxiv.org/html/2505.08677v2
