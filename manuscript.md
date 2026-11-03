---
title: "Equal mean diversity, unequal risk"
subtitle: "An exact comparison of partial and bundled cultural transmission"
author: "Anonymous"
date: "3 November 2026"
---

## Abstract

When people learn from one another, a single encounter can pass on many pieces of knowledge at once, or each piece can travel separately. Does this bundling matter if every individual trait is passed on at the same average rate? We answer the question exactly for the continuous-time model of cultural trait accumulation studied by Wakano, Ohtsuki, Kobayashi and Baake (2026). In that model individuals invent traits, learn from one another and die. At each encounter every trait is copied independently with probability $b$, and encounters occur $1/b$ times as often, so the per-trait transmission rate is fixed. We prove that the expected number of traits, and the expected number of traits held by exactly $k$ individuals, do not depend on $b$. However, every such count becomes more spread out as $b$ increases, in the precise sense of *convex order*. Here, spread means variation across possible population histories. As a result, the variance, the expected shortfall below any target repertoire size and the chance of holding no traits at all all increase with bundling. The same ordering holds when forgetting also comes in bundles. With Poisson innovation, it also holds for the discrete-generation versions 1–3 for which Wakano and Aoki (2024) posed an open distributional question. We characterise exactly when the increase is strict and give exact moment equations. In supercritical, well-connected populations the effect is large. For example, at the population size of 100 used in the source simulations, copying 93% of a donor's traits per encounter leaves a third of the whole-repertoire variance. The proofs combine established probability tools; the contribution lies in the model-specific orderings, their strictness and their consequences. No empirical data are used.

**In brief.** One big lesson exposes many traits to the same fortunes, while many small lessons spread those fortunes out. On average the two arrangements preserve the same amount of culture. Bundled learning makes good and bad histories more extreme, including the history in which everything is lost.

## 1. Introduction

### 1.1 The question

Models of cultural evolution often count traits: how many distinct techniques, words or designs a population holds, and how widely each is shared (Strimling et al., 2009; Lehmann et al., 2011; Aguilar & Ghirlanda, 2015; Nakamura et al., 2020). A trait is lost when everyone who knows it dies or forgets it before passing it on, and debates about cultural loss turn on such events (Henrich, 2004; Kolodny et al., 2015).

Learning, however, does not happen trait by trait. A learner who spends time with a teacher may pick up many things at once. Traits that travel through the same encounters share their fortunes. They spread together when the encounters happen, and they are lost together when the encounters do not happen or the carriers die. A model that tracks each trait's average transmission rate cannot see this dependence.

This paper asks a sharp version of the question. Suppose we hold fixed the rate at which each individual trait is transmitted, and change only how bundled the transmission is. What changes, and what stays the same?

### 1.2 The answer in one paragraph

Means stay the same; spread and risk increase with bundling. Section 4 proves that every count of the form "number of traits whose carriers have property X" has the same expected value under every degree of bundling. Such counts include the repertoire of any group, the number of traits known to exactly $k$ people and the number known to at least $k$ people. Each count is larger in convex order when bundling is greater: more spread out around the same mean (Section 4 defines this precisely). Figure 1 shows the mechanism, and Section 7 shows that the size of the effect can be very large.

![The mechanism. Left: with bundled copying ($b=1$), encounters are rare and each passes every trait the donor holds, so the two traits travel through the same encounters and share their fate. Right: with partial copying ($b=1/4$), encounters are four times as frequent and each passes a given trait with probability 1/4. The two traits now mostly travel through different encounters. A cross marks a reset (death), which removes every trait at that site under both arrangements.](figures/schematic.png){width=100%}

### 1.3 What is held fixed

Let $s$ be the rate of encounters per learner and $b$ the probability that a particular trait is copied at an encounter. We compare arrangements with the same product
$$\beta=sb,$$
the per-trait effective learning rate. Lower $b$ is compensated by more encounters. A single trait, viewed on its own, has exactly the same transmission process under every arrangement; what changes is how often different traits share the same opportunities. This is the comparison used by Wakano et al. (2026, §4), who replace their encounter rate $s$ by $s/b$. **Notation dictionary:** our $\beta$ is their $s$, and our encounter rate $\beta/b$ is their $s/b$. The comparison does not hold teaching time, attention or encounter costs fixed, and it does not advise copying less accurately. Remark 1 treats the other natural comparison, a fixed encounter rate.

### 1.4 Contributions

1. **Equal means, ordered spread (Theorems 1–2, Corollaries 1–2).** Given the shared history of encounters and deaths, every trait count is Poisson with a random mean. These random means, and hence the counts, increase in convex order with $b$. Variance, every expected shortfall and the probability of an empty repertoire are therefore nondecreasing. The covariance matrix of the popularity spectrum increases in the positive-semidefinite order.
2. **Exactly when the effect is strict (Proposition 1).** For a group $A$, the increase is strict if and only if some inventing site has a directed chain of learning links of length at least one into $A$. Otherwise the law does not depend on $b$ at all.
3. **Bundled forgetting (Theorem 3).** If deaths are replaced by loss events that erase each trait with probability $d$, the same ordering holds in $b$ and in $d$ separately. When both are fully unbundled the count is exactly Poisson.
4. **The original discrete-generation question (Theorem 4).** In the discrete-generation model of Kobayashi et al. (2018), with Poisson innovation, Versions 2, 1 and 3 of Wakano and Aoki (2024) are ordered in convex order at equal means, for every sample and every spectrum class.
5. **Size of the effect (Section 7).** Exact computations show a sharp contrast between regimes. Below the critical learning rate, bundling changes variance by at most about a quarter. Above it, the variance ratio between partial and whole-repertoire copying collapses as the population grows, and most of the change happens close to $b=1$.
6. **Further consequences.** We give exact temporal consequences for complete-loss episodes of the whole population (Section 8), exact moment equations and a two-person closed form (Sections 9 and 11), and an identification result: unlike the expected spectrum, repertoire variance carries information about bundling (Section 10).

### 1.5 Relation to prior work

The general tools are established, and we credit them where they are used. They are Poisson marking, conditional Jensen ("information garbling"), the lifting of convex order to Cox processes (Błaszczyszyn & Yogeshwaran, 2009), Poisson chaos and Mehler formulas (Last & Penrose, 2011; Last et al., 2016), renewal identities and the association inequalities of Esary et al. (1967) and Christofides and Vaggelatou (2004). Once a conditional Poisson representation is available, the ordering in Theorem 2 takes a few lines. The model-specific content lies elsewhere: in establishing that representation for these genealogies, in the strictness characterisation, in the extensions to every spectrum class, to bundled forgetting and to the discrete versions, and in the computations.

The phenomenon itself has close relatives.

- **Gene transfer.** Baumdicker and Pfaffelhuber (2014, Remark 2.4) note that, in the infinitely-many-genes model, transferring or losing several genes in one event preserves first moments but changes second moments. Their §6.3 computes variances from paired genealogies, as our Section 9 does.
- **Population genetics.** Recombination plays the role of partial copying there. Linked sites share one genealogy, so the number of segregating sites is mixed Poisson with excess variance (Watterson, 1975), and recombination reduces that variance at fixed mean (Hudson, 1983; Pluzhnikov & Donnelly, 1996). Recombination does not literally garble one common history, so a convex ordering is not automatic there.
- **Cultural evolution.** Wakano and Aoki (2024, §§1.6.1–1.6.2) anticipate the cultural phenomenon. Without proof, they remark that a version in which all traits share the same cultural parents should fluctuate greatly, and that package transmission may increase the variance of the cultural spectrum. They list distributions beyond means, in the presence of correlated genealogies, as an unresolved task. Our theorems give this a rigorous form: in the continuous-time model that Wakano et al. (2026, §4 and Figure 16) explored by simulation, and in the discrete model itself (Section 6).
- **Linked traits.** Yeh et al. (2019) study the transmission of linked trait packages in a model with a different intervention.

Scheie et al. (2026) study sequential accumulation of independent traits; their published abstract concerns Gaussian limits of the popularity spectrum, not comparisons across bundling.

### 1.6 Roadmap

Section 2 states the model and notation. Section 3 explains why a trait count is conditionally Poisson. Section 4 proves the main ordering and its strict form. Section 5 adds bundled forgetting, and Section 6 treats the discrete-generation versions. Section 7 measures how large the effect is. Sections 8–11 cover complete-loss episodes, exact moment equations, identification and a two-person example with independent checks. Section 12 discusses interpretation and limits. Longer proofs are in the appendices.

## 2. The model

There are $N$ labelled sites, representing individuals or explicitly defined population units. At site $i$, *resets* occur at rate $u_i>0$; a reset removes every trait there, as when an individual dies and is replaced by a newborn with nothing. *Innovations* arrive at rate $\mu_i\ge0$; every innovation creates a new, distinct trait. For each directed donor-to-learner pair $j\to i$ ($j\ne i$) fix an *effective rate* $a_{ji}\ge0$. With *completeness* $0<b\le1$, potential encounters $j\to i$ occur at rate $a_{ji}/b$. At each encounter, every trait held by the donor is copied independently with probability $b$, and the learner keeps the traits it already has. All encounter, reset and innovation processes are independent. Copying decisions are independent across traits and encounters, given the common history. Traits are neutral: they have no prerequisites, capacity limits or effects on any rate.

In the *homogeneous* model, $u_i=u$, $\mu_i=\mu$ and $a_{ji}=\beta/N$ for $i\neq j$. This includes a silent self-choice in the per-learner encounter rate $s=\beta/b$, as in Wakano et al. (2026), so $\beta/N$ is the effective rate of one directed pair.

For a fixed nonempty set $A$ of sites, $C_{A,b}$ is the number of distinct traits held by at least one member of $A$ at a stationary observation time. More generally, for a trait present at the observation time, its *carrier set* is the set of sites holding it. For any family $\mathcal S$ of nonempty sets of sites, $N_{\mathcal S,b}$ counts the traits whose carrier set lies in $\mathcal S$. Examples:

- $C_{A,b}$: carrier set meets $A$;
- the popularity class $P_k$: carrier set has exactly $k$ members;
- $P_{\ge k}$: carrier set has at least $k$ members.

| Symbol | Meaning |
|:--|:--|
| $N$; $u_i$; $\mu_i$ | number of sites; reset (death) rate; innovation rate at site $i$ |
| $a_{ji}$ | effective per-trait learning rate from $j$ to $i$ |
| $b$ | completeness: chance a given trait is copied at an encounter |
| $\beta=sb$ | homogeneous per-learner effective rate (Wakano et al.'s $s$) |
| $\eta=\mu/u$, $\tau=\beta/u$ | dimensionless innovation and learning rates |
| $C_{A,b}$; $P_k$ | repertoire of group $A$; traits held by exactly $k$ sites |
| $R_A$, $M_{A,b}$ | innovation-weighted ancestral volume; its conditional mean (Section 3) |

**The independent-clock limit $b\to0$.** We write $b=0$ for the limit of ever more frequent, ever narrower encounters. In this limit each trait has its own successful-contact clocks at rates $a_{ji}$, but resets are still shared. It does not mean "no learning". Because traits still share deaths, the limiting repertoire is generally not Poisson. Appendix B justifies the limit.

**Correspondence with the discrete versions.** The endpoints match the three versions of Wakano and Aoki (2024, §1.6.1) as follows:

| Here | Wakano–Aoki version | What traits share |
|:--|:--|:--|
| $b=1$ | Version 3 | the successful cultural parents |
| $0<b<1$ | Version 1 | the potential encounters, but not their success |
| $b\to0$ | Version 2 (up to shared deaths) | only deaths |

Overlapping generations and shared deaths make the analogy imperfect, which is why Section 6 treats the discrete model directly.

For finite $N$ and positive reset rates, the backward ancestral process along successful arrows is a finite absorbing Markov chain. Its absorption time therefore has an exponential tail. This yields a stationary construction from innovations in the infinite past, with finite counts and finite moments of every order. A time interval containing no arrows and a reset of every site erases all dependence on the initial state.

## 3. Why a trait count is conditionally Poisson

*Idea.* A trait is present in the sample if it was invented at some site and then copied along a chain of successful encounters that ends in the sample, with no reset breaking the chain. Once the shared history of encounters and deaths is fixed, different traits succeed or fail independently. A count of independent rare successes is Poisson.

Let $D$ be the reset history and $G_b$ the potential-encounter history. Retain each arrow of $G_b$ independently with probability $b$ to obtain a successful-arrow graph $H$; its unconditional edge rates are $a_{ji}$, whatever $b$ is. Trace the sample $A$ backwards through $(H,D)$. When a current ancestral site is the learner at a retained arrow, add the donor; when a reset occurs at an ancestral site, remove it. Write $B_A(t)$ for the ancestral set at backward time $t$, and define the *innovation-weighted ancestral volume*
$$R_A(H,D)=\int_0^\infty\sum_{i\in B_A(t)}\mu_i\,dt,\qquad M_{A,b}=\mathbb E[R_A(H,D)\mid G_b,D].$$
$R_A$ is the expected number of traits the sample would hold if every trait used the single successful graph $H$; $M_{A,b}$ averages over the trait-specific thinnings. With equal innovation rates, $R_A=\mu L_A$, where $L_A$ is the total ancestral length. At $b=1$, $M_{A,1}=R_A$. For $b=0$ we put $M_{A,0}=\mathbb E[R_A\mid D]$. Since $R_A\le N\max_i\mu_i\cdot T_{\rm abs}$ for the absorption time $T_{\rm abs}$, all moments are finite.

More generally, an innovation at site $j$ and backward time $t$ spreads forwards along $H$ and is removed by $D$. This gives a deterministic carrier set at the observation time. Let $R_{\mathcal S}(H,D)=\int_0^\infty\sum_j\mu_j\mathbf 1\{\text{carrier set of an innovation at }(j,-t)\in\mathcal S\}\,dt$ and $M_{\mathcal S,b}=\mathbb E[R_{\mathcal S}\mid G_b,D]$. Since a trait reaches $A$ exactly when its carrier set meets $A$, $R_A$ is the special case $\mathcal S=\{S:S\cap A\neq\varnothing\}$.

**Theorem 1 — conditional Poisson law.** Given $(G_b,D)$, the counts $N_{\mathcal S,b}$ for disjoint families are independent and
$$\boxed{N_{\mathcal S,b}\mid G_b,D\sim\operatorname{Poisson}(M_{\mathcal S,b}).}$$
Consequently $\mathbb E N_{\mathcal S,b}=\mathbb E R_{\mathcal S}$, which does not depend on $b$; $\operatorname{Var}N_{\mathcal S,b}=\mathbb ER_{\mathcal S}+\operatorname{Var}M_{\mathcal S,b}$; and the $k$th factorial moment is $\mathbb E[M_{\mathcal S,b}^k]$. In the $b=0$ limit the conditioning is on $D$ alone.

*Proof.* Given $(G_b,D)$, innovations form a Poisson random measure in space and time. Give each innovation an independent collection of Bernoulli($b$) marks, one per potential arrow. Its carrier set at the observation time is a function of its own marks and of $D$. By the marking theorem, the innovations whose carrier sets fall in disjoint families form independent Poisson processes. The mean number in $\mathcal S$ is the integrated probability of landing in $\mathcal S$, which equals $M_{\mathcal S,b}$ by Tonelli's theorem. Exponential moments of $T_{\rm abs}$ justify the integrals. $\square$

The $b=1$ case of Theorem 1 is the mixed-Poisson ancestry formula of Wakano et al. (2026, equation 15). The new ingredient is the conditional mean over trait-specific thinnings. Theorem 1 characterises the law, but an unevaluated random mean is not a closed-form probability mass function; Section 9 gives finite systems for exact second moments.

## 4. Bundling orders diversity risk

### 4.1 Convex order

For integrable random variables, $X\le_{\rm cx}Y$ means $\mathbb E\phi(X)\le\mathbb E\phi(Y)$ for every convex $\phi$ with defined expectations (Shaked & Shanthikumar, 2007, Chapter 3). Then $X$ and $Y$ have equal means, and $Y$ is more spread out. Equivalently, $\mathbb E(X-q)_+\le\mathbb E(Y-q)_+$ for all real $q$. Convex order does not make $Y$ larger; it makes $Y$ more variable around the same centre.

### 4.2 The main comparison

*Idea.* Observing frequent, less informative encounters is like observing the rarer, fully informative encounters through extra noise. Averaging over extra noise can only reduce variability (Jensen's inequality).

**Theorem 2 — matched-rate comparison.** Fix the innovation, reset and effective edge rates, and a family $\mathcal S$. For $0\le b_1\le b_2\le1$,
$$\boxed{M_{\mathcal S,b_1}\le_{\rm cx}M_{\mathcal S,b_2},\qquad N_{\mathcal S,b_1}\le_{\rm cx}N_{\mathcal S,b_2}.}$$
In particular $C_{A,b_1}\le_{\rm cx}C_{A,b_2}$ for every group $A$, and $P_{k}$ and $P_{\ge k}$ are ordered in the same way for every $k$.

*Proof.* Let $0<b_1<b_2$. Generate $H$ with rates $a_{ji}$. Add independent Poisson arrows with rates $a_{ji}(1/b_2-1)$ to form $G_2$, and then further independent arrows with rates $a_{ji}(1/b_1-1/b_2)$ to form $G_1$. Use the same reset history $D$ throughout. By Poisson superposition and thinning, given $G_2$ the graph $H$ is a $b_2$-thinning, and given $G_1$ it is a $b_1$-thinning. So $(G_1,D)$ and $(G_2,D)$ have the laws of the potential histories at $b_1$ and $b_2$. The extra arrows in $G_1$ are independent of $(H,G_2,D)$, so $H\to(G_2,D)\to(G_1,D)$ is a Markov chain, and
$$M_{\mathcal S,b_1}=\mathbb E[R_{\mathcal S}\mid G_1,D]=\mathbb E\big[\mathbb E[R_{\mathcal S}\mid G_2,D]\mid G_1,D\big]=\mathbb E[M_{\mathcal S,b_2}\mid G_1,D].$$
Conditional Jensen gives the intensity ordering. If $b_1=0$, condition on $D$ instead. The same argument applies at any finite observation time when all sites start empty, restricting innovations and encounters to the interval since that start.

To pass to counts, let $\phi$ be convex on the nonnegative integers and $\psi(x)=\mathbb E\phi(P_x)$ with $P_x\sim{\rm Poisson}(x)$. Then $\psi''(x)=\mathbb E[\phi(P_x+2)-2\phi(P_x+1)+\phi(P_x)]\ge0$, first for bounded $\phi$ and then by monotone approximation. So $\psi$ is convex, and Theorem 1 finishes the proof. $\square$

The Poisson step is the scalar case of Błaszczyszyn and Yogeshwaran (2009, Proposition 3.6). The conditional-expectation step is standard information-garbling reasoning; Dendoncker et al. (2026, Proposition 4.2) give a recent additive-noise analogue in a working paper.

**Corollary 1 — risk measures.** For $b_1\le b_2$ and any family $\mathcal S$ (in particular for $C_A$):

- the expected shortfall $\mathbb E(q-N_{b})_+$ below any target $q$ is nondecreasing in $b$, and so is the expected excess $\mathbb E(N_b-q)_+$;
- the variance is nondecreasing in $b$;
- the probability that the count is zero, $\mathbb E e^{-M_b}$, is nondecreasing in $b$.

**Corollary 2 — the whole spectrum.** For any finite collection of disjoint families (for example $P_1,\dots,P_N$), the covariance matrix increases in the positive-semidefinite order as $b$ increases. Hence $\operatorname{Var}\sum_k w_kP_k$ is nondecreasing in $b$ for every real weight vector $w$, while every $\mathbb EP_k$ is fixed.

*Proof.* By Theorem 1, ${\rm Cov}(N)={\rm diag}(\mathbb E R)+{\rm Cov}(M_b)$, and the diagonal term does not depend on $b$. Along the coupling, ${\rm Cov}(M_{b_2})-{\rm Cov}(M_{b_1})=\mathbb E\,{\rm Cov}(M_{b_2}\mid G_1,D)\succeq0$. $\square$

Corollary 2 is the rigorous form of the remark by Wakano and Aoki (2024, §1.6.2) that package transmission may increase the variance of the cultural spectrum. It is a statement about each linear count and about covariances. It is not a multivariate convex order of the whole spectrum vector.

### 4.3 Exactly when bundling matters

**Proposition 1 — strictness.** Assume all $u_i>0$ and fix a group $A$. The following are equivalent.

(i) There is a site $j$ with $\mu_j>0$ and a directed path $j=v_0\to v_1\to\dots\to v_m\in A$ with $m\ge1$ and $a_{v_{l-1}v_l}>0$ for each step. (The path may start inside $A$ and may return to its starting site.)

(ii) $\operatorname{Var}C_{A,b}$ is strictly increasing on $[0,1]$.

(iii) $\mathbb P(C_{A,b}=0)$ is strictly increasing on $[0,1]$.

If (i) fails, the law of $C_{A,b}$ does not depend on $b$.

The proof (Appendix D) constructs an explicit history of positive probability. On it, a trait invented at $j$ reaches $A$ only through the path's encounters, and whether those encounters succeed remains uncertain given the coarser history. In the homogeneous model, condition (i) holds for every nonempty group when $N\ge2$ and $\beta,\mu>0$. Condition (i) is weaker than requiring a direct link from an inventor into $A$. For example, in the chain $0\to1\to2$ with innovation only at site 0, the repertoire of site 2 has variance strictly increasing in $b$ (Section 11.2).

### 4.4 What convex order does not say

Convex order does not order the probability $\mathbb P(C\le q)$ of falling to or below an arbitrary threshold. This failure is real, not just a gap in the proof. In the two-person example of Section 11 ($u=1$, $\beta=2$, $\mu=1/2$, mean $3/2$), the exact stationary probabilities are:

| $q$ | $b=1/4$ | $b=1/2$ | $b=1$ |
|:--|--:|--:|--:|
| 0 | 0.336751 | 0.346256 | 0.363636 |
| 1 | 0.610286 | 0.614643 | 0.628099 |
| 2 | 0.785021 | 0.784029 | 0.787378 |
| 3 | 0.885961 | 0.882856 | 0.879175 |

Bundling makes very low counts more likely ($q\le1$) and moderately low counts less likely ($q\ge3$), and $q=2$ is not monotone at all. What *is* ordered for every $q$ is the expected shortfall $\mathbb E(q-C)_+$, which weighs how far the count falls below the target.

**Remark 1 — fixed encounter rate instead.** Suppose instead, as in Wakano et al.'s own simulations, that the encounter rate is held fixed and $b$ varies, so effective rates rise with $b$. More generally, suppose effective rates, innovation rates and completeness all weakly increase while reset rates weakly decrease. Then every count $C_A$, and more generally every count for an upward-closed family $\mathcal S$ (such as $P_{\ge k}$), increases in the *increasing convex* order. The proof first applies Theorem 2 at the old effective rates and then a monotone coupling: sharing the old arrows' marks, extra arrows, fewer resets and more innovations can only enlarge each trait's carrier set and add traits, so the conditional mean increases pointwise. This does not extend to exact-frequency counts such as $P_k$, which can fall when carrier sets grow. Thus bundling increases spread under either comparison. Only under matched rates does the mean stay fixed.

### 4.5 The shape of the variance

Write $V(b)=\operatorname{Var}(C_{A,b})$. The Poisson isometry and Mehler's thinning formula (Appendix C) give
$$\boxed{V(b)=V(0)+\sum_{k\ge1}a_kb^k,\qquad a_k\ge0,\qquad \sum_{k\ge1}a_k=V(1)-V(0)<\infty.}$$
Intuitively, $a_k$ is the variance contributed by interactions among $k$ shared encounters. Partial copying damps the $k$th term by $b^k$. So $V$ is increasing and convex, and has nonnegative derivatives of all orders on $[0,1)$. It satisfies the chord bound
$$V(0)\le V(b)\le(1-b)V(0)+bV(1).$$
Under condition (i) of Proposition 1, $a_1>0$. The baseline $V(0)=\mathbb ER_A+\operatorname{Var}(\mathbb E[R_A\mid D])$ is not Poisson: shared deaths alone already create excess variance. When high-order terms dominate, the increase concentrates near $b=1$, which is what Section 7 finds in large supercritical populations.

## 5. Bundled forgetting

Deaths are an extreme form of bundled loss: a reset removes everything at once. Suppose instead that loss events occur at site $i$ at rate $u_i/d$, and each removes each trait present independently with probability $d\in(0,1]$. The effective per-trait loss rate is then $u_i$ for every $d$, and $d=1$ is the original model. Write $C_{A,b,d}$ and $N_{\mathcal S,b,d}$ for the counts. We use $d=0$ for the independent-loss-clock limit, as with $b=0$.

**Theorem 3 — bundled learning and bundled loss.** Let $0\le b_1\le b_2\le1$ and $0\le d_1\le d_2\le1$. Then
$$N_{\mathcal S,b_1,d_1}\le_{\rm cx}N_{\mathcal S,b_2,d_2}$$
for every family $\mathcal S$, with equal means. At $b=d=0$ every such count is exactly Poisson.

*Proof.* Let $D'$ be the effective per-trait loss history (rates $u_i$), and let $D_2$ and $D_1$ add independent loss events at rates $u_i(1/d_2-1)$ and $u_i(1/d_1-1/d_2)$. Build $G_2$ and $G_1$ from $H$ as in Theorem 2. Given $(G_2,D_2)$, the pair $(H,D')$ is obtained by independently thinning arrows with probability $b_2$ and loss events with probability $d_2$. Theorem 1 holds with marks on both arrows and loss events: each trait's carrier set is a function of its own marks. The intensity is the conditional mean of $R_{\mathcal S}(H,D')$, where the single-graph functional now uses the effective loss events as resets. The chain $(H,D')\to(G_2,D_2)\to(G_1,D_1)$ is Markov, and the proof of Theorem 2 applies unchanged. When $b=d=0$, the conditioning is on nothing, so the intensity is the constant $\mathbb ER_{\mathcal S}$. The $d\to0$ limit is justified exactly as the $b\to0$ limit in Lemma B1, with loss events in place of encounters. $\square$

In the two-person example the forgetting channel is the stronger one. With $u=1$, $\beta=2$ and $\mu=1/2$, the variance of the total repertoire (mean $3/2$) is:

| rows $b$, columns $d$ | $d=0$ | $d=1/4$ | $d=1/2$ | $d=1$ |
|:--|--:|--:|--:|--:|
| $0$ | 1.5000 | 1.6611 | 1.8775 | 2.6500 |
| $1/4$ | 1.5082 | 1.6757 | 1.9038 | 2.7500 |
| $1/2$ | 1.5172 | 1.6922 | 1.9342 | 2.8750 |
| $1$ | 1.5385 | 1.7320 | 2.0111 | 3.2500 |

With fully independent loss ($d=0$), copying completeness has almost no effect. Shared deaths are what make shared learning matter. Appendix A gives the closed form in $(b,d)$.

**Remark 2 — gene transfer.** The same two-channel argument applies whenever gains are Poisson along lineages and each gene's fate is an independent mark on shared transfer and loss events. It therefore indicates the direction of the second-moment difference noted by Baumdicker and Pfaffelhuber (2014, Remark 2.4) for multi-gene transfer and loss at fixed per-gene rates: the more genes an event moves or removes, the larger the variance. We do not formalise their coalescent model here.

## 6. The discrete-generation versions of Wakano and Aoki

Wakano and Aoki (2024, §§1.3 and 1.6) pose the distributional question for the discrete-generation model of Kobayashi, Wakano and Ohtsuki (2018). We state that model as they describe it. The 2018 full text could not be obtained, so we rely on the 2024 description; the innovation law is the one detail we must assume.

*Model.* A population of $N$ individuals has non-overlapping generations. Each newborn is assigned $K$ potential models, drawn uniformly without replacement from the previous generation. Each link from a model to the newborn is successful independently with probability $\beta\in(0,1)$. A newborn's repertoire consists of every trait carried by at least one successful model, together with its own new traits. **We assume the number of new traits per newborn is Poisson($\mu$)**, independently of everything else; the sources specify only its mean. The three versions are:

- *Version 1* (the original model): potential models are shared by all traits, and each trait has its own success indicators;
- *Version 2*: each trait has its own potential models;
- *Version 3*: all traits share the successful models.

A single trait has the same law in all three versions. For a sample $A$ of current newborns, let $C^{(v)}_A$ be the number of distinct traits under version $v$, either at stationarity or at a finite horizon starting from an empty population. Let $L_A$ be the total ancestral volume along successful links, $G$ the potential-model assignment and $m=\mu\mathbb EL_A$.

**Theorem 4 — ordering of the versions.** In both regimes, $C^{(2)}_A$ is exactly Poisson($m$), $C^{(1)}_A$ is mixed Poisson with intensity $\mu\,\mathbb E[L_A\mid G]$, and $C^{(3)}_A$ is mixed Poisson with intensity $\mu L_A$. Hence
$$\boxed{C^{(2)}_A\le_{\rm cx}C^{(1)}_A\le_{\rm cx}C^{(3)}_A,}$$
with common mean $m$. The same holds for every count of traits defined by their age and carrier set within the sample, including every entry of the age-frequency spectrum. The corresponding covariance matrices are ordered in the positive-semidefinite order.

*Proof.* Condition on the information that traits share: nothing in Version 2, $\sigma(G)$ in Version 1, and $\sigma(G,H)$ in Version 3, where $H$ is the set of successful links. These $\sigma$-fields are nested. Given them, traits are independent marks of the Poisson innovations, so Theorem 1's argument gives the stated mixed Poisson laws. The intensities form a reverse martingale along the nested $\sigma$-fields, and conditional Jensen and the Poisson lifting of Theorem 2 give the ordering. Version 2's intensity is a constant because, with deterministic generational turnover, no other randomness is shared. $\square$

Appendix E characterises strictness. Versions 1 and 3 differ whenever the horizon is at least one generation. Versions 2 and 1 differ whenever $K<N$, except for a single sampled individual with $K=1$, and except at very short horizons. Two further remarks complete the picture.

- **Interpolation.** A version in which a random share of links are "bundled" (they teach everything or nothing) interpolates monotonically between Versions 1 and 3.
- **Role of the innovation law.** The Poisson assumption matters only for Version 2. Without it, Version 1 $\le_{\rm cx}$ Version 3 still holds, by association (Esary et al., 1967; Christofides & Vaggelatou, 2004). But Version 2's position can change: with at most one innovation per newborn there is an exact counterexample in which Version 2 is the most variable (Appendix E).

The hardest case, Version 1, is thus bracketed by the two versions Wakano and Aoki call potentially solvable. For their illustrative setting ($N=150$, $K=2$, $\beta=1/2$, $\mu=1$, sample of 10), the mean is 327.0. The standard deviation is 18.1 under Version 2 and 315.6 under Version 3.

## 7. How large is the effect?

The theorems give a direction, not a size. We computed $V_N(b)=\operatorname{Var}C_{V,b}$ for the whole population of the homogeneous model ($u=1$, $\mu=0.1$) from the exact pair equations of Section 9 (Figure 2). The values are numerical solutions of exact linear systems with very small residuals; Section 11.2 describes the checks.

![Variance at partial copying relative to whole-repertoire copying, $V_N(b)/V_N(1)$, for four learning rates and population sizes $N=10$ to $150$. Each curve uses the same expected repertoire at every $b$. Below the critical rate ($\beta/u<1$) the ratio stays above about 0.75. Above it, the ratio collapses as $N$ grows and the change concentrates near $b=1$. Exact solutions of the pair-replica equations; $\beta/u=2$, $N=150$ omitted because it could not be solved reliably.](figures/variance_scaling.png){width=100%}

| $\beta/u$ | $N$ | $\mathbb EC$ | $V(0)/V(1)$ | $V(1/2)/V(1)$ | $V(0.9)/V(1)$ | $V(0.99)/V(1)$ |
|:--|--:|--:|--:|--:|--:|--:|
| 0.7 | 50 | 8.27 | 0.783 | 0.843 | 0.950 | 0.994 |
| 0.7 | 150 | 25.4 | 0.749 | 0.814 | 0.937 | 0.993 |
| 1.0 | 150 | 47.6 | 0.232 | 0.320 | 0.616 | 0.931 |
| 1.3 | 50 | 39.1 | 0.117 | 0.195 | 0.482 | 0.892 |
| 1.3 | 100 | 251.2 | 0.039 | 0.080 | 0.271 | 0.773 |
| 1.3 | 150 | 1364 | 0.015 | 0.035 | 0.154 | 0.659 |
| 2.0 | 100 | $6.27\times10^{8}$ | $2\times10^{-6}$ | $4\times10^{-5}$ | 0.015 | 0.466 |

Three patterns stand out. They are numerical findings, not theorems.

1. **Regime dependence.** Below the critical learning rate, the matched-rate comparison changes variance by at most about a quarter at every size computed. Above it, whole-repertoire copying is far more variable, and the ratio falls rapidly with $N$.
2. **A boundary layer at $b=1$.** In supercritical populations, most of the increase occurs at high completeness. At $\beta/u=1.3$ and $N=150$, 86% of the total increase $V(1)-V(0)$ occurs between $b=0.9$ and $b=1$. Slightly incomplete copying from shared encounters already decouples the long shared genealogies.
3. **The source simulations.** Wakano et al. (2026, Figure 16) simulated $N=100$, $\beta=1.3$, $b=0.99$ and $b=13/14$, and observed that the trait count no longer mirrors the genealogy. Exactly, $V(0.99)/V(1)=0.773$ and $V(13/14)/V(1)=0.336$: copying 93% of traits per encounter already removes two-thirds of the whole-repertoire variance.

## 8. Complete-loss episodes

The probability of an empty repertoire is a stationary quantity. For the *whole population*, however, the empty state supplies extra structure. When nobody holds any trait, copying and resets change nothing, and the next event that matters is an innovation. Let $\alpha=\sum_i\mu_i>0$ and $\pi_{0,b}=\mathbb P(C_{V,b}=0)$. The empty interval, of exponential length with mean $1/\alpha$, followed by the next nonempty episode $B_b$, forms a regenerative cycle.

**Corollary 3 — frequency and persistence.** Complete losses occur at long-run rate $r_b=\alpha\pi_{0,b}$. The mean nonempty episode is $\mathbb EB_b=(1-\pi_{0,b})/(\alpha\pi_{0,b})$, and the mean cycle length is $1/(\alpha\pi_{0,b})$. As $b$ increases, losses become more frequent and nonempty episodes become shorter on average, strictly so under Proposition 1's condition. Moreover, $\mathbb Ee^{-zB_b}$ is nondecreasing in $b$ for each $z>0$.

*Proof.* Stationary flow into the empty state equals flow out, $\alpha\pi_{0,b}$. The renewal-reward identity gives $\pi_{0,b}=(1/\alpha)/(1/\alpha+\mathbb EB_b)$. For $b>0$, a block with no innovations or encounters and a reset of every site clears any configuration with fixed positive probability, so returns to the empty state have finite mean. At $b=0$, the chain of carrier-set counts is nonexplosive, because only innovations increase the number of traits and they arrive at constant rate. The empty state is reachable from every configuration. The stationary construction gives an invariant law with positive mass at the empty state, so that state is positive recurrent (Norris, 1997, Theorem 3.5.3). For the transform, start empty at time zero and let $p_b(t)=\mathbb P(C_{V,b}(t)=0)$. Theorem 2 applied at finite horizons gives $p_{b_1}(t)\le p_{b_2}(t)$. Alternating renewal gives $\int_0^\infty e^{-zt}p_b(t)\,dt=[z+\alpha-\alpha\,\mathbb Ee^{-zB_b}]^{-1}$, which is increasing in $\mathbb Ee^{-zB_b}$. $\square$

The renewal identities are standard. Daley et al. (2026, §3.2) record the same cycle formulas for another birth–death setting. Wakano et al. (2026, §§3.3 and 3.5) and Kobayashi et al. (2021) analyse extinction of single traits and genealogies, which are different events: after a single trait disappears, other traits and new innovations can remain.

**How rare is complete loss?** Learning only adds traits, so $\pi_{0,b}\le\prod_i u_i/(u_i+\mu_i)$, the probability with no learning. In large populations this event is extremely rare. In the homogeneous model with $\beta=1.3$, $\mu=0.1$, $u=1$ and $b=1$, the exact values are:

| $N$ | 2 | 10 | 20 | 50 | 100 |
|:--|--:|--:|--:|--:|--:|
| $\pi_{0,1}$ | 0.784 | 0.161 | 0.0181 | $2.3\times10^{-5}$ | $3.2\times10^{-10}$ |

Corollary 3 therefore matters mainly for small groups or for specialised trait domains with few traits. The sawtooth collapses in Wakano et al.'s Figure 16 are falls to low but nonzero levels. They are covered by the expected-shortfall ordering of Corollary 1, not by Corollary 3. The episode results concern literal absence of the modelled trait domain. They do not cover partial losses, arbitrary initial repertoires or recovery after a shock.

## 9. Exact second moments

Run two independent trait-specific thinnings of the same potential history. On an edge of effective rate $a_{ji}$, the two ancestral replicas share successful arrows at rate $a_{ji}b$ and have separate ones at rate $a_{ji}(1-b)$ each; resets are common. For an ancestral set $B$, let $g(B)=\sum_{i\in B}\mu_i$ and let $m(B)$ be its expected integrated reward, so $-Qm=g$ for the single-set generator $Q$. Let $F_b(B_1,B_2)$ be the expected product of the two integrated rewards. First-step analysis gives
$$\boxed{-Q_b^{(2)}F_b(B_1,B_2)=g(B_1)m(B_2)+g(B_2)m(B_1),}$$
with $F_b=0$ if either set is empty. Conditional independence of the thinnings gives
$$\operatorname{Var}(C_{A,b})=m(A)+F_b(A,A)-m(A)^2.$$
This is an exact finite linear system, not a moment closure. Paired ancestral rewards have a close precedent in Baumdicker and Pfaffelhuber (2014, equations 6.6–6.11).

**Homogeneous lumping.** Let $(i,j,k)$ count sites in the first ancestral set only, the second only, and both. Put $n_0=N-i-j-k$ and $a=\beta/N$. The transitions are:

| New state | Rate | Event |
|:--|--:|:--|
| $(i+1,j,k)$ | $an_0[i+(1-b)k]$ | first replica alone gains an outside donor |
| $(i,j+1,k)$ | $an_0[j+(1-b)k]$ | second replica alone gains an outside donor |
| $(i,j,k+1)$ | $abn_0k$ | both gain the same outside donor |
| $(i-1,j,k+1)$ | $ai(j+k)$ | a first-only site joins the second set |
| $(i,j-1,k+1)$ | $aj(i+k)$ | a second-only site joins the first set |
| $(i-1,j,k)$, $(i,j-1,k)$, $(i,j,k-1)$ | $ui$, $uj$, $uk$ | a reset |

A sample of size $n$ starts at $(0,0,n)$. The system has $\binom{N+3}{3}-2N-1$ states. General directed networks need labelled pairs of sets.

**A subtraction-free form.** For equal innovation rates, write $\ell_n=\mathbb EL_n$ and $\Delta_n=\ell_n-\ell_{n-1}$. These satisfy $\Delta_N=1/u$ and $\Delta_n=1/u+\beta(N-n)\Delta_{n+1}/(Nu)$. Computing $F_b-\ell^2$ directly cancels catastrophically in supercritical regimes. Instead, solve for the replica covariance $c=F-\ell\otimes\ell$:
$$-Q_b^{(2)}c=\Gamma_b,\qquad \Gamma_b(i,j,k)=uk\,\Delta_{n_1}\Delta_{n_2}+\tfrac{\beta}{N}b\,n_0k\,\Delta_{n_1+1}\Delta_{n_2+1},$$
with $n_1=i+k$, $n_2=j+k$ and $\Delta_{N+1}=0$. Then $\operatorname{Var}C_{n,b}=\mu\ell_n+\mu^2c(0,0,n)$. Here $\Gamma_b\ge0$ is the rate of joint jumps of the two replicas' expected rewards: shared resets plus shared successful arrows. The excess variance is therefore $\mu^2$ times the expected integrated rate of such shared events, and at $b=0$ only shared resets contribute.

## 10. What variance can identify

Let $P_k$ count the traits held by exactly $k$ of all $N$ individuals in the homogeneous model. One trait's carrier number is a birth–death chain with birth rate $\beta k(N-k)/N$ and death rate $uk$, whatever $b$ is. Balance gives the expected popularity spectrum (Nakamura et al., 2020; Wakano et al., 2026)
$$\mathbb EP_k=\frac{N\eta}{k}\Big(\frac{\tau}{N}\Big)^{k-1}(N-1)^{\underline{k-1}},\qquad \eta=\frac{\mathbb EP_1}{N},\quad\tau=\frac{2N}{N-1}\frac{\mathbb EP_2}{\mathbb EP_1}.$$

**Corollary 4 — structural identification.** For known $N\ge2$ and positive parameters, the population quantities $(\mathbb EP_1,\mathbb EP_2,\operatorname{Var}C_{V,b})$ determine $(\eta,\tau,b)$.

*Proof.* The spectrum means determine $\eta$ and $\tau$. At fixed $\eta$, $\tau$ and $N$, the variance is continuous in $b$, by the power series or the finite resolvent. By Proposition 1 it is also strictly increasing. So a compatible variance determines $b$, and the relative encounter rate is $s/u=\tau/b$. $\square$

**What this does and does not license.**

- A good fit to the *expected* spectrum says nothing about bundling, because the same means occur for every $b$ (Theorem 1).
- "Variance" here means variation of the total repertoire *across population histories*, not variation among individuals within one population. Estimating it requires comparable, independently evolving populations, or a long stationary time series with its dependence accounted for. A large trait catalogue is not a substitute for replication.
- Identification is sharp near $b=1$ and weak near $b=0$. At $\beta=1.3$, $N=100$, the relative sensitivity $V'(b)/V(b)$ is 1.25 at $b=0$, 1.79 at $b=1/2$ and 29.6 at $b=1$.
- The result assumes whole-site resets ($d=1$). If forgetting is also partly bundled, the same variance is compatible with different pairs $(b,d)$ (Section 5), so unknown forgetting confounds the identification of $b$.
- Absolute time is not identified, since multiplying $u$, $\mu$ and $s$ by a constant leaves the stationary law unchanged. A variance outside the model's feasible range signals that the model is wrong, not that $b$ takes some value.
- Differences between populations in environment, innovation or measurement also inflate variance. Whether statistics from a single population (for example, co-carriage of pairs of traits) identify $b$ is left open.

## 11. A two-person example and independent checks

### 11.1 Two individuals

For $N=2$, let $\rho=\beta/(2u)$, the effective rate of each directed edge relative to turnover, and $\eta=\mu/u$. Then $\mathbb EC_{V,b}=\eta(2+\rho)$ and
$$\boxed{\operatorname{Var}C_{V,b}=\eta(2+\rho)+\eta^2\,\frac{\rho^3+(8-2b)\rho^2+10\rho+4}{2+(3-2b)\rho},\qquad \frac{\partial\operatorname{Var}C_{V,b}}{\partial b}=\frac{2\eta^2\rho(\rho+1)(\rho+2)^2}{[2+(3-2b)\rho]^2}>0.}$$
With $v=(\operatorname{Var}C-\mathbb EC)/\eta^2$ and $K_\rho=\rho^3+8\rho^2+10\rho+4$, the inverse is $b=[v(2+3\rho)-K_\rho]/[2\rho(v-\rho)]$. Appendix A derives these formulas.

For $u=1$, $\beta=2$ and $\mu=1/2$, the mean is $3/2$ at every $b$, and the exact variances are

| $b$ | $0$ (limit) | $1/4$ | $1/2$ | $1$ |
|:--|--:|--:|--:|--:|
| $\operatorname{Var}C$ | $53/20$ | $11/4$ | $23/8$ | $13/4$ |
| $\operatorname{Var}P_1$ | $29/20$ | $19/12$ | $7/4$ | $9/4$ |
| $\operatorname{Var}P_2$ | $9/10$ | $1$ | $9/8$ | $3/2$ |

with $\mathbb EP_1=1$ and $\mathbb EP_2=1/2$ throughout. At $b=1$ the empty-state probability is exactly $4/11$.

### 11.2 Checks

The computations deliberately use different representations, and every check script fails loudly (exit status 1) if any comparison fails, also under `python -O`.

- **Forward versus backward.** A forward moment closure on counts of traits with each carrier set uses no ancestral reasoning. In exact rational arithmetic, it agrees with the backward pair system for homogeneous $N=3,4$, a heterogeneous directed three-site network and the chain $0\to1\to2$. In every case means are constant in $b$, and the variances of $C_A$, of every $P_k$ and of every $P_{\ge k}$ increase strictly. In the chain, with innovation only at site 0, the repertoire of site 2 has variance $0.7192$, $0.8089$, $0.9296$, $1.0935$ and $1.3195$ at $b=0,1/4,1/2,3/4,1$.
- **Bundled forgetting.** The same closure extended to partial loss confirms Theorem 3 on a grid of $(b,d)$ for $N=2,3$, including the Poisson law at $b=d=0$.
- **Boundaries.** With no learning, each site holds a geometric number of traits and the sample count is negative binomial. The case $N=1$ is geometric at every $b$. At $b=1$ the pair system reproduces the single-genealogy recursion.
- **Forward simulation.** An independent simulation tracks actual trait identities through 100,000 complete regenerative cycles at each of $b=0.1,0.25,0.5,1$ for $N=2$. The largest standardised deviations from the exact mean and variance were $0.71$ and $1.31$ standard errors. The simulated empty probabilities, $0.3310$, $0.3365$, $0.3479$ and $0.3634$, increase with $b$.
- **Large populations.** Section 7's values come from the subtraction-free system of Section 9. It is solved by preconditioned iterative solves inside iterative refinement, with residuals evaluated in extended precision from exact rates. Final relative residuals are below $10^{-26}$. The values agree with direct sparse factorisation for $N\le50$ (relative difference at most $9\times10^{-11}$) and with the single-genealogy recursion at $b=1$. For example, at $N=50$, $\beta=1.3$, the mean is $39.14076494$ and the variance is $177.97071$ at $b=1/2$ and $914.95061$ at $b=1$. These are floating-point computations with residual checks, not interval certificates.
- **Discrete versions.** Exact rational pair systems for $N=3$ and $K=1,2,3$ agree with brute-force enumeration and with trait-identity simulations. They confirm Theorem 4 and the strictness pattern of Appendix E.

These are computational checks of the mathematics, not observations of human populations.

## 12. Discussion

**What the result says.** Each trait's own chance of surviving is the same under every arrangement; bundling changes only whether traits are gained and lost together. When every trait is passed on at the same average rate, organising learning into bigger bundles leaves the expected amount and expected distribution of culture unchanged. It makes the fortunes of different traits more strongly linked. Bundled encounters are common shocks: they carry many traits together, and the same traits are later lost together when their few carriers die. Partial encounters spread a trait's fate over many independent chances. The consequence is wider variation across possible histories. Both prosperous and impoverished histories become more likely, and so does complete loss. With bundled forgetting the effect compounds, and in well-connected, growing-repertoire (supercritical) populations it is large.

**Analogues.** The same structure appears elsewhere: in linkage versus recombination in population genetics (Section 1.5), in the variance of pangenome size under multi-gene transfer (Remark 2), and in epidemics, where concentrating transmission in fewer events increases variability at a fixed mean reproduction number (Lloyd-Smith et al., 2005). The present paper supplies an exact ordering for cultural accumulation, not only a second-moment observation.

**What it does not say.** A trait count is neither a measure of cultural worth nor a complete measure of complexity. The theorems assume neutral traits that are independent given the shared history. Prerequisites, mutually exclusive variants, limited memory, prestige that depends on repertoire, or contact rates that respond to culture can break the Poisson marking argument. With prerequisites, bundling might even protect a dependent cluster. The matched-rate comparison is a normalisation, not a description of a real trade-off of time or attention; Remark 1 shows that bundling also increases spread when encounters are held fixed. The result does not explain any specific archaeological loss, and it does not order the probability of falling below an arbitrary threshold (Section 4.4).

**Open questions.** Three questions remain open: a large-population theory of the boundary layer in Section 7; single-population statistics that identify bundling; and non-Poisson innovation in the discrete model.

**Data and code.** All code, machine-readable results and figures accompany this paper. No empirical data are used. Version 0.1.0-candidate is archived at <https://doi.org/10.5281/zenodo.23933211> (all versions: <https://doi.org/10.5281/zenodo.23933210>), with the code at <https://github.com/ipitchford/bundled-cultural-transmission-risk>. Before release, the exact checks of Section 11.2 were replayed from a clean copy of the package in normal and optimised Python modes; the large-population computations of Section 7 were not rerun. Assurance is limited to internal AI-assisted review: the paper has not had human peer review, and the proofs are written arguments, not machine-verified.

## Appendix A. Two-person formulas

Rescale time so that $u=1$. Let $X,Y,Z$ count the traits held only by person 1, only by person 2, and by both.

- An encounter from person 1 to person 2, at rate $\rho/b$, moves a Binomial$(X,b)$ number of traits from $X$ to $Z$. The reverse encounter acts on $Y$ in the same way.
- A death of person 1 sends $(X,Y,Z)$ to $(0,Y+Z,0)$, and a death of person 2 acts symmetrically.
- Innovations add one to $X$ or $Y$ at rate $\eta$ each.

The stationary means are $\mathbb EX=\mathbb EY=\eta$ and $\mathbb EZ=\rho\eta$. Put $f=\mathbb E[X(X-1)]$, $c=\mathbb E[XY]$, $w=\mathbb E[XZ]$ and $g=\mathbb E[Z(Z-1)]$. The generator gives
$$c=\frac{\eta^2}{1+\rho},\quad [1+2\rho(1-b)]f-2(1+\rho)w=2\eta^2,\quad (2+\rho)w=\rho\eta^2+\rho c+\rho(1-b)f,\quad g=\rho(2w+bf).$$
Solving,
$$f=\frac{2\eta^2(1+\rho)(2+\rho)}{2+(3-2b)\rho},\qquad w=\frac{\rho\eta^2}{1+\rho}+\frac{\rho(1-b)f}{2+\rho}.$$
Since $\mathbb E[C(C-1)]=2f+g+2c+4w$, Section 11's formula follows. The inverse has a positive denominator because $v-\rho=(\rho+1)(\rho+2)^2/[2+(3-2b)\rho]>0$.

With bundled forgetting at completeness $d$, the same closure (solved symbolically) gives
$$\operatorname{Var}C_V=\eta(2+\rho)+\eta^2\frac{d\,\mathcal A(d,\rho)-2b\rho(\rho+1)[d(\rho+1)-1]}{(2-d)\mathcal B(d,\rho)-2b\rho(\rho+2-d)},$$
where
$$\mathcal A=(\rho^2+2\rho+2)d^2+(\rho^4+\rho^3-7\rho^2-16\rho-10)d+4(\rho+1)(2\rho^2+4\rho+3),$$
$$\mathcal B=d^2+(\rho^2-3\rho-5)d+2(\rho+1)(\rho+3).$$
It reduces to the formula above at $d=1$ and agrees with the exact grid of Section 5.

## Appendix B. Boundary cases and the $b\to0$ limit

With no learning, each of $n$ sampled sites holds the innovations made since its last reset: a geometric count with mean $\mu/u$. The sample count is negative binomial, with $\mathbb EC=n\mu/u$ and $\operatorname{Var}C=n\{\mu/u+(\mu/u)^2\}$. Shared reset uncertainty survives even when copying is completely separated across traits. For $N=1$ this law holds at every $b$, and the strict comparison correctly disappears.

**Lemma B1 — the independent-clock limit.** As $b\downarrow0$, $C_{A,b}$ converges in law to the independent-clock model, and $V(b)\to V(0)$.

*Proof.* Fix a window $[-T,0]$ and condition on the resets and on the $K$ innovations made in it. Split the potential encounters on each edge according to the set of these traits they copy successfully. Encounters that copy two or more of them form a Poisson process whose expected number in the window is at most $T\binom K2 b\sum_{j\ne i}a_{ji}\to0$. Encounters that copy exactly one given trait form independent Poisson processes with rates $a_{ji}(1-b)^{K-1}$. Couple each with a rate-$a_{ji}$ process by adding independent points of rate $a_{ji}[1-(1-b)^{K-1}]$; the probability that any added point falls in the window also tends to zero. Off these two events, the traits' fates coincide with those of the independent-clock model. The expected number of present traits older than $T$ does not depend on $b$ and tends to zero as $T\to\infty$. This proves convergence in law. Convergence of the variance follows from the power series of Section 4.5, which is continuous on $[0,1]$. Theorem 2 at $b_1=0$ uses only the definition $M_{A,0}=\mathbb E[R_A\mid D]$, so it does not depend on this limit. $\square$

## Appendix C. The variance series

Let $\Lambda$ be the effective-arrow intensity measure on directed edges and backward time. Conditional on almost every reset history $D$, $R_D(H)=R_A(H,D)$ is square integrable. With the add-one operator $\Delta_xR_D(H)=R_D(H+\delta_x)-R_D(H)$, put
$$T_{k,D}(x_1,\dots,x_k)=\mathbb E_H[\Delta_{x_1}\cdots\Delta_{x_k}R_D(H)],\qquad a_k=\frac1{k!}\mathbb E_D\|T_{k,D}\|^2_{L^2(\Lambda^k)}.$$
Two independent $b$-thinnings $H_1,H_2$ of $G_b$ share arrows of intensity $b\Lambda$, and each has private intensity $(1-b)\Lambda$. Given $H_1$, $H_2$ is therefore a $b$-thinning of $H_1$ plus independent Poisson immigration of intensity $(1-b)\Lambda$. This is the Mehler operator $P_b$ of Last, Peccati and Schulte (2016, equation 3.1), and their Lemma 3.1 gives $T_k(P_bR_D)=b^kT_k(R_D)$. The isometry of Last and Penrose (2011, Theorem 1.1) then yields
$$\mathbb E[M_{A,b}^2\mid D]=\mathbb E_H[R_DP_bR_D]=(\mathbb E_HR_D)^2+\sum_{k\ge1}\frac{b^k}{k!}\|T_{k,D}\|^2.$$
Integrating over $D$ and adding the $b$-invariant Poisson term $\mathbb ER_A$ proves Section 4.5. Summability follows from the $b=1$ case.

Adding an arrow cannot shrink the ancestral set, so $\Delta_xR_D\ge0$. Under condition (i) of Proposition 1, the witness of Appendix D gives $\mathbb E_H\Delta_xR_D>0$ on a set of positive $\mathbb P_D\otimes\Lambda$ measure. Hence $a_1=\mathbb E_D\int(\mathbb E_H\Delta_xR_D)^2\,\Lambda(dx)>0$. The identities are established; their use here is model-specific.

## Appendix D. Proof of Proposition 1

*(i) fails $\Rightarrow$ no dependence on $b$.* Call a site *fed* if some site with positive innovation rate has a positive-rate path of length at least one to it. If (i) fails, no member of $A$ is fed. Every site that enters the backward ancestral set through an arrow then has $\mu=0$: otherwise that site would have a path of length at least one into $A$. A member of $A$ that has been reset and later re-enters through an arrow also contributes nothing, for the same reason. Hence $R_A=\sum_{i\in A}\mu_i\,\sigma_i$, where $\sigma_i$ is the backward time to the first reset at $i$. This depends on $D$ alone, so $M_{A,b}=R_A$ for every $b$, and the law of $C_{A,b}$ does not depend on $b$.

*(i) $\Rightarrow$ strictness.* Choose a shortest path $v_0\to\dots\to v_m$ satisfying (i), with $m\ge1$, $\mu_{v_0}>0$ and $v_m\in A$. By minimality, the intermediate sites $v_1,\dots,v_{m-1}$ are not in $A$ and have zero innovation rate. Let $0<b_1<b_2\le1$, and build $H\subseteq G_2\subseteq G_1$ as in Theorem 2. Fix backward times $T_0>T>0$ and forward times $-T<s_1<\dots<s_m<0$. Let $E$ be the event, measurable with respect to $(G_1,D)$, on which the following hold:

- in $[-T_0,-T]$ every site is reset and $G_1$ has no arrows (so nothing before $-T_0$ reaches the sample);
- in $(-T,0]$, $G_1$ has exactly the $m$ arrows $v_{l-1}\to v_l$ at times near $s_l$, and no others;
- $v_0$ has no reset in $(-T,s_1)$;
- each $v_l$ has no reset between receiving at $s_l$ and sending at $s_{l+1}$;
- $v_m$ has no reset in $(s_m,0]$;
- if $v_0\in A$, then $v_0$ has a reset in $(s_1,0)$, and in $(s_1,s_m)$ if $v_0=v_m$.

These conditions constrain finitely many independent Poisson processes on bounded intervals, so $\mathbb P(E)>0$.

Let $t_1$ be the actual time of the first path arrow and $r_0\in[-T_0,-T]$ the last reset of $v_0$ before it. On $E$, every contribution to $R_A$ other than the inventions of $v_0$ in $(r_0,t_1)$ is a function of $D$, because those contributions involve no arrows. The inventions of $v_0$ in $(r_0,t_1)$ reach $A$ if and only if all $m$ path arrows are in $H$. Thus on $E$,
$$M_{A,b_2}=c(D)+b_2^m\,\mu_{v_0}(t_1-r_0)\,\mathbf 1\{\text{all path arrows in }G_2\}.$$
Given $(G_1,D)$, $G_2$ keeps each arrow of $G_1$ independently with probability $b_1/b_2\in(0,1)$. The indicator is therefore Bernoulli with parameter $(b_1/b_2)^m\in(0,1)$. So $\operatorname{Var}(M_{A,b_2}\mid G_1,D)>0$ on $E$.

Hence
$$\operatorname{Var}C_{A,b_2}-\operatorname{Var}C_{A,b_1}=\mathbb E\operatorname{Var}(M_{A,b_2}\mid G_1,D)>0,$$
and strict conditional Jensen for $x\mapsto e^{-x}$ gives (iii). For $b_1=0$, condition on $D$ alone. Given $D\in E_D$ (the reset part of $E$), compare two events for $G_2$ of positive conditional probability: no arrows in $[-T_0,0]$, and exactly the prescribed path arrows. Then $M_{A,b_2}$ equals $c(D)$ on the first and exceeds it on the second, so $\operatorname{Var}(M_{A,b_2}\mid D)>0$. Finally, (ii) and (iii) each fail when (i) fails, by the first part. $\square$

## Appendix E. The discrete versions: strictness and innovation law

Let $T$ be the horizon in generations. Let $\mu>0$, $0<\beta<1$ and sample size $n\ge1$.

**Strictness.**

1. $C^{(1)}_A$ and $C^{(3)}_A$ differ in law if and only if $T\ge1$. Given $G$, all of the sample's links fail with probability $(1-\beta)^{nK}$, and all succeed with probability $\beta^{nK}$. These give different ancestral volumes, so $\operatorname{Var}(L_A\mid G)>0$ for every $G$, and $\operatorname{Var}C^{(3)}_A-\operatorname{Var}C^{(1)}_A=\mu^2\,\mathbb E\operatorname{Var}(L_A\mid G)>0$.
2. $C^{(2)}_A$ and $C^{(1)}_A$ differ in law if and only if $K<N$ and either (a) $n\ge2$ and $T\ge1$, or (b) $n=1$, $K\ge2$ and $T\ge2$.
   - In case (a), compare an assignment in which two sampled newborns share all their potential models with one in which they differ in one model. The expected number of distinct ancestors one generation back differs. That number feeds a strictly increasing expected future volume.
   - Case (b) is the same argument one generation further back.
   - Necessity: if $K=N$, then $G$ is deterministic. If $n=1$ and $K=1$, then for every $G$ the ancestry is a single path surviving to depth $a$ with probability $\beta^a$. If $n=1$ and $T=1$, the number of ancestors is Binomial$(K,\beta)$ for every $G$. In each case the Version 1 intensity is constant.

**Non-Poisson innovation.** Let the innovation numbers have any law independent of transmission, with bounded means. Given $G$ and the innovation numbers, each trait's presence indicator is a nondecreasing function of independent link indicators. The vector of these indicators is therefore associated (Esary et al., 1967). Version 1 uses independent copies of the link indicators for each trait, so its vector is the independent duplicate of Version 3's. An associated vector dominates its independent duplicate in the supermodular order (Christofides & Vaggelatou, 2004). Applying this to the convex function of the sum gives Version 1 $\le_{\rm cx}$ Version 3, conditionally and hence unconditionally.

The comparison with Version 2 can fail. Take $N=3$, $K=1$, $\beta=9/10$, $T=3$, a single sampled newborn, and at most one innovation per newborn, occurring with probability $1/2$. Then the variances are $0.970115$ (Version 1), $1.116320$ (Version 3) and $1.303038$ (Version 2), with common mean $1.7195$.

## References

Aguilar, E., & Ghirlanda, S. (2015). Modeling the genealogy of a cultural trait. *Theoretical Population Biology, 101*, 1–8. https://doi.org/10.1016/j.tpb.2014.12.007

Baumdicker, F., & Pfaffelhuber, P. (2014). The infinitely many genes model with horizontal gene transfer. *Electronic Journal of Probability, 19*, 1–27. https://doi.org/10.1214/EJP.v19-2642

Błaszczyszyn, B., & Yogeshwaran, D. (2009). Directionally convex ordering of random measures, shot noise fields, and some applications to wireless communications. *Advances in Applied Probability, 41*(3), 623–646. https://doi.org/10.1239/aap/1253281057

Christofides, T. C., & Vaggelatou, E. (2004). A connection between supermodular ordering and positive/negative association. *Journal of Multivariate Analysis, 88*(1), 138–151. https://doi.org/10.1016/S0047-259X(03)00064-2

Daley, D. J., Nazarathy, Y., & Wang, J. (2026). Second order properties of thinned counts in finite birth–death processes. *Queueing Systems, 110*(1), Article 17. https://doi.org/10.1007/s11134-026-09974-5

Dendoncker, V., Denuit, M., & Robert, C. Y. (2026). *Reverse-time scores and convex-order monotonicity for additive processes with Lévy clocks* (LIDAM Discussion Paper ISBA 2026/23). UCLouvain. https://hdl.handle.net/2078.5/276938

Esary, J. D., Proschan, F., & Walkup, D. W. (1967). Association of random variables, with applications. *The Annals of Mathematical Statistics, 38*(5), 1466–1474. https://doi.org/10.1214/aoms/1177698701

Henrich, J. (2004). Demography and cultural evolution: How adaptive cultural processes can produce maladaptive losses—the Tasmanian case. *American Antiquity, 69*(2), 197–214. https://doi.org/10.2307/4128416

Hudson, R. R. (1983). Properties of a neutral allele model with intragenic recombination. *Theoretical Population Biology, 23*(2), 183–201. https://doi.org/10.1016/0040-5809(83)90013-8

Kobayashi, Y., Kurokawa, S., Ishii, T., & Wakano, J. Y. (2021). Time to extinction of a cultural trait in an overlapping generation model. *Theoretical Population Biology, 137*, 32–45. https://doi.org/10.1016/j.tpb.2021.01.002

Kobayashi, Y., Wakano, J. Y., & Ohtsuki, H. (2018). Genealogies and ages of cultural traits: An application of the theory of duality to the research on cultural evolution. *Theoretical Population Biology, 123*, 18–27. https://doi.org/10.1016/j.tpb.2018.04.007

Kolodny, O., Creanza, N., & Feldman, M. W. (2015). Evolution in leaps: The punctuated accumulation and loss of cultural innovations. *Proceedings of the National Academy of Sciences, 112*(49). https://doi.org/10.1073/pnas.1520492112

Last, G., & Penrose, M. D. (2011). Poisson process Fock space representation, chaos expansion and covariance inequalities. *Probability Theory and Related Fields, 150*(3–4), 663–690. https://doi.org/10.1007/s00440-010-0288-5

Last, G., Peccati, G., & Schulte, M. (2016). Normal approximation on Poisson spaces: Mehler's formula, second order Poincaré inequalities and stabilization. *Probability Theory and Related Fields, 165*(3–4), 667–723. https://doi.org/10.1007/s00440-015-0643-7

Lehmann, L., Aoki, K., & Feldman, M. W. (2011). On the number of independent cultural traits carried by individuals and populations. *Philosophical Transactions of the Royal Society B: Biological Sciences, 366*(1563), 424–435. https://doi.org/10.1098/rstb.2010.0313

Lloyd-Smith, J. O., Schreiber, S. J., Kopp, P. E., & Getz, W. M. (2005). Superspreading and the effect of individual variation on disease emergence. *Nature, 438*(7066), 355–359. https://doi.org/10.1038/nature04153

Nakamura, M., Wakano, J. Y., Aoki, K., & Kobayashi, Y. (2020). The popularity spectrum applied to a cross-cultural question. *Theoretical Population Biology, 133*, 104–116. https://doi.org/10.1016/j.tpb.2019.10.003

Norris, J. R. (1997). *Markov chains*. Cambridge University Press. https://doi.org/10.1017/CBO9780511810633

Pluzhnikov, A., & Donnelly, P. (1996). Optimal sequencing strategies for surveying molecular genetic diversity. *Genetics, 144*(3), 1247–1262. https://doi.org/10.1093/genetics/144.3.1247

Scheie, F., Jansson, F., & Sjöstrand, J. (2026). Sequential accumulation of independent cultural traits: A Gaussian limit for the popularity spectrum. *Theoretical Population Biology, 172*, 107–121. https://doi.org/10.1016/j.tpb.2026.10.003

Shaked, M., & Shanthikumar, J. G. (2007). *Stochastic orders*. Springer. https://doi.org/10.1007/978-0-387-34675-5

Strimling, P., Sjöstrand, J., Enquist, M., & Eriksson, K. (2009). Accumulation of independent cultural traits. *Theoretical Population Biology, 76*(2), 77–83. https://doi.org/10.1016/j.tpb.2009.04.006

Wakano, J. Y., & Aoki, K. (2024). Forward and backward modeling of cultural evolutionary processes. *Proceedings of the National Academy of Sciences, 121*(48), e2322880121. https://doi.org/10.1073/pnas.2322880121

Wakano, J. Y., Ohtsuki, H., Kobayashi, Y., & Baake, E. (2026). Evolving genealogies in cultural evolution, the descendant process, and the number of cultural traits. *Theoretical Population Biology, 168*, 1–18. https://doi.org/10.1016/j.tpb.2025.12.003

Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology, 7*(2), 256–276. https://doi.org/10.1016/0040-5809(75)90020-9

Yeh, D. J., Fogarty, L., & Kandler, A. (2019). Cultural linkage: The influence of package transmission on cultural dynamics. *Proceedings of the Royal Society B: Biological Sciences, 286*(1916), 20191951. https://doi.org/10.1098/rspb.2019.1951
