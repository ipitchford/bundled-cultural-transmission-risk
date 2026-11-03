---
title: "The discrete-generation versions of Wakano and Aoki"
subtitle: "Supplementary section: an exact convex ordering of Versions 1–3"
date: "3 November 2026 · Draft section for lead review"
---

## D.1 The question in its original model

Wakano and Aoki (2024, §1.6.2) ask for the *distribution*, and not only the expected value, of the number of cultural traits in the discrete-generation model of Kobayashi, Wakano and Ohtsuki (2018). They attribute the task to Kobayashi et al. (2018). They describe three ways of applying one ancestral process to many traits (§1.6.1). They suggest that Version 3, with one genealogy, "may be solvable" but fluctuates greatly; that Version 2, with independent genealogies, may also be solvable; and that Version 1, the original model, is probably the hardest because its genealogies are correlated. This section answers the qualitative part of that question exactly, in their model. At any fixed sample and observation generation, the three repertoire counts have the same mean and are ordered in convex order, Version 2 $\leq_{\rm cx}$ Version 1 $\leq_{\rm cx}$ Version 3. Version 2 is exactly Poisson. Versions 1 and 3 are mixed Poisson with explicit random means.

**Source access.** The full text of Kobayashi et al. (2018) could not be obtained. ScienceDirect refused automated access, and no preprint or author copy was found. The model below therefore follows the detailed description in Wakano and Aoki (2024, §§1.3–1.6), the 2018 abstract (as indexed by RePEc), and the summary of the 2018 model in Wakano et al. (2026, §2). The places where these sources leave a detail open are listed in §D.2.

## D.2 The model, as defined in the sources

*Population and generations.* A population of fixed size $N$ has discrete, non-overlapping generations "as in the Wright–Fisher process" (Wakano & Aoki 2024, §1.3.1). Label generations by $t\in\mathbb Z$, with the observation generation at $t=0$, and individuals by $[N]=\{1,\ldots,N\}$.

*Potential models and cultural parents.* Each newborn $(t,i)$ is assigned $1\le K\le N$ *potential behavioural models*, drawn at random without replacement from generation $t-1$; different newborns may receive different sets (§1.3.1). Write $G_t(i)\subseteq[N]$ for this $K$-set. The sets $G_t(i)$ are independent across newborns and uniformly distributed. Each link from a potential model to the newborn is *successful* independently with probability $0<\beta<1$. The successful models are the newborn's *cultural parents*, and their number $P$ is binomial with mean $K\beta$ (§1.3.1, Fig. 1). A cultural parent always transmits the trait in question if it carries it, and transmission is redundant when several parents carry it.

*Innovation.* Each newborn "innovates on average $\mu$ traits", and every innovation is novel (§1.3.2). Wakano et al. (2026, §2) summarise the 2018 model as one in which "every newborn invents a random number of new traits". Neither secondary source states the law of this number. **We assume that innovation numbers are independent Poisson($\mu$) variables, independent of all transmission randomness.** Section D.8 shows that this assumption carries weight.

*Timing.* Equation 1 of Wakano and Aoki (2024) gives $\mathbb E[C_{a,n}]=\mu\,\mathbb E[Y_a]$, with $Y_0=n$, so traits invented by the sampled newborns themselves (age $0$) are counted. A newborn's repertoire therefore consists of the traits learned from its cultural parents together with its own innovations, and the whole repertoire is available to the next generation.

*Versions* (§1.6.1).

- **Version 1** (the 2018 model; §1.3.2 and §1.6.1). All traits share each newborn's potential-model set $G_t(i)$. Cultural parents are chosen independently for different traits: every trait has its own independent Bernoulli($\beta$) success indicator on every link.
- **Version 2.** For each trait, every newborn draws its own potential-model set, and success indicators are independent across traits.
- **Version 3.** All traits share both the potential models and the cultural parents, and "all traits possessed by each cultural parent are transmitted simultaneously". Here $P$ is again binomial with mean $K\beta$.

*What is observed.* A sample of $n$ individuals is taken from the current generation. We fix it as a labelled set $A\subseteq[N]$ with $|A|=n$. Exchangeability makes the choice of labels immaterial. The sources report the number $C_n$ of distinct traits in the sample (eq. 2), the age spectrum $C_{a,n}$ (eq. 1) and the age-frequency spectrum $C_{a,n,i}$, which counts traits of age $a$ carried by exactly $i$ sample members (§1.4.2; 2018 abstract).

*Horizon.* We consider two observation regimes. In the *finite-horizon* regime, generation $-T-1$ carries no traits and innovation occurs in generations $-T,\ldots,0$, with $0\le T<\infty$. In the *stationary* regime, innovation has occurred in every generation $t\le0$ ($T=\infty$).

**Ancestral sets.** A *transmission graph* $h=(h_t(i))$ assigns to each newborn a set $h_t(i)\subseteq[N]$ of cultural parents. For a sample $A$, put $A_0(h)=A$ and
$$A_{t-1}(h)=\bigcup_{i\in A_t(h)}h_t(i)\qquad(t\le0).$$
Induction over generations shows that a trait invented by $(s,j)$ and transmitted along $h$ is carried by some member of $A$ at time $0$ if and only if $j\in A_s(h)$. More precisely, member $i\in A$ carries it if and only if $j\in A_s^{\{i\}}(h)$. Then $Y_a=|A_{-a}(h)|$ is the ancestor count of Wakano and Aoki (§1.3.1). Write
$$L_A(h)=\sum_{s=-T}^{0}|A_s(h)|=\sum_{a=0}^{T}Y_a$$
for the total ancestral volume.

*Generic graphs.* Let $G$ be the potential-model assignment and let $H$ be its $\beta$-thinning: $H_t(i)$ keeps each element of $G_t(i)$ independently with probability $\beta$. In Version 1, trait $\tau$ is transmitted along its own thinning $H^\tau$ of the shared $G$. In Version 2, it uses a thinning $H^\tau$ of its own assignment $G^\tau$. In Version 3, every trait uses $H$.

*Finiteness.* Whatever the starting set, the event that every newborn's links all fail has probability $(1-\beta)^{NK}>0$ in each generation. The ancestor chain $Y$ is therefore absorbed at $0$ after a time with a geometric tail. In particular $L_A(H)\le N\cdot(\text{absorption time}+1)$ has finite moments of every order.

## D.3 The single-trait marginals coincide

**Lemma D0.** In all three versions, the transmission graph used by any one trait has the law of $H$. The fate of a single trait, including its carrier set at every later generation, therefore has the same law in all three versions.

*Proof.* In Version 1, $H^\tau$ is a $\beta$-thinning of a uniform independent assignment $G$. In Version 2, it is a $\beta$-thinning of the uniform independent assignment $G^\tau$. In Version 3, it is $H$. In each case the newborns' successful sets are independent, uniformly random $K$-sets thinned with probability $\beta$. $\square$

Version 3 as written by Wakano and Aoki therefore matches the other two versions; no reparametrisation is needed. The three versions differ only in how much of this common single-trait randomness is shared between traits:
$$\mathcal F_2=\{\varnothing,\Omega\}\subseteq\mathcal F_1=\sigma(G)\subseteq\mathcal F_3=\sigma(G,H).$$
Given $\mathcal F_v$, the traits of version $v$ use independent copies of the conditional law of $H$ given $\mathcal F_v$. This is also why Wakano and Aoki observe that the expected values of Kobayashi et al. hold in every version (§1.6.1).

## D.4 Conditional Poisson representation

Write $m=\mu\,\mathbb E L_A(H)=\mu\sum_{a=0}^T\mathbb E Y_a$, the expected count of eq. 2 of Wakano and Aoki. Write $C^{(v)}_A$ for the number of distinct traits in $A$ at generation $0$ under version $v$, and define the random intensities
$$\Lambda_3=\mu L_A(H),\qquad \Lambda_1=\mathbb E[\mu L_A(H)\mid G],\qquad \Lambda_2=m.$$

**Theorem D1 — conditional repertoire laws.** In both regimes,
$$\boxed{C^{(v)}_A\mid\mathcal F_v\sim\operatorname{Poisson}(\Lambda_v),\qquad v=1,2,3.}$$
In particular, $C_A^{(2)}$ is exactly Poisson with mean $m$, and
$$\mathbb E C^{(v)}_A=m,\qquad \operatorname{Var}C^{(v)}_A=m+\operatorname{Var}\Lambda_v,\qquad \mathbb E\big[(C^{(v)}_A)_k\big]=\mathbb E[\Lambda_v^k].$$
The conditional law holds jointly for the age-frequency classes as well (Corollary D4).

*Proof.* Fix $v$ and condition on $\mathcal F_v$. The innovation numbers $\Pi_{s,j}$ of the newborns $(s,j)$, $-T\le s\le0$, are independent Poisson($\mu$) variables, independent of $\mathcal F_v$. Each trait of $(s,j)$ carries an independent mark $H^\tau$ drawn from the conditional law $\kappa_v$ of $H$ given $\mathcal F_v$. For Version 3, $\kappa_3$ is a point mass, and all traits are marked with $H$. The trait reaches $A$ if and only if $j\in A_s(H^\tau)$, an event of conditional probability $p_{s,j}=\kappa_v\{h:j\in A_s(h)\}$. Poisson thinning makes the number of traits of $(s,j)$ that reach $A$ Poisson with mean $\mu p_{s,j}$, independently over $(s,j)$. A countable sum of independent Poisson variables whose means have a finite sum is Poisson, as the product of their generating functions shows. The total mean is
$$\sum_{s,j}\mu p_{s,j}=\mu\int\sum_{s,j}\mathbf 1\{j\in A_s(h)\}\,\kappa_v(dh)=\mu\int L_A(h)\,\kappa_v(dh)=\Lambda_v,$$
by Tonelli's theorem. It is finite almost surely because $\mathbb E L_A(H)<\infty$. The moment formulas follow from those of mixed Poisson laws and the finite moments of $L_A$. $\square$

In the non-overlapping model, generational turnover is deterministic, and no random event is shared between traits other than those in $\mathcal F_v$. Version 2 is therefore an ordinary Poisson variable: given innovation and turnover structure, nothing else correlates traits' fates. This contrasts with the continuous-time model of the main text, where shared resets keep even the limiting independent-contact model away from Poisson (§2). Version 2 therefore has the closed-form law $\mathbb P(C^{(2)}_A=k)=e^{-m}m^k/k!$, and in particular $\mathbb P(C^{(2)}_A=0)=e^{-m}$.

Version 3 is also explicitly computable. Its mixing variable is $\mu$ times the total length of the single chain $Y$, so
$$\mathbb E z^{C^{(3)}_A}=\mathbb E\exp\{-\mu(1-z)L_A(H)\}.$$
The right-hand side solves a finite linear system on the $N+1$ states of $Y$. This confirms Wakano and Aoki's suggestion that Version 3 "may be solvable". Version 1 has the intensity $\mu\,\mathbb E[L_A(H)\mid G]$, a conditional expectation over a random graph. No finite closed form is claimed for it here.

## D.5 The ordering

For integrable $X,Y$, $X\le_{\rm cx}Y$ means $\mathbb E\phi(X)\le\mathbb E\phi(Y)$ for every convex $\phi$ for which both expectations exist. Equivalently, $\mathbb EX=\mathbb EY$ and $\mathbb E(X-a)_+\le\mathbb E(Y-a)_+$ for every real $a$ (Shaked & Shanthikumar 2007, Theorem 3.A.1).

**Lemma D2 — Poisson lifting.** If $\Lambda\le_{\rm cx}\Lambda'$ are nonnegative and $C,C'$ are mixed Poisson with these intensities, then $C\le_{\rm cx}C'$.

*Proof.* The means agree. For an integer $q\ge0$, let $h_q(x)=\mathbb E(P_x-q)_+$, where $P_x$ is Poisson($x$). Differentiating the Poisson series term by term gives $h_q'(x)=\mathbb P(P_x\ge q)$. This derivative is nondecreasing in $x$, because Poisson laws increase stochastically in their mean. Hence $h_q$ is convex, and $\mathbb E(C-q)_+=\mathbb E h_q(\Lambda)\le\mathbb E h_q(\Lambda')=\mathbb E(C'-q)_+$. For integer-valued variables, $a\mapsto\mathbb E(C-a)_+$ is linear between consecutive integers, and it equals $\mathbb EC-a$ for $a\le0$. Integer thresholds therefore suffice. $\square$

This is the scalar case of the Cox-process ordering of Błaszczyszyn and Yogeshwaran (2009, Proposition 3.6), as used in the main text.

**Theorem D2 — ordering of the three versions.** Fix $N,K,\beta,\mu$, a sample $A$ and either regime. Then
$$\boxed{\Lambda_2\le_{\rm cx}\Lambda_1\le_{\rm cx}\Lambda_3,\qquad C^{(2)}_A\le_{\rm cx}C^{(1)}_A\le_{\rm cx}C^{(3)}_A,}$$
with common mean $m=\mu\sum_a\mathbb EY_a$. Consequently, for every threshold $q$,
$$\mathbb E(q-C^{(2)}_A)_+\le\mathbb E(q-C^{(1)}_A)_+\le\mathbb E(q-C^{(3)}_A)_+,$$
and likewise for $(C-q)_+$. The variances satisfy
$$\operatorname{Var}C^{(2)}_A=m,\quad \operatorname{Var}C^{(1)}_A=m+\operatorname{Var}\,\mathbb E[\mu L_A\mid G],\quad \operatorname{Var}C^{(3)}_A=\operatorname{Var}C^{(1)}_A+\mathbb E\operatorname{Var}(\mu L_A\mid G),$$
and the probabilities of an empty sample satisfy
$$e^{-m}=\mathbb P(C^{(2)}_A=0)\le\mathbb P(C^{(1)}_A=0)\le\mathbb P(C^{(3)}_A=0).$$

*Proof.* By definition, $\Lambda_1=\mathbb E[\Lambda_3\mid\mathcal F_1]$ and $\Lambda_2=\mathbb E[\Lambda_1\mid\mathcal F_2]$, where $\mathcal F_2\subseteq\mathcal F_1\subseteq\mathcal F_3$. A conditional expectation lies below its argument in convex order, by conditional Jensen. That gives the intensity chain, and Theorem D1 with Lemma D2 gives the count chain. The variance identities are the law of total variance applied to $\Lambda_3$ along $\mathcal F_1$. The empty-state probabilities are $\mathbb E e^{-\Lambda_v}$, and $x\mapsto e^{-x}$ is convex. $\square$

The garbling structure is the same as in the main text (§4), with two separate steps. Version 3 to Version 1 forgets which links succeeded and keeps who the potential models were. Version 1 to Version 2 forgets the potential models as well. The per-trait transmission law $(K,\beta)$ is held fixed throughout, as in Wakano and Aoki. No compensating change of encounter frequency is involved, unlike the $s=\beta/b$ convention of the continuous model.

**Long-run interpretation.** In Version 3, from any configuration the next generation is empty with probability at least $(1-\beta)^{NK}e^{-N\mu}>0$, so that chain is uniformly ergodic (Doeblin). In Versions 1 and 2 this bound fails for large repertoires (traits fail separately), so uniform ergodicity is not claimed; there, convergence from each finite initial configuration follows because every finite collection of initial traits goes extinct almost surely and the backward stationary construction then applies. The stationary law above is the limit of the finite-horizon law as $T\to\infty$, since $C^{(T)}\uparrow C^{(\infty)}$ almost surely under the backward coupling. (Correction, 3 Nov 2026: an earlier draft of this paragraph asserted uniform ergodicity for all three versions; an audit showed that this holds only for Version 3. The ordering theorem does not use it.)

## D.6 When the inequalities are strict

Strictness here means that the laws differ. For mixed Poisson counts this is equivalent to unequal laws of the intensities, since $\mathbb Ez^C=\mathbb Ee^{-\Lambda(1-z)}$. Along a conditional-expectation coupling it is equivalent to a strict variance inequality. Strictness then also holds for $\mathbb P(C=0)$, and for $\mathbb E\phi(C)$ whenever $\phi$ is strictly convex with finite expectations.

**Proposition D3.** Let $\mu>0$, $0<\beta<1$ and $n\ge1$.

1. $C^{(1)}_A$ and $C^{(3)}_A$ have different laws if and only if $T\ge1$.
2. $C^{(2)}_A$ and $C^{(1)}_A$ have different laws if and only if $K<N$ and either (i) $n\ge2$ and $T\ge1$, or (ii) $n=1$, $K\ge2$ and $T\ge2$.

In particular, in the stationary regime both inequalities are strict whenever $K<N$, except for a single sampled individual with $K=1$. In that exceptional case Versions 1 and 2 coincide in law.

*Proof.* (1) If $T=0$, then $L_A=n$ and all three counts are Poisson($n\mu$). Let $T\ge1$. For every $G$, all links of the sample fail with conditional probability $(1-\beta)^{nK}>0$, giving $L_A=n$. All succeed with probability $\beta^{nK}>0$, giving $L_A\ge n+K$. Thus $\operatorname{Var}(L_A\mid G)>0$ for every $G$.

(2, sufficiency.) Let $\ell(b)$ be the expected ancestral volume of generations $\le-1$ when $b$ ancestors are present at generation $-1$. By exchangeability, this depends only on $b$. It is strictly increasing in $b$, because ancestral sets are monotone in their initial set and generation $-1$ contributes $b$ itself. Let $G_0$ denote the sample's potential-model sets. Then
$$\mathbb E[L_A\mid G_0]=n+\mathbb E[\ell(|A_{-1}|)\mid G_0].$$
Given $G_0$, the indicators $\{k\in A_{-1}\}$ are independent with probabilities $f(c_k)$, where $f(c)=1-(1-\beta)^c$ and $c_k=\#\{i\in A:k\in G_0(i)\}$.

In case (i), compare two events of positive probability: $E_1$, on which every $G_0(i)=S$, and $E_2$, which is the same except that $G_0(i_2)=S\setminus\{x\}\cup\{y\}$ for one $i_2\in A$. The event $E_2$ requires $K<N$ and $n\ge2$. On $E_2$, the two indicators for $x,y$ have probabilities $f(n-1)$ and $f(1)$. Their union occurs with probability $1-(1-\beta)^n=f(n)$, while both occur with positive probability. Hence, under the obvious coupling, $|A_{-1}|$ is at least as large on $E_2$ as on $E_1$, and strictly larger with positive probability. Strict monotonicity of $\ell$ makes $\mathbb E[L_A\mid G_0]$ nonconstant. Finally, $\operatorname{Var}\Lambda_1\ge\operatorname{Var}\mathbb E[\Lambda_1\mid G_0]>0$.

In case (ii), condition on $(G_0,G_{-1})$ with $G_0$ fixed. Compare the events that all generation $-1$ potential sets equal $S$, and that they equal $S$ except at one individual $k^*\in G_0(i)$. Ancestor sets $B=A_{-1}$ with $k^*\notin B$ or $|B|=1$ give equal conditional means. Sets with $k^*\in B$ and $|B|\ge2$ have probability $\beta\{1-(1-\beta)^{K-1}\}>0$, and for them the argument of case (i) gives a strict increase. This requires $T\ge2$.

(2, necessity.) If $K=N$, then $G$ is deterministic and $\Lambda_1=m$. If $n=1$ and $K=1$, the ancestry is a single path whose survival to depth $a$ has probability $\beta^a$ for every $G$. If $n=1$ and $T=1$, then $|A_{-1}|\sim\mathrm{Bin}(K,\beta)$ for every $G$. In each case $\Lambda_1$ is constant. $\square$

## D.7 The whole age-frequency spectrum

For an age $0\le a\le T$ and a nonempty $S\subseteq A$, let $N_{a,S}$ count the traits of age $a$ whose carriers in $A$ are exactly $S$. Then $C_{a,n,i}=\sum_{|S|=i}N_{a,S}$ is the age-frequency spectrum of Kobayashi et al. (2018). The popularity class $P_i=\sum_aC_{a,n,i}$ counts traits held by exactly $i$ sample members. For a family $\mathcal S$ of pairs $(a,S)$, put $N_{\mathcal S}=\sum_{(a,S)\in\mathcal S}N_{a,S}$.

**Corollary D4 — every spectrum count.** Fix any family $\mathcal S$. In each version, conditional on $\mathcal F_v$, the counts $(N_{a,S})$ are independent Poisson variables. The intensities are $\Lambda_{v}(a,S)=\mathbb E[\mu\,\nu_{a,S}(H)\mid\mathcal F_v]$, where $\nu_{a,S}(h)=\#\{j:\{i\in A:j\in A^{\{i\}}_{-a}(h)\}=S\}$. Consequently:

1. $N^{(2)}_{\mathcal S}\le_{\rm cx}N^{(1)}_{\mathcal S}\le_{\rm cx}N^{(3)}_{\mathcal S}$, with equal means. This holds for each age-frequency entry $C_{a,n,i}$, each popularity class $P_i$, the age spectrum $C_{a,n}$, the number of traits held by at least $q$ members, and the number held by one subsample and not another.
2. In Version 2 all $N_{a,S}$ are independent Poisson variables.
3. For any finite collection of classes, the covariance matrices are ordered in the positive-semidefinite order, $\operatorname{Cov}^{(2)}\preceq\operatorname{Cov}^{(1)}\preceq\operatorname{Cov}^{(3)}$. Hence $\operatorname{Var}\sum_xw_xN_x$ is ordered in the same direction for every real weight vector $w$.

*Proof.* Mark each trait by its age and sample carrier set. The proof of Theorem D1 applies to the marked Poisson process and gives (2) and the conditional independence. Item (1) is Theorem D2's proof with $L_A$ replaced by $\sum_{\mathcal S}\nu_{a,S}$. For (3), $\operatorname{Cov}^{(v)}=\operatorname{diag}(\mathbb E\Lambda)+\operatorname{Cov}(\Lambda_v)$, where the diagonal term is the same in every version. Along the martingale $\Lambda_1=\mathbb E[\Lambda_3\mid\mathcal F_1]$, $\operatorname{Cov}(\Lambda_3)-\operatorname{Cov}(\Lambda_1)=\mathbb E\operatorname{Cov}(\Lambda_3\mid\mathcal F_1)\succeq0$, and $\operatorname{Cov}(\Lambda_1)\succeq0=\operatorname{Cov}(\Lambda_2)$. $\square$

Item (1) is a scalar statement for each fixed linear count. It does not assert a multivariate convex order for the whole spectrum vector. The Poisson-lifting step need not be jointly convex for general convex functions of several counts, and we make no such claim. Item (3) is the precise form of the "greater variance in the cultural spectrum" that Wakano and Aoki anticipate for Version 3 (§1.6.2). It also shows that Version 1 lies between Versions 2 and 3.

## D.8 Partial bundling and the role of the innovation law

**A one-parameter interpolation.** Let $U_\ell$ be independent Uniform$(0,1)$ labels attached to the links $\ell$ of $G$, independent of everything else. Fix $\rho\in[0,1]$. A link with $U_\ell\le\rho$ is a *bundled* link: one success indicator $\xi_\ell$ is used for all traits, so the model teaches everything or nothing. The other links decide success trait by trait. Each trait still has success probability $\beta$ on every link. With $\mathcal F^{\rho}=\sigma(G,U,(\xi_\ell\mathbf 1\{U_\ell\le\rho\})_\ell)$, the intensity is $\Lambda^{\rho}=\mathbb E[\mu L_A(H)\mid\mathcal F^\rho]$. The σ-fields are nested, $\mathcal F^{\rho}\subseteq\mathcal F^{\rho'}$ for $\rho\le\rho'$, and $\rho=0,1$ give Versions 1 and 3. A similar family with per-newborn labels $V_{t,i}$, where the potential set is shared when $V_{t,i}\le\theta$ and redrawn for each trait otherwise, joins Version 2 ($\theta=0$) to Version 1 ($\theta=1$).

**Proposition D5.** $C_A$ increases in convex order along both families, from Version 2 through Version 1 to Version 3, with fixed mean.

*Proof.* Theorem D1 holds with $\mathcal F^\rho$ in place of $\mathcal F_v$, because traits are conditionally independent copies of $H$ given $\mathcal F^\rho$. The independent labels do not alter the conditional law of $H$ beyond the revealed indicators. Nesting gives $\Lambda^\rho=\mathbb E[\Lambda^{\rho'}\mid\mathcal F^\rho]$. Conclude as in Theorem D2. $\square$

In the discrete model, $\rho$ plays the role of bundle completeness: the share of teaching relationships that pass the whole repertoire at once.

**Non-Poisson innovation.** The Poisson assumption is needed for the comparison with Version 2, but not for the comparison of Versions 1 and 3.

**Proposition D6.** Let the innovation numbers $(\Pi_{s,j})$ have any joint law with uniformly bounded means, independent of all transmission randomness. Then $C^{(1)}_A\le_{\rm cx}C^{(3)}_A$, and the same holds for every count of traits carried by at least $q$ sample members. However, $C^{(2)}_A\le_{\rm cx}C^{(1)}_A$ can fail, and so can $C^{(2)}_A\le_{\rm cx}C^{(3)}_A$.

*Proof.* First take a finite horizon, and condition on $G$ and $\Pi$. The link indicators of $H$ are independent, and each trait's presence indicator $\mathbf 1\{j_\tau\in A_{s_\tau}(H)\}$ is a nondecreasing function of them. By Esary, Proschan and Walkup (1967), the finite vector of indicators is associated. Version 1 uses independent copies of the link indicators for each trait, so its vector is the *independent duplicate* of Version 3's vector. Christofides and Vaggelatou (2004) prove that an associated vector dominates its independent duplicate in supermodular order. Since $x\mapsto\phi(\sum_ix_i)$ is supermodular for convex $\phi$, the sums are ordered in convex order conditionally, and hence after mixing. The stationary case follows by monotone $L^1$ convergence as $T\to\infty$; convex order is closed under $L^1$ limits. The same argument applies to "carried by at least $q$ members", which is also a nondecreasing event.

For the counterexample, take $N=3$, $K=1$, $\beta=9/10$, $T=3$, $n=1$, and let each newborn invent one trait with probability $\mu=1/2$ and none otherwise. With $p_a=\mu\beta^a$:

- Version 1 follows one potential-model path, fixed by $G$, and its depth-$a$ indicators are independent Bernoulli($p_a$) variables. Thus $\operatorname{Var}C^{(1)}=\sum_ap_a(1-p_a)$.
- In Version 2, for each depth $a\ge1$, every one of the three innovators lies on the trait's own random path with probability $1/3$. Thus $\operatorname{Var}C^{(2)}=p_0(1-p_0)+\sum_{a\ge1}3\frac{p_a}{3}(1-\frac{p_a}{3})$.

The exact values (`bernoulli_counterexample.py`) are $\operatorname{Var}C^{(1)}=0.970115<\operatorname{Var}C^{(3)}=1.116320<\operatorname{Var}C^{(2)}=1.303038$. The common mean is $1.7195$. $\square$

Thus the ranking of Version 2 is a property of Poisson innovation, not of the genealogies alone. A Poisson total with a given mean is the convex-order minimum among mixed Poisson laws. A sum of Bernoulli variables is less spread when its probability is concentrated on fewer innovators, as it is in Version 1. If Kobayashi et al. (2018) specify a non-Poisson innovation number, then only the Version 1 $\le_{\rm cx}$ Version 3 part of Theorem D2 transfers as stated.

## D.9 Exact moments and numerical checks

**Exact variance.** Take two conditionally independent copies $H,H'$ of the single-trait graph given $\mathcal F_v$. Then $\operatorname{Var}C^{(v)}_A=m+\mu^2(\mathbb E[L_A(H)L_A(H')]-(\mathbb EL_A)^2)$. Both ancestral sets evolve together as a finite Markov chain on pairs of subsets of $[N]$:

- in Version 1, a newborn ancestral to both replicas uses one potential set and two independent thinnings;
- in Version 3, it uses one thinning;
- in Version 2, the two replicas are independent.

First-step analysis then gives a finite linear system. For each horizon it is a backward recursion, as in §6 of the main text. Exchangeability reduces the system to counts of sites in the first set only, the second only, and both. For Version 3, the single chain $Y$ on $\{0,\ldots,N\}$ suffices.

**Checks** (all single-threaded; scripts and outputs in this folder).

1. `exact_versions.py`, with exact rational arithmetic, for $N=3$, $K\in\{1,2,3\}$, $\beta\in\{1/2,4/5,9/10\}$, $\mu=1/2$, samples $n=1,2,3$, horizons $T=0,\ldots,4$ and the stationary regime. In every case the means agree, $\operatorname{Var}C^{(2)}=m$, and $\operatorname{Var}C^{(2)}\le\operatorname{Var}C^{(1)}\le\operatorname{Var}C^{(3)}$. Equalities occur exactly where Proposition D3 predicts them:
   - all versions coincide at $T=0$;
   - Versions 1 and 2 coincide for $n=1,T=1$, for $n=1,K=1$ (all $T$), and for $K=N=3$;
   - every other case is strict.

   An independent brute-force enumeration of all $3^{3T}$ potential-model assignments ($T\le3$) reproduces the Version 1 variance to $10^{-10}$. It also gives $\mathbb P(C=0)$ ordered as in Theorem D2.
2. `ya_chain.py` recomputes Version 2 and Version 3 from the lumped count chain $Y$ alone and agrees with (1).
3. `mc_versions.py` is a forward simulation of explicit trait identities, generation by generation, with no use of the theorem. It was run with $1.5\times10^6$ replicates per version for $(N,K,\beta,\mu,T)=(3,2,4/5,1/2,3)$ and $(3,1,9/10,1/2,3)$. The second case was rerun with a second seed.

Exact and simulated values for $(N,K,\beta,\mu,T)=(3,2,4/5,1/2,3)$, with Monte Carlo standard errors in parentheses:

| Sample | Mean (exact) | Var V2 exact / MC | Var V1 exact / MC | Var V3 exact / MC |
|:--|--:|--:|--:|--:|
| $n=1$ | 3.4103 | 3.4103 / 3.413 (0.004) | 3.4377 / 3.438 (0.004) | 4.3075 / 4.302 (0.005) |
| $n=2$ | 4.6292 | 4.6292 / 4.626 (0.006) | 4.7214 / 4.727 (0.006) | 5.1268 / 5.117 (0.006) |
| $n=3$ | 5.4093 | 5.4093 / 5.402 (0.007) | 5.4887 / 5.492 (0.007) | 5.7517 / 5.739 (0.007) |

For the second case, $(3,1,9/10,1/2,3)$ with $n=2$, the exact variances are $2.8953<3.0980<3.3589$. The two simulations gave $2.891,3.087,3.355$ and $2.894,3.097,3.358$, each with standard errors of about $0.004$. For the full population, the exact values of $\mathbb P(C=0)$ are $0.0219<0.0252<0.0284$ (simulated $0.0218,0.0252,0.0282$).

All simulated means lie within $3$ standard errors of the exact means. Most simulated variances lie within $2$ standard errors of the exact values, and none lies more than $3$ away. The only excursions near $3$ standard errors are in Version 1 of the first $K=1$ run, where they are correlated across the nested samples. The rerun with a second seed lies within $1$ standard error of the exact values. The simulated popularity classes $P_1,\ldots,P_n$ have equal means across versions and variances ordered $V2<V1<V3$ in every case (Corollary D4).

With at most one innovation per newborn, simulation reproduces the counterexample of Proposition D6, giving $0.971$, $1.303$ and $1.118$ for Versions 1, 2 and 3. For the $K=2$ case it also shows $\operatorname{Var}C^{(1)}<\operatorname{Var}C^{(2)}$ clearly at $n=1,2$, and by about $2.7$ standard errors at $n=3$. Version 1 $<$ Version 3 in every case, as Proposition D6 requires.

**Size of the gaps.** In these tiny populations, Version 1 lies close to the Poisson floor, and Version 3 can be far above it. For $N=3$, $K=2$, $\beta=1/2$, $\mu=1/2$ and the stationary full population, the variances are $8.08$, $8.17$ and $43.0$. For Wakano and Aoki's illustrative setting ($K=2$, $\beta=1/2$, $N=150$, $n=10$, $\mu=1$; their §1.6.2), `ya_chain.py` gives $\mathbb EC=327.0$, consistent with their figure-read "~300". The standard deviation of Version 2 is $18.1$ and that of Version 3 is $315.6$. Theorem D2 places Version 1 between these, but its value is not computed here; it requires the $O(N^3)$ pair system or simulation. No large-$N$ statement is inferred from these numbers.

## D.10 What this settles and what it does not

*Settled, within the stated model.* The question of Wakano and Aoki (2024, §1.6.2) has an exact qualitative answer: at equal per-trait transmission, sharing more of the transmission randomness across traits increases the repertoire count in convex order. The ordering is Version 2 $\le_{\rm cx}$ Version 1 $\le_{\rm cx}$ Version 3. It covers every sample size, every finite horizon from empty and the stationary state, and every age-frequency class count. The hardest case, Version 1, is bracketed by the two cases they call potentially solvable:

- Version 2 is exactly Poisson with mean $\mu\sum_a\mathbb EY_a$;
- Version 3 is mixed Poisson with intensity $\mu\sum_aY_a$, and its law is computable from the single chain $Y$.

The expected spectra are identical across versions. Spectrum means alone therefore cannot distinguish them, although the spectrum's covariance can. Strictness is characterised completely (Proposition D3). Partial bundling interpolates monotonically (Proposition D5).

*Not settled.*

- There is no closed form for the law of Version 1. Its intensity is the conditional expectation $\mu\,\mathbb E[L_A\mid G]$ over a random graph, and only its moments are given by finite systems.
- There is no asymptotic result for large $N$, and no quantification of how close Version 1 is to the Poisson floor in realistic populations.
- No statement about temporal autocorrelation or episode durations is made beyond the stationary marginal law and Kac's mean recurrence time to complete absence.
- The Version 2 comparisons depend on Poisson innovation (Proposition D6). Whether Kobayashi et al. (2018) assume Poisson innovation could not be checked against their full text.
- The model keeps the neutral, independent-transmission assumptions of the source: no prerequisites, no selection, no biased choice of potential models, and fixed $N$. Nonrandom choice of models, which Wakano and Aoki raise in their §1.7 and §3.2, is outside the result unless it can still be represented as a shared graph with conditionally independent trait marks.

## References

Błaszczyszyn, B., & Yogeshwaran, D. (2009). Directionally convex ordering of random measures, shot noise fields and some applications to wireless communications. *Advances in Applied Probability, 41*(3), 623–646. https://doi.org/10.1239/aap/1253281057

Christofides, T. C., & Vaggelatou, E. (2004). A connection between supermodular ordering and positive/negative association. *Journal of Multivariate Analysis, 88*(1), 138–151. https://doi.org/10.1016/S0047-259X(03)00064-2

Esary, J. D., Proschan, F., & Walkup, D. W. (1967). Association of random variables, with applications. *The Annals of Mathematical Statistics, 38*(5), 1466–1474. https://doi.org/10.1214/aoms/1177698701

Kobayashi, Y., Wakano, J. Y., & Ohtsuki, H. (2018). Genealogies and ages of cultural traits: An application of the theory of duality to the research on cultural evolution. *Theoretical Population Biology, 123*, 18–27. https://doi.org/10.1016/j.tpb.2018.04.007 (full text not obtained; see §D.1)

Shaked, M., & Shanthikumar, J. G. (2007). *Stochastic Orders*. Springer Series in Statistics. Springer. https://doi.org/10.1007/978-0-387-34675-5

Wakano, J. Y., & Aoki, K. (2024). Forward and backward modeling of cultural evolutionary processes. *Proceedings of the National Academy of Sciences, 121*, e2322880121. https://doi.org/10.1073/pnas.2322880121

Wakano, J. Y., Ohtsuki, H., Kobayashi, Y., & Baake, E. (2026). Evolving genealogies in cultural evolution, the descendant process, and the number of cultural traits. *Theoretical Population Biology, 168*, 1–18. https://doi.org/10.1016/j.tpb.2025.12.003

## Files in this folder

- `exact_versions.py`: exact pair-chain moments and brute-force enumeration.
  - Outputs: `exact_N3K2_b1o2_mu1o2.txt`, `exact_N3K2_b4o5_mu1o2.txt`, `exact_N3K1_b9o10_mu1o2.txt`, `exact_N3K3_b1o2_mu1o2.txt`.
  - Usage: `exact_versions.py N K beta mu`.
- `ya_chain.py`: lumped $Y$-chain for Versions 2 and 3. Output: `ya_chain_N150_K2_b1o2_mu1_n10.txt`.
- `mc_versions.py`: forward trait-identity simulation.
  - Usage: `mc_versions.py N K beta mu T reps poisson|bernoulli seed`.
  - Outputs: `mc_*_poisson*.txt` and `mc_*_bernoulli.txt`.
- `bernoulli_counterexample.py`: the exact counterexample of Proposition D6. Output: `bernoulli_counterexample_output.txt`.
- `tools_html2txt.py` and `sources_dl/`: the downloaded Wakano–Aoki (2024) and Wakano et al. (2026) texts used for the locators.
