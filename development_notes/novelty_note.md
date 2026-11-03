# Prior-art gate: partial cultural transmission and distributional order

**Audit date:** 2 November 2026. **Status:** conditional proceed for a narrowly stated cultural-model result; no certification of priority. This note separates what the retrieved papers prove from the candidate result supplied by the research team. No empirical validation is implied.

## 1. Object identity and the genuinely relevant question

Wakano, Ohtsuki, Kobayashi and Baake (2026), *Evolving genealogies in cultural evolution, the descendant process, and the number of cultural traits*, TPB 168, 1–18, [DOI](https://doi.org/10.1016/j.tpb.2025.12.003), [full version 2](https://arxiv.org/html/2505.08677v2), specify:

- Fixed population size (N); death/replacement by an empty individual at rate (u>0) per individual.
- Directed learning encounters at rate (s/N) per ordered donor/learner pair; self-encounters are silent.
- Unique, nonrecurrent innovations at rate μ per individual.
- Baseline transmission copies the donor’s entire repertoire, merging it with the learner’s repertoire.

Their equation (15) is already the mixed-Poisson representation (C_n\mid L_n\sim\mathrm{Poisson}(\mu L_n)). Section 4 changes copying to independent Bernoulli((b)) success for **each trait within each shared encounter**, and replaces (s) by (s/b). It explores this modification through simulations, not an analytical distributional comparison.

For the proposed extension, use a fresh symbol β for the **effective per-trait learning rate** and write the raw encounter rate as β/(b). Holding β fixed is essential. Independent per-trait copying outcomes do **not** make complete trait histories independent: their potential encounters and carrier deaths remain shared.

Candidate object: the stationary number (C_b) of distinct, present traits in a specified sample. Candidate claim supplied by the team: (b_1<b_2) at fixed ((N,u,\mu,\beta)) implies (C_{b_1}\leq_{cx}C_{b_2}), while expected counts and expected popularity spectra coincide. This is a statement about dependence between traits, not a change in the marginal transmission process of one trait.

## 2. Explicit open-problem context

Wakano and Aoki (2024), *Forward and backward modeling of cultural evolutionary processes*, PNAS 121, e2322880121, [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11621620/), §§1.6.1–1.6.2, distinguish three discrete-generation versions:

| Version | Potential models and actual successful parents |
|---|---|
| 1 | Potential models shared across traits; successes independent by trait. |
| 2 | Potential models independently sampled for each trait. |
| 3 | Actual successful cultural parents shared across all traits. |

They identify distributions beyond means as an important unresolved task and describe Version 1 as the difficult case because genealogies are correlated. This is unusually direct motivation. However, the proposed continuous-time model is an analogue of Version 1, not the identical Wright–Fisher model. A theorem for it must not be advertised as solving the entire discrete-generation problem. Expected genealogical lengths and expected trait spectra were already known in this literature.

## 3. Closest collisions: techniques that cannot be claimed as new

### 3.1 Shared transmission changes higher moments: explicitly anticipated in 2014

Baumdicker and Pfaffelhuber (2014), *The infinitely many genes model with horizontal gene transfer*, [full text](https://arxiv.org/html/1301.6547v2), [DOI](https://doi.org/10.1214/EJP.v19-2642), **Remark 2.4**, already propose common transfer/loss events with per-gene probabilities (p_\gamma,p_\rho). They state that preserving γ=γ′(p_\gamma) and ρ=ρ′(p_\rho) preserves their first-moment results, while second moments change. This is closer prior art than the baseline model alone and rules out a broad novelty claim about unchanged means but altered variability.

Their §6.3, equations (6.6)–(6.11), also use two genealogical lengths and overlap counts (x) (gene 1 only), (y) (both), (z) (gene 2 only). Equation (6.8) supplies coupled first-step moment equations and (6.11) a small-transfer-rate covariance expansion.

**Boundary:** their analysed baseline includes clonal reproduction and gene-specific loss/transfer; the multi-gene variation is described, not solved there. Our possible advance is the exact order throughout the transmission family and any proved identification consequence. Paired replicas, overlap triples and the statement that shared events affect variance are established.

Earlier related conditioning: Baumdicker, Hess and Pfaffelhuber (2010), *The diversity of a distributed genome in bacterial populations*, Annals of Applied Probability 20, 1567–1606, [full paper](https://arxiv.org/pdf/0907.2572), [DOI](https://doi.org/10.1214/09-AAP657), Proposition 5.5 and subsequent propositions, derive conditional Poisson gene-presence patterns on a shared genealogical tree and count covariances. Conditional Poisson marking itself is not novel.

### 3.2 Poissonisation preserves the relevant order: already explicit in 2009

Błaszczyszyn and Yogeshwaran (2009), *Directionally convex ordering of random measures, shot noise fields and some applications to wireless communications*, Advances in Applied Probability 41, 623–646, [full paper](https://arxiv.org/pdf/0806.3180), Proposition 3.6, show that directionally convex order of random intensity measures passes to the associated Cox processes. Its scalar case gives

\[
W_1\leq_{cx}W_2
\quad\Longrightarrow\quad
\mathrm{Poisson}(\mu W_1)\leq_{cx}\mathrm{Poisson}(\mu W_2).
\]

This is a standard lifting lemma, even if independently reproved in the manuscript. It does not establish the model-specific ordering of the cultural mixing variables.

### 3.3 Conditional means under added noise: standard reasoning, with a close 2026 source

Dendoncker, Denuit and Robert (2026), *Reverse-time scores and convex-order monotonicity for additive processes with Lévy clocks*, working paper ISBA 2026/23, [record](https://hdl.handle.net/2078.5/276938), [full PDF](https://research.dial.uclouvain.be/bitstreams/8f09d482-51e4-4bce-8669-d9ab9bef70a5/download), Proposition 4.2 and Corollary 4.3, establish reverse-martingale and convex-order relations for conditional means under increasing additive noise. Their proof uses conditional expectation and Jensen’s inequality.

The proposed contact history can similarly be viewed as a latent effective Poisson graph plus independent extra Poisson encounters. Thus an information-garbling argument must be presented as an application of standard conditional-expectation order, not a new general probability principle. This working paper treats real-valued additive observations and does not state the cultural genealogical theorem.

## 4. Other close models and why they do not presently settle the candidate

| Source inspected | Result or model | Exact distinction |
|---|---|---|
| Yeh, Fogarty and Kandler (2019), *Cultural linkage: the influence of package transmission on cultural dynamics*, [DOI](https://doi.org/10.1098/rspb.2019.1951), [full XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6939269/fullTextXML) | Agent-based model with traits linked into transmissible packages; links form and break; each attempted copy succeeds independently. Intermediate linkage can reduce pairwise diversity. | Finite, fixed trait dimensions with variants; no carrier birth/death. Package size affects marginal acquisition rates. Its nonmonotone diversity findings neither prove nor refute a fixed-β count-order theorem. This is important anthropological precedent for package transmission and observational equifinality. |
| Jahnel, Köppl, Lüchtrath and Vu (2026), *Convex order and faster transmission in first contact percolation*, [full text](https://arxiv.org/html/2605.28766v1), Definition 2.1 and Theorem 2.3 | Convex ordering of encounter point processes yields ordered infection-path counts and existence probabilities. | No recovery; changes the single-process contact law. Our comparison preserves each trait’s effective Poisson law and changes intertrait dependence. Do not claim first use of convex order in graphical transmission models. |
| McVinish (2025), *Re-examining the drivers of variation in parasite loads across hosts in the Tallis–Leyton model*, J Math Biol 91:21, [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12296923/), Theorems 2–4 | Convex/Lorenz comparisons under compound-Poisson contact rates, batch sizes and binomial thinning. | Parasite arrivals have exogenous independent batch sizes; parasites do not reproduce within hosts. Cultural repertoires feed back into later copying opportunities. Its results are close methodological analogues, not the same object. |
| Scheie, Jansson and Sjöstrand (2026), *Sequential accumulation of independent cultural traits: A Gaussian limit for the popularity spectrum*, [PubMed abstract](https://pubmed.ncbi.nlm.nih.gov/42918355/), [DOI](https://doi.org/10.1016/j.tpb.2026.10.003) | Published 23 October 2026. Repeated copying and forgetting among long-lived individuals; stationary popularity spectrum and Gaussian asymptotics. | **Full text not obtained.** The abstract explicitly distinguishes its mechanism from turnover through birth/death and reports coincidence with Nakamura et al.’s expected spectrum. It does not mention a bundling/convex-order theorem, but abstract silence cannot exclude one. Priority remains subject to this access limitation. |

## 5. Narrow claim that survives this search

The potentially original claim is an **exact comparison across partial-transmission probabilities in the specified cultural accumulation model**, at a fixed marginal transmission rate: the mixing variable and total count become more dispersed in convex order as encounters carry traits more jointly. No retrieved cultural source states that comparison. A conditional-Poisson construction for Bernoulli transmission supplies the bridge from a shared history to this result; the stochastic-order principles and genealogical replica machinery are established.

Anthropological value would lie in showing that identical expected trait-frequency summaries can coexist with systematically different dispersion and complete-absence probabilities. This is a precise answer to a current model question about the consequences of correlated genealogies. It does not by itself show that the effect explains any historical cultural loss, or that repertoire counts measure cultural value or complexity.

Required wording limits:

1. At fixed β, lower (b) is compensated by more encounters. Do not describe the theorem as a benefit of poorer copying at unchanged opportunities.
2. (P(C_b=0)) itself measures stationary absence in the specified sample/population. The additional full-population regeneration argument now proves the precisely defined temporal consequences in Section 7 below. These do not cover arbitrary initial repertoires, arbitrary loss thresholds, sample-empty episodes, or post-shock recovery.
3. Convex order controls convex loss functions and stop-loss transforms. It does not order arbitrary upper/lower tail probabilities or establish uniformly smaller collapse hazards.
4. Distinguish a conditional-Poisson characterization from a tractable closed-form distribution. An unevaluated random mixing measure is not a complete numerical solution.
5. No quality judgement follows from novelty alone. The short generic proof may be mathematically elegant, but significance depends on the scope of the cultural conclusion and evidence that it changes what researchers can infer or calculate.

## 6. Search and access record

Searches used both search engines on 2 November 2026, including exact titles and combinations of cultural transmission, correlated genealogies, partial transmission, variance, mixed Poisson, convex order, contact process, package transmission, Poisson thinning and information garbling. The core 2026 article, 2024 review, 2014 HGT article, 2009 ordering article and 2026 contact-percolation paper were opened as full text. The 2019 cultural-linkage and 2025 parasite papers were additionally retrieved through Europe PMC’s public XML service; scratch copies in `source_checks` are for verification, not an intended redistribution bundle.

The 2018 Kobayashi–Wakano–Ohtsuki paper was located and its abstract plus its model’s detailed treatment in the 2024 paper inspected, but its full text was not obtained. The October 2026 Scheie et al. article was searched by title, DOI, authors, arXiv and institutional-repository terms; publisher retrieval failed and no lawful full-text copy was found in this audit. A follow-up inspected Jansson’s [current author publications page](https://www.fredrik.name/publications), where the new paper appears without a PDF/preprint link; Crossref metadata indicate CC BY 4.0, but OpenAlex supplies no repository full-text location. Public Elsevier endpoints returned a site-unavailable response. For Kobayashi (2018), the author publication pages at [Kochi](https://www.kochi-tech.ac.jp/profile/en/kobayashi-yutaka.html), [Wakano’s site](https://joewakano.sakura.ne.jp/research/) and [Ohtsuki’s site](https://rcies.soken.ac.jp/labs/ohtsuki/publications_e.html) lead only to the publisher; Crossref/OpenAlex indicate closed access without a repository copy. Metadata responses are saved under `source_checks`. These are substantive scope limits. This search is an informed prior-art screen, not an exhaustive proof of absence.

## 7. Targeted rerun: variance shape and complete-empty episodes

**Date:** 2 November 2026. This rerun covers the two subsequent corollaries, not a new exhaustive review. Their exact statements and exclusions are frozen in `target.yaml`, `math_results.md` §§9–10 and manuscript §§4.2–5. The main convex-order finding and the two previously recorded access gaps remain unchanged.

### 7.1 Variance shape: established Poisson identities, a model-specific application

Conditional on the reset history $D$, the two effective arrow histories have common intensity $b\Lambda$ and independent private intensities $(1-b)\Lambda$. Their covariance is therefore the Poisson noise-stability function $\langle R_D,P_bR_D\rangle$. Last and Penrose (2011), **Theorem 1.1, equation (1.5)**, supply its chaos isometry; Last, Peccati and Schulte (2016), **§3, equation (3.1), Lemma 3.1 and equation (3.6)**, give the thinning/immigration operator and its multiplier $b^k$ on order-$k$ kernels. Both author manuscripts were opened and the displayed equations inspected.

Consequently the nonnegative expansion

$$V(b)=V(0)+\sum_{k\geq1}a_kb^k,\qquad
a_k=\frac1{k!}\mathbb E_D\|\mathbb E_H\Delta^kR_D\|^2\geq0$$

is a direct application of established probability. The coefficient normalisation includes innovation rates inside $R_D$. At $N=2,\beta=2,u=1,\mu=1/2$, its exact fingerprint is

$$V(b)=\frac74+\frac{9}{2(5-2b)}
=\frac{53}{20}+\frac9{10}\sum_{k\geq1}(2b/5)^k.$$

This yields $V(0),V(1/4),V(1/2),V(1)=53/20,11/4,23/8,13/4$. The possible contribution is identifying this noise-stability function inside the cultural genealogy and establishing its model assumptions and strictness. The nonnegative-series principle, the Mehler operator, and the resulting convexity are **not new general theorems**.

### 7.2 Temporal results: known renewal formulas and a narrower new comparison

Use the **whole population**, total innovation rate $\alpha=\sum_i\mu_i>0$, stationary empty probability $\pi_{0,b}$, an exponential empty interval $J$ of rate $\alpha$, and the ensuing nonempty interval $B_b$. Innovations continue throughout $B_b$. The formulas

$$\text{entry rate}=\alpha\pi_{0,b},\qquad
\mathbb EB_b=\frac{1-\pi_{0,b}}{\alpha\pi_{0,b}},\qquad
\mathbb E(J+B_b)=\frac1{\alpha\pi_{0,b}}$$

are standard stationary-flow and renewal-reward identities. An inspected primary locator is Daley, Nazarathy and Wang (2026), **§3.2, equations (3.14)–(3.15)**, which explicitly records the cycle and busy-period means under the substitutions $\lambda_0=\alpha$ and $\tau_1=B_b$. Their object is a thinned count of transitions of a birth–death chain, rather than the stationary repertoire count; this paper does not supply the cultural fixed-$\beta$ ordering. It illustrates the established identities and is not claimed to originate them.

Likewise,

$$\widehat p_{00,b}(z)=
\frac1{z+\alpha-\alpha\mathbb E e^{-zB_b}}$$

is elementary alternating renewal. Pointwise order of empty-initial $p_{00,b}(t)$ then gives the stated Laplace-transform order of $B_b$. The prospective contribution is the cultural model's finite-horizon empty ordering and its particular consequences across $b$, not the renewal identity or the algebraic implication. Laplace-transform order does not establish stochastic order or ordered instantaneous hazards. The initial innovation has site probabilities $\mu_i/\alpha$; an arbitrary established repertoire is outside this experiment.

### 7.3 Closest temporal cultural precedents inspected

Wakano et al. (2026), **§3.3 equation (29) and §3.5 equations (32)–(35)**, already analyse ancestral/descendant extinction and interpret sawtooth genealogy lengths. They explicitly identify the descendant count with a single-trait SIS law. A descendant cohort's extinction and an innovating population's return to no traits are different events: subsequent innovations may remain. Their §4 still treats compensated partial copying through simulations. Thus this paper must be credited for the existing temporal analysis; the proposed addition is the comparison of the specifically defined complete-empty episodes.

Kobayashi, Kurokawa, Ishii and Wakano (2021), **§2.1; equations (5)–(8); §4.1**, study one initially present trait under best-of-$K$ transmission, with zero absorbing. Their published full text was inspected through a public [ResearchGate copy](https://www.researchgate.net/publication/348612891_Time_to_extinction_of_a_cultural_trait_in_an_overlapping_generation_model). It gives extinction distributions and means, including a branching approximation. This is neither a continuing-innovation multi-trait excursion nor a comparison changing only intertrait dependence. No matching fixed-marginal bundle result was found in these inspected passages.

### 7.4 Queries, adjudication and permitted wording

Both search engines were used. Exact query families included: `"cultural transmission" "busy period"`; `"cultural traits" "bundling" "extinction"`; `"cultural evolution" "Poisson" "Mehler"`; `"cultural traits" "extinction" "Laplace"`; `"cultural transmission" "empty" "renewal"`; `"cultural" "package transmission" variance Poisson`; `"cultural transmission" ("regenerative" OR "busy period")`; `"cultural traits" "extinction" "partial transmission"`; `"cultural" "Poisson" "variance" "Mehler"`; `"horizontal gene transfer" "busy period"`; `"cultural" "absolute monotonicity" variance`; `"cultural transmission" "variance" "53/20"`; and `"23-2b" "5-2b" variance`. Relevant title queries then retrieved the 2021 extinction and 2026 queueing papers. Unrelated search results were not treated as evidence. No directly matching cultural corollary was located; retrieval silence is not proof of priority.

| Contribution unit | Gate outcome | Defensible description |
|---|---|---|
| General nonnegative Poisson variance series | **COLLISION / established** | Apply credited isometry and Mehler identities. |
| Cultural variance's dependence on copying completeness | **BRIDGE, with bounded priority uncertainty** | Identify the ancestry functional and derive its variance shape and strictness. |
| Empty-entry and mean-excursion formulas | **COLLISION / established** | Apply stationary flow and renewal reward to the full population. |
| General alternating-renewal transform identity | **COLLISION / established** | Supply the short identity for clarity, without a novelty claim. |
| Fixed-marginal cultural empty-initial and excursion-transform comparisons | **BRIDGE, with bounded priority uncertainty** | Corollaries of the model-specific graphical ordering; no matching result found in this bounded search. |

Permitted: “The graphical comparison also yields variance-shape and complete-empty-episode consequences through established Poisson and renewal identities.” Forbidden: “A new renewal theorem”, “the first theory of cultural extinction”, or an unqualified priority claim. These corollaries strengthen the model's interpretation; they do not independently establish empirical significance or a REF rating.

### References added in this rerun

- Daley, D. J., Nazarathy, Y., & Wang, J. (2026). Second order properties of thinned counts in finite birth–death processes. *Queueing Systems, 110*, Article 17. https://doi.org/10.1007/s11134-026-09974-5
- Kobayashi, Y., Kurokawa, S., Ishii, T., & Wakano, J. Y. (2021). Time to extinction of a cultural trait in an overlapping generation model. *Theoretical Population Biology, 137*, 32–45. https://doi.org/10.1016/j.tpb.2021.01.002
- Last, G., & Penrose, M. D. (2011). Poisson process Fock space representation, chaos expansion and covariance inequalities. *Probability Theory and Related Fields, 150*(3–4), 663–690. https://doi.org/10.1007/s00440-010-0288-5. [Inspected manuscript](https://arxiv.org/pdf/0909.3205).
- Last, G., Peccati, G., & Schulte, M. (2016). Normal approximation on Poisson spaces: Mehler's formula, second order Poincaré inequalities and stabilization. *Probability Theory and Related Fields, 165*(3–4), 667–723. https://doi.org/10.1007/s00440-015-0643-7. [Inspected manuscript](https://arxiv.org/pdf/1401.7568).
