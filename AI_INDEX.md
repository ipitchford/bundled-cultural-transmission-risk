# AI index: equal mean diversity, unequal risk (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the
paper [output/cultural_transmission_paper.pdf](output/cultural_transmission_paper.pdf) (source [manuscript.md](manuscript.md)).
DOI: https://doi.org/10.5281/zenodo.23933211 (concept DOI 10.5281/zenodo.23933210).
Release page: https://evidencepress.org/releases/bundled-cultural-transmission-risk/

## Question

Continuous-time cultural trait accumulation model of Wakano, Ohtsuki, Kobayashi and Baake (2026): N sites, resets
(deaths) at rate u_i, unique innovations at rate mu_i, directed potential encounters at rate a_ji/b; at each encounter
every donor trait is copied independently with probability b. The per-trait effective rate a_ji (beta = s*b in the
homogeneous model; their s) is held fixed while b varies in [0, 1]; b = 0 is the independent-clock limit, not "no
learning". Question: how does the law of trait counts (repertoire C_A of a group A, popularity counts P_k, P_{>=k})
depend on b?

## Claims and where they are proved

Numbering follows [manuscript.md](manuscript.md). Prior-art audit: [development_notes/novelty_note.md](development_notes/novelty_note.md).

| ID | Claim | Location | Evidence type |
|---|---|---|---|
| C1 | Theorem 1: given the potential-encounter and reset history, counts for disjoint carrier-set families are independent Poisson with mean M_{S,b}; means E N_{S,b} do not depend on b | Section 3 | written proof (Poisson marking) |
| C2 | Theorem 2: for b1 <= b2, M_{S,b1} <=cx M_{S,b2} and N_{S,b1} <=cx N_{S,b2} (every family S, including C_A, P_k, P_{>=k}); stationary and empty-start finite horizon | Section 4.2 | written proof (superposition coupling, conditional Jensen, Poisson lifting of Blaszczyszyn-Yogeshwaran 2009 Prop. 3.6) |
| C3 | Corollaries 1-2: variance, expected shortfall E(q-N)_+, expected excess and P(N=0) nondecreasing in b; covariance matrix of disjoint counts nondecreasing in the PSD order | Section 4.2 | written proof |
| C4 | Section 4.4: threshold probabilities P(C <= q) are NOT ordered (exact N=2 counterexample table) | Section 4.4 | exact computation ([code/exact_extras.py](code/exact_extras.py), [results/exact_extras.json](results/exact_extras.json)) |
| C5 | Proposition 1: Var C_{A,b} (equivalently P(C_{A,b}=0)) strictly increasing iff some inventing site has a directed path of length >= 1 into A; otherwise the law of C_A does not depend on b | Section 4.3, Appendix D | written proof; exact chain example |
| C6 | Remark 1: at fixed encounter rate, counts of upward-closed families increase in increasing convex order | Section 4.2 | written proof sketch |
| C7 | Section 4.5: V(b) = V(0) + sum a_k b^k with a_k >= 0 (increasing, convex, chord bound); Lemma B1 b -> 0 limit | Section 4.5, Appendices B-C | written proof |
| C8 | Theorem 3: bundled forgetting (loss events at rate u/d deleting each trait with probability d): ordered in b and d separately; Poisson at b = d = 0 | Section 5 | written proof; exact grid ([code/partial_forgetting_check.py](code/partial_forgetting_check.py)) |
| C9 | Theorem 4: with Poisson innovation, discrete-generation Versions of Wakano and Aoki (2024): C^(2) <=cx C^(1) <=cx C^(3), equal means, all age/carrier counts; strictness in Appendix E; Version 2 position fails without Poisson innovation (exact counterexample) | Section 6, Appendix E | written proof; exact rational chains ([code/discrete/](code/discrete/)) |
| C10 | Section 7: size of the effect, V_N(b)/V_N(1) for beta/u in {0.7, 1, 1.3, 2}, N up to 150; e.g. N = 100, beta = 1.3: V(13/14)/V(1) = 0.336 | Section 7, Figure 2 | NUMERICAL (exact linear systems solved in floating point with residual checks; not interval-certified) |
| C11 | Corollary 3: complete-loss episodes (whole population): loss rate alpha*pi_0 increases and mean nonempty episode decreases with b | Section 8 | written proof (renewal) |
| C12 | Section 9: exact pair-replica linear system for second moments; subtraction-free covariance form | Section 9 | written derivation; independent forward closure agreement |
| C13 | Corollary 4: (E P_1, E P_2, Var C_V) identify (eta, tau, b) for known N; the expected spectrum alone carries no information about b | Section 10 | written proof |
| C14 | Section 11.1: N = 2 closed form; u = 1, beta = 2, mu = 1/2: mean 3/2 at every b, Var C = 53/20, 11/4, 23/8, 13/4 at b = 0, 1/4, 1/2, 1 | Section 11, Appendix A | exact rational computation |

## Verification artefacts

- [code/run_checks.py](code/run_checks.py): 87 checks (N = 2 rational pair system, labelled vs lumped generators,
  boundaries, b = 1 genealogy, re-check of four 100,000-cycle forward simulations). Output [results/verification.json](results/verification.json).
- [agent_math/pair_ancestry.py](agent_math/pair_ancestry.py): independently written backward implementation, 50 checks.
- [code/forward_closure_check.py](code/forward_closure_check.py): independent forward moment closure in exact rationals, 146 checks.
- [code/partial_forgetting_check.py](code/partial_forgetting_check.py): bundled forgetting grid, 68 checks.
- [code/compute_illustrations.py](code/compute_illustrations.py), [code/scaling_study.py](code/scaling_study.py) and
  [results/scaling.json](results/scaling.json): the Section 7 numbers and Figure 2.
- [code/discrete/](code/discrete/): exact Version 1-3 chains, the N = 150 illustration and the non-Poisson counterexample, with recorded outputs.
- All certificate scripts use explicit checks (no `assert`) and exit 1 on failure, also under `python -O`.
- Working derivations and sceptical checks: [development_notes/](development_notes/).

## Exclusions and limits

Neutral, idealised model; no empirical data and no fit to any population. Traits independent given the shared history
(no prerequisites, memory limits or culture-dependent contact). Threshold probabilities are not ordered. Identifying b
needs replicated populations and assumes whole-site resets. Complete loss is extremely rare in large populations.
Discrete-model results assume Poisson innovation; the full text of Kobayashi et al. (2018) was not obtained. Section 7
values are floating-point computations. Written proofs, not machine-verified; internal AI-assisted review only, no human
peer review. Prior-art screen bounded; it does not certify priority.

## Credit

Wakano, Ohtsuki, Kobayashi and Baake (2026): the model, the b-comparison by simulation and the b = 1 mixed-Poisson formula.
Wakano and Aoki (2024): the discrete versions and the open distributional question. Baumdicker and Pfaffelhuber (2014,
Remark 2.4): transferring several genes in one event preserves first moments but changes second moments. Standard tools:
Poisson marking, conditional Jensen, Blaszczyszyn and Yogeshwaran (2009) for convex order of Cox processes, Last-Penrose
chaos expansions, association inequalities. Licences: [LICENSE.md](LICENSE.md).
