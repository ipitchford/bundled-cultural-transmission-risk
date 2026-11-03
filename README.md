# Equal mean diversity, unequal risk

Version 0.1.0-candidate (unrefereed). DOI: https://doi.org/10.5281/zenodo.23933211 (concept DOI for all versions: https://doi.org/10.5281/zenodo.23933210). Release page: https://evidencepress.org/releases/bundled-cultural-transmission-risk/

Research package for the paper *Equal mean diversity, unequal risk: an exact comparison of partial and bundled cultural transmission* (Anonymous, 3 November 2026). Text and data CC0-1.0; code MIT.

**The result in plain words.** In a standard model of cultural trait accumulation (Wakano, Ohtsuki, Kobayashi & Baake 2026), hold fixed the rate at which each individual trait is passed on, and change only how many traits each learning encounter carries. The expected number of traits, and the expected number known to exactly k people, do not change. But every such count becomes more variable across possible histories (it increases in convex order), so the chance of losing everything and the expected shortfall below any target both rise with bundling. The same holds for bundled forgetting and, with Poisson innovation, for the discrete-generation versions 1-3 of Wakano & Aoki (2024). In large, well-connected populations the effect is large.

## Read first

- `manuscript.md` and `output/cultural_transmission_paper.pdf`: the paper (theorems, proofs, numerics, limits, references).
- `figures/`: schematic (Figure 1), variance scaling (Figure 2) and the original illustration.
- `code/`, `agent_math/`, `results/`: all computations and their machine-readable outputs.
- `development_notes/`: working derivations (`math_results.md`), the prior-art audit (`novelty_note.md`), sceptical checks, positioning note, the fuller discrete-generation section (`discrete_versions.md`) and the numerics notes. These are working records; the paper supersedes them where they differ.

## Rate conventions

There are N sites, reset (death) rate u per site and innovation rate mu per site. At copying probability b, the potential encounter rate per learner is beta/b, including a silent self-choice; a distinct ordered pair has rate beta/(N*b). Each donor trait is copied independently with probability b, so beta is the effective per-trait rate (Wakano et al.'s s). For N=2 the dimensionless edge rate is rho=beta/(2*u). The value b=0 denotes the limit of independent successful-contact clocks for different traits at fixed effective rate; resets remain shared. It does not denote no learning.

## Reproduction

Executed environment for this revision: Python 3.11.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, SymPy 1.14.0, single-threaded (`OPENBLAS_NUM_THREADS=1`). All analyses are offline. Every script inserts its own directory into `sys.path`, so it runs under isolated mode. From this directory:

```bash
PY="uv run --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 --with sympy python3 -I -B"
$PY code/run_checks.py                     # ~40 s (includes the 400,000-cycle simulation)
$PY code/run_checks.py --reuse-forward results/verification.json --out DIR   # exact checks only, ~2 s
$PY agent_math/pair_ancestry.py            # ~1 s
$PY code/forward_closure_check.py          # ~1 min
$PY code/partial_forgetting_check.py       # ~2 s
$PY code/exact_extras.py                   # tens of minutes (two truncated chains per b at cmax 60 and 45)
$PY code/compute_illustrations.py --large  # ~15 s
$PY code/scaling_study.py                  # ~26 min single core; writes results/scaling.json
$PY code/plot_scaling.py                   # reads results/scaling.json
$PY code/plot_transmission.py [--paper]
$PY code/plot_schematic.py
# discrete-generation versions (Section 6, Appendix E):
$PY code/discrete/exact_versions.py 3 2 4/5 1/2
$PY code/discrete/ya_chain.py 150 2 1/2 1 10
$PY code/discrete/bernoulli_counterexample.py
$PY code/build_pdfs.py                     # needs pandoc and pdflatex on PATH
```

Each script that writes results accepts `--out DIR` (default `results/`, or `figures/`, `output/`, `agent_math/` as appropriate). Files in the default directory are refreshed; an existing file elsewhere is never overwritten. Every certificate script uses explicit checks rather than `assert`, sets its status string from the recorded results and exits with status 1 if any check fails, including under `python -O`. All scripts were run once normally and once with `-O`; the simulation was run only once (the `-O` run re-checked the saved simulations with `--reuse-forward`).

## Verification details

- `code/run_checks.py`: rational N=2 pair system, closed form and inverse; 24 labelled-versus-lumped backward comparisons (N=2,3,4, including b=0); 24 no-learning boundary checks; 15 b=1 single-genealogy checks; four forward simulations of 100,000 regenerative cycles each (realised |z| at most 0.71 for means and 1.31 for variances). 87 checks.
- `agent_math/pair_ancestry.py`: independently written backward ancestral implementation, N=2 formula, lumped/labelled agreement, a heterogeneous directed network. 50 checks.
- `code/forward_closure_check.py`: an independent FORWARD moment closure on carrier-set counts (`code/forward_closure.py`), exact rationals, compared by exact equality with the bundle's lumped backward generator and with an independent labelled backward generator: homogeneous N=3 and N=4 (b = 0, 1/4, 1/2, 1), a heterogeneous directed N=3 network (b = 0, 1/3, 2/3, 1) and a directed chain 0->1->2 with innovation only at 0 and sample {2}. Also checks constant means and nondecreasing variances of C_A, every popularity count P_k and every count of traits held by at least k sites. 146 checks.
- `code/partial_forgetting_check.py`: the same closure with loss events at rate u_i/d deleting each trait with probability d; exact grid over (b,d) in {0,1/4,1/2,1}^2 for N=2,3; Poisson law at b=d=0; symbolic N=2 closed form reducing to the paper's formula at d=1. 68 checks.
- `code/exact_extras.py`: N=2 stationary law P(C<=q) from a truncated forward chain with a regenerative truncation bound; exact b=1 empty-state probabilities by Feynman-Kac; exact Var P_1, Var P_2 at N=2.
- `code/scaling_study.py`: V_N(b) for beta/u in {0.7,1,1.3,2}, N up to 150, by a cancellation-free covariance form of the pair system, Jacobi-BiCGSTAB inner solves and iterative refinement with double-double residuals from exact rational rates. Cross-checked against direct sparse LU (N<=50), the original Section 6.1 formulation (N<=30) and the single-genealogy recursion at b=1; reports residuals and a floating-point a-posteriori error bound per point.
- `code/compute_illustrations.py`: N=50 figure points (direct) and N=50/100 b=1/2 Jacobi-BiCGSTAB solves. 12 checks.


Computational precision and Monte Carlo error are distinct from empirical validity. The model uses no empirical data.

## Scope and provenance

The contribution is the exact matched-rate ordering (for every carrier-pattern count), its strictness characterisation, the bundled-forgetting and discrete-generation extensions, and their variance, identification and regenerative consequences in these models. Poisson marking, stochastic-order principles, Poisson noise expansions, association inequalities, renewal identities and paired-genealogy methods are credited to prior work. The full text of Kobayashi et al. (2018) was not obtained; the discrete model follows Wakano & Aoki (2024). The prior-art screen does not certify priority. Checks were internal and AI-assisted; the paper has not had human peer review.

No retrieved article PDFs, third-party code or web-page copies are redistributed.

## How to cite

Anonymous (2026). *Equal mean diversity, unequal risk: an exact comparison of partial and bundled cultural transmission* (Version 0.1.0-candidate) [Unrefereed research candidate]. Evidence Press. https://doi.org/10.5281/zenodo.23933211

Machine-readable citation: `CITATION.cff`. Claims map for automated readers: `AI_INDEX.md`. Licences: `LICENSE.md`.
