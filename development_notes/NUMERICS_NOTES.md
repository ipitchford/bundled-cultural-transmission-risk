# Numerics notes for the revision (code/numerics engineer, 3 Nov 2026)

Everything below was produced in this working copy with Python 3.11.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 and SymPy 1.14.0, single-threaded. Each script was run once normally (into `results/` / `figures/`) and once under `python -O` (into a scratch directory). All certificate scripts exited 0 in both modes. A deliberate failing check under `-O` exited 1. Parameters: u = 1 unless stated; the edge rate is beta/N (silent self-choice).

## 1. Certification layer (reviewer W1)

- There are no `assert` statements left in `code/*.py` or `agent_math/*.py`. Checks are explicit (`code/certify.py`, plus a self-contained class in `agent_math/pair_ancestry.py`). Each script sets its status from the recorded results and exits with status 1 on any failure, including under `-O`.
- Every script runs as `python3 -I -B path/to/script.py`. Scripts that write results take `--out DIR`. A script refuses to overwrite an existing file outside its default directory.
- `run_checks.py --reuse-forward results/verification.json` re-checks the saved simulations without rerunning them. The 400,000-cycle simulation was run once. It reproduced the earlier seeded values exactly, in 36 s.
- Realised simulation z-scores, for the §8.2 text (reviewer asked for these instead of "within six SE"):

| b | mean z | variance z | P̂(C=0) |
|--:|--:|--:|--:|
| 0.1 | 0.64 | 1.05 | 0.3310 |
| 0.25 | 0.71 | 1.27 | 0.3365 |
| 0.5 | 0.38 | 1.31 | 0.3479 |
| 1 | −0.64 | −1.04 | 0.3634 |

  Maximum |z| is 0.71 for the mean and 1.31 for the variance. The Monte Carlo SE of P̂ is about 0.0011.

## 2. Independent forward closure (reviewer W2): `code/forward_closure_check.py`, 146 checks passed

The forward moment closure is on carrier-set counts n_S, in exact rationals, with the b→0 limit taken analytically. It equals the backward pair system **exactly**, with rational equality, in every case. For homogeneous cases it was compared both with the bundle's lumped generator (solved in Fractions) and with an independent labelled backward generator. It also reproduces the methodology reviewer's values: N=3 Var = 401845979721590906/304364967034359375 at b=1/2 and 1854829/1265625 at b=1; N=4 Var = 1361/64 at b=1.

In every case, **means are constant in b, and the variances of C_A, every P_k and every G_k (traits held by ≥ k sites) are strictly increasing in b.**

| Case (b grid) | E C (whole population) | Var C_V across the grid |
|:--|--:|:--|
| Homogeneous N=3, β=13/10, μ=1/5 (0, ¼, ½, 1) | 1052/1125 | 1.23947, 1.27497, 1.32028, 1.46554 |
| Homogeneous N=4, β=2, μ=1/2 (0, ¼, ½, 1) | 39/8 | 11.3294, 12.4570, 14.0648, 21.2656 |
| Directed heterogeneous N=3 (reviewer's network) (0, ⅓, ⅔, 1) | 10145/3144 | 7.60665, 8.28062, 9.39598, 11.5684 |
| Same network, A={1} | 3205/3144 | 2.3577, 2.9961, 3.9463, 5.5886 |
| **Chain 0→1→2, μ only at 0, a=1 per edge, A={2}** (0, ¼, ½, ¾, 1) | 43/99 | 0.719205, 0.808921, 0.92956, 1.09347, 1.31946 |

The chain row shows that strictness holds beyond Proposition 1's direct-edge condition. Sample {2}'s only in-neighbour has μ=0, yet the variance increases strictly. This supports the reviewer's conjectured path characterisation (a site with μ_j>0 and a directed path into A).

Popularity counts (homogeneous N=4, b = 0, ¼, ½, 1):
- Var P_1: 2.767, 3.364, 4.303, 9.232
- Var P_2: 2.397, 2.818, 3.631, 9.177
- Var P_3: 1.920, 2.285, 2.937, 7.259
- Var P_4: 0.818, 1.014, 1.332, 3.099

E P_k = 2, 3/2, 1, 3/8, matching the §7 formula.

**Suggested text:** the proof of Theorems 1 and 2 applies unchanged to P_k, G_k and any count of traits whose carrier configuration lies in a fixed set. Each innovation's configuration is a function of its own marked path. So each such count is conditionally Poisson with mean E[·|G_b,D] and increases in convex order. The numerics confirm this; the lead may want to state it as a remark.

## 3. Partial forgetting (`code/partial_forgetting_check.py`, 68 checks passed)

The extended model has loss events at rate u_i/d, each deleting each trait at the site with probability d. Over the exact grid (b,d) ∈ {0 (limit), ¼, ½, 1}² for N=2 (β=2, μ=½) and N=3 (β=13/10, μ=1/5):
- the means of C_V and of every P_k and G_k are invariant;
- every variance is strictly increasing in b at each d, and in d at each b;
- at b=d=0 every count is exactly Poisson (Var = mean).

N=2, Var C_V (rows b, columns d = 0, ¼, ½, 1):

| b \ d | 0 | ¼ | ½ | 1 |
|--:|--:|--:|--:|--:|
| 0 | 1.500000 | 1.661104 | 1.877451 | 2.650000 |
| ¼ | 1.508197 | 1.675743 | 1.903846 | 2.750000 |
| ½ | 1.517241 | 1.692186 | 1.934211 | 2.875000 |
| 1 | 1.538462 | 1.732014 | 2.011062 | 3.250000 |

The d=1 column is the paper's 53/20, 11/4, 23/8, 13/4. **Partial forgetting changes the variance much more than partial copying does.** At d=0 the b-effect almost vanishes: 1.5 → 1.538.

N=3, Var C_V (same layout):

| b \ d | 0 | ¼ | ½ | 1 |
|--:|--:|--:|--:|--:|
| 0 | 0.935111 | 0.978007 | 1.035384 | 1.239466 |
| ¼ | 0.938739 | 0.983938 | 1.045429 | 1.274971 |
| ½ | 0.942698 | 0.990547 | 1.056947 | 1.320277 |
| 1 | 0.951847 | 1.006404 | 1.086151 | 1.465544 |

**N=2 exact closed form**, from sympy and checked against the grid, with ρ = β/(2u) and η = μ/u:

$$\operatorname{Var}C_V=\eta(2+\rho)+\eta^2\,\frac{d\,A(d,\rho)-2b\rho(\rho+1)\,[d(\rho+1)-1]}{(2-d)\,B(d,\rho)-2b\rho(\rho+2-d)},$$

where A = (ρ²+2ρ+2)d² + (ρ⁴+ρ³−7ρ²−16ρ−10)d + 4(ρ+1)(2ρ²+4ρ+3) and B = d² + (ρ²−3ρ−5)d + 2(ρ+1)(ρ+3).

At d=1 the numerator and denominator share the factor (ρ+1), and the formula reduces exactly to the paper's. At d=0 the excess is bρ(ρ+1)/[2(ρ+1)(ρ+3) − bρ(ρ+2)]. The expressions for ∂/∂b and ∂/∂d are in `results/partial_forgetting.json`.

## 4. Scaling study (the main new numerics): `code/scaling_study.py`, `results/scaling.json`, `figures/variance_scaling.png`

**Method change (important).** The §6.1 formulation Var = μℓ + μ²(h − ℓ²) cancels catastrophically in supercritical cases. At β=2, N=100, b=0, the old direct or iterative h-solve gave 7.83e11; the correct value is 8.18e11. I therefore solve an equivalent cancellation-free system for the replica covariance c = F − ℓ⊗ℓ:

$$-Q_b^{(2)}c=\Gamma_b,\qquad \Gamma_b(i,j,k)=u\,k\,\Delta_{n_1}\Delta_{n_2}+\tfrac{\beta}{N}\,b\,(N-i-j-k)\,k\,\Delta_{n_1+1}\Delta_{n_2+1},$$

with n1=i+k, n2=j+k, Δ_n = ℓ_n − ℓ_{n−1} (Δ_{N+1}=0), and Var C_V = μℓ_N + μ² c(0,0,N). Γ is the carré-du-champ of the joint jumps: shared resets plus shared successful arrows. It is ≥ 0, so the computation involves no subtraction. **Suggested remark for §6:** the excess variance equals μ² × the expected integrated Γ along the pair chain. At b=0 only shared resets contribute.

The solver uses Jacobi-preconditioned BiCGSTAB inner solves, warm-started along b. They sit inside iterative refinement whose residuals are computed in double-double from exact rational rates. Final relative residuals are 1e-27 to 1e-30.

Cross-checks, all passed:
- against direct sparse LU for every N ≤ 50 point (max relative difference 8.6e-11 at β=2, N=50; ≤ 4e-14 elsewhere);
- against the original h-formulation for N ≤ 30;
- against the single-genealogy recursion at b=1 for every (β,N) (≤ 6e-15 relative).

An a-posteriori bound T_max‖residual‖, evaluated in floating point, is ≤ 1e-26 relative for every reported point except β=2, N=100, where it is 3e-16.

**Not computed: β=2, N=150.** The pair chain's exit time is about 1.6e12, and the double-precision inner solves diverge. The script skips this case and says so. The full sweep took 26 min on one core. The slowest single solves are about 1 min, at N=150.

Grid: u=1, μ=0.1, b ∈ {0, .1, .25, .5, .75, .9, .95, .99, .999, 1}, plus 13/14 for β=1.3. Full tables are in the JSON. Key columns follow; "slope" is (V(1)−V(.99))/.01, and "share > 0.9" is (V(1)−V(.9))/(V(1)−V(0)).

| β/u | N | E C | V(0) | V(½) | V(1) | V(0)/V(1) | V(½)/V(1) | V(.9)/V(1) | V(.99)/V(1) | slope | share > 0.9 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 0.7 | 10 | 1.4958 | 1.7325 | 1.8161 | 1.9801 | 0.875 | 0.917 | 0.977 | 0.998 | 0.49 | 0.18 |
| 0.7 | 50 | 8.2734 | 9.7817 | 10.535 | 12.495 | 0.783 | 0.843 | 0.950 | 0.994 | 7.3 | 0.23 |
| 0.7 | 100 | 16.847 | 20.005 | 21.683 | 26.366 | 0.759 | 0.822 | 0.941 | 0.993 | 18.5 | 0.25 |
| 0.7 | 150 | 25.436 | 30.254 | 32.873 | 40.393 | 0.749 | 0.814 | 0.937 | 0.993 | 30.5 | 0.25 |
| 1.0 | 10 | 1.9130 | 2.3269 | 2.5318 | 3.0671 | 0.759 | 0.826 | 0.944 | 0.994 | 2.0 | 0.23 |
| 1.0 | 50 | 13.246 | 18.231 | 22.414 | 43.188 | 0.422 | 0.519 | 0.780 | 0.968 | 138 | 0.38 |
| 1.0 | 100 | 29.790 | 44.189 | 58.311 | 150.61 | 0.293 | 0.387 | 0.678 | 0.947 | 800 | 0.46 |
| 1.0 | 150 | 47.614 | 74.332 | 102.70 | 320.50 | 0.232 | 0.320 | 0.616 | 0.931 | 2211 | 0.50 |
| 1.3 | 10 | 2.5933 | 3.4151 | 3.9086 | 5.5902 | 0.611 | 0.699 | 0.889 | 0.986 | 7.8 | 0.29 |
| 1.3 | 20 | 7.4480 | 11.553 | 14.822 | 31.969 | 0.361 | 0.464 | 0.748 | 0.963 | 118 | 0.39 |
| 1.3 | 30 | 14.656 | 27.106 | 38.392 | 116.85 | 0.232 | 0.329 | 0.637 | 0.939 | 712 | 0.47 |
| 1.3 | 50 | 39.141 | 107.34 | 177.97 | 914.95 | 0.117 | 0.195 | 0.482 | 0.892 | 9924 | 0.59 |
| 1.3 | 75 | 103.91 | 485.76 | 918.13 | 7649.1 | 0.064 | 0.120 | 0.358 | 0.832 | 1.28e5 | 0.69 |
| 1.3 | 100 | 251.18 | 1982.3 | 4074.2 | 51156.7 | 0.039 | 0.080 | 0.271 | 0.773 | 1.16e6 | 0.76 |
| 1.3 | 150 | 1364.2 | 26772 | 61655 | 1.7429e6 | 0.015 | 0.035 | 0.154 | 0.659 | 5.95e7 | 0.86 |
| 2.0 | 10 | 6.6071 | 12.132 | 16.196 | 39.776 | 0.305 | 0.407 | 0.717 | 0.958 | 165 | 0.41 |
| 2.0 | 20 | 62.421 | 281.79 | 505.55 | 3646.2 | 0.077 | 0.139 | 0.433 | 0.885 | 4.2e4 | 0.62 |
| 2.0 | 30 | 497.47 | 5692.4 | 12898 | 2.4339e5 | 0.023 | 0.053 | 0.269 | 0.811 | 4.6e6 | 0.75 |
| 2.0 | 50 | 29131 | 1.5988e6 | 6.0422e6 | 8.4813e8 | 0.0019 | 0.0071 | 0.110 | 0.682 | 2.7e10 | 0.89 |
| 2.0 | 75 | 4.3775e6 | 1.2450e9 | 1.0175e10 | 1.9163e13 | 6.5e-5 | 5.3e-4 | 0.039 | 0.560 | 8.4e14 | 0.96 |
| 2.0 | 100 | 6.2678e8 | 8.1786e11 | 1.5481e13 | 3.9286e17 | 2.1e-6 | 3.9e-5 | 0.015 | 0.466 | 2.1e19 | 0.985 |

β=1.0 and β=0.7 at N=20, 30 and 75 are in the JSON.

**Wakano et al. (2026) Fig. 16 settings (β=1.3, N=100):** V(0.99) = 39534.0 (77.3% of V(1)), V(13/14) = 17179.19 (33.6% of V(1)), V(1) = 51156.66, E C = 251.18. At b = 13/14, therefore, the variance is only a third of the whole-repertoire value, despite a 93% copying probability.

**Sensitivity V′(b)/V(b)**, from the exact derivative solve, β=1.3: N=100 gives 1.25 at b=0, 1.79 at b=½ and 29.6 at b=1; N=50 gives 0.81, 1.33 and 12.4. Identification of b (Corollary 3) is therefore sharp near b=1 and weak at small b.

**Findings to state (W3):**
1. Regime dependence. Subcritically (β/u=0.7) V(0)/V(1) stays at 0.75–0.88 and saturates with N. At β/u=1 it falls slowly (0.23 at N=150). Supercritically it collapses with N: 0.015 at β=1.3, N=150, and 2e-6 at β=2, N=100.
2. Boundary layer at b=1. The share of the total increase above b=0.9 rises with N to 0.86 (β=1.3, N=150) and 0.985 (β=2, N=100). Endpoint slopes grow roughly like V(1) times a factor that increases with N; V′(1)/V(1) is 30 at β=1.3, N=100.
3. Every computed curve is convex and increasing, as §4.2 predicts. Nothing here is a large-N theorem.

## 5. Small exact extras (`code/exact_extras.py`, 20 checks passed)

**(a) N=2, u=1, β=2 (ρ=1), μ=½ (η=½).** The stationary P(C ≤ q) comes from the forward (X,Y,Z) chain truncated at X+Y+Z ≤ 45. The truncation bound comes from regenerative coupling at the empty state (method in the JSON): ≤ 5.7e-13 (b=¼), 1.3e-12 (b=½) and 8.5e-10 (b=1). Values changed by ≤ 1.1e-10 against truncation at 35. The truncated means and variances reproduce 3/2 and 11/4, 23/8, 13/4, and P(C=0) at b=1 = 4/11.

| q | b=¼ | b=½ | b=1 | pattern |
|--:|--:|--:|--:|:--|
| 0 | 0.336751 | 0.346256 | 0.363636 (=4/11) | increasing |
| 1 | 0.610286 | 0.614643 | 0.628099 | increasing |
| **2** | **0.785021** | **0.784029** | **0.787378** | **not monotone** |
| 3 | 0.885961 | 0.882856 | 0.879175 | decreasing |
| 4 | 0.941089 | 0.937938 | 0.931456 | decreasing |
| 5 | 0.970144 | 0.967698 | 0.961134 | decreasing |
| 6 | 0.985086 | 0.983418 | 0.977965 | decreasing |

So P(C ≤ q) rises with b for q ≤ 1, falls for q ≥ 3, and is not monotone at q = 2. This is a concrete instance of §4.1's warning. The exact P(C=0) values at b=¼ and ½ (0.336751, 0.346256) agree with the simulated 0.3365 and 0.3479 within 1 SE.

**(b) b=1, exact P(C_V=0)** by Feynman–Kac on the full-population ancestral birth–death chain, in rationals (β=1.3, μ=0.1, u=1). Corollary 2 gives the loss frequency r = Nμπ₀ and the mean nonempty episode (1−π₀)/(Nμπ₀).

| N | P(C_V=0) | no-learning bound (10/11)^N | ratio | mean nonempty episode (1/u units) |
|--:|--:|--:|--:|--:|
| 2 | 40/51 = 0.784314 | 0.826446 | 0.949 | 1.375 |
| 10 | 0.161111 | 0.385543 | 0.418 | 5.21 |
| 20 | 0.0181452 | 0.148644 | 0.122 | 27.1 |
| 50 | 2.29096e-5 | 8.51855e-3 | 2.7e-3 | 8.73e3 |
| 100 | 3.24962e-10 | 7.25657e-5 | 4.5e-6 | 3.08e8 |

The bound prod u/(u+μ) is a valid upper bound for every b (couple the inventors' own innovations), but it is very loose for large N. Exact fractions for N ≤ 20 are in the JSON. The N=2, β=2, μ=½ check returns 4/11.

**(c) N=2, u=1, β=2, μ=½, exact.** E P_1 = 1 and E P_2 = ½ at every b.

| b | Var P_1 | Var P_2 | Var C |
|--:|--:|--:|--:|
| 0 | 29/20 | 9/10 | 53/20 |
| ¼ | 19/12 | 1 | 11/4 |
| ½ | 7/4 | 9/8 | 23/8 |
| 1 | 9/4 | 3/2 | 13/4 |

## 6. Figures

- `figures/variance_scaling.png`: V_N(b)/V_N(1) against b, with panels for β/u = 0.7, 1, 1.3, 2 and viridis-ordered curves for N. Caption suggestion: "Exact solutions of the pair-replica equations (u=1, μ=0.1); points at the computed b values. β/u=2, N=150 omitted (not numerically certifiable)."
- `figures/schematic.png`: two panels. Bundled copying (b=1) shows 3 encounters with both traits passing together. Partial copying (b=¼) shows 12 encounters with the traits mostly passing separately. A cross marks a reset. Okabe–Ito blue and vermillion; matplotlib only.

## 7. Failures and surprises

- **The h-formulation of §6.1 is numerically unsafe in strongly supercritical cases.** It loses about 6 digits at β=2, N=100, and the old solver gave a variance about 4% too low at b=0. Every number the manuscript currently quotes (β=1.3, N ≤ 100) agrees with the new method to the printed digits. The new method is now the reference.
- β=2, N=150 could not be computed in double-precision Krylov, even with refinement. It is reported as not computed.
- Under `--large`, `compute_illustrations.py` now uses Jacobi-BiCGSTAB in place of ILU-GMRES. N=100, b=½ gives 4074.17250296 (relative residual 1.1e-11; the old ILU path was slow).
- `build_pdfs.py` could not run because xelatex is not on PATH. `build_bundle.py` fails in this working copy because `reviewer_report.md` is absent. Both are packaging steps for the publication run.
- `compute_illustrations.py` rewrites `results/large_N_moments.json` with the six N=50 points. The exploratory N=75 values it held are superseded by `scaling.json`.
- The first attempt of `exact_extras.py` at cmax=40 failed its own check at b=1: the truncation bound was 1.3e-8 and the mean was off by 2e-10. The final version uses cmax=45 with tolerances matched to the bound.
