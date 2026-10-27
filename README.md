# Discrete-deductible simulation study (Section 3)

Files

| file | purpose |
|---|---|
| `endogenous_deductible_pipeline.R` | estimation pipeline of the empirical analysis (per-unit likelihoods, `fit_stage`, `sandwich_se`, ...). Two changes relative to the version delivered with Section 4: (i) `pbvnorm()` now implements the Genz (2004) algorithm for |rho| >= 0.925 as well (vectorised; the previous row-by-row `mvtnorm` fallback made the optimiser extremely slow whenever a trial value of rho was large), (ii) the k-sum of eq. (10) evaluates the copula once on the grid k = n-1, ..., n+kmax and differences it (same values, half the cost), (iii) `minimize()` reports convergence = 0 when either L-BFGS-B or the nlminb polish converged. |
| `sim_functions.R` | design, data generation (`simulate_case`), stage constructors with the Pareto restriction (sigma = alpha1 = 1), starting values, `fit_models`, true/estimated conditional means given the deductible (`cond_sev`, `cond_freq`), `predict_all` (MAE and premium MAE), `run_replication` |
| `run_sim.R` | driver; one process per core, e.g. `for w in 1 2 3 4 5 6 7 8; do nohup Rscript run_sim.R 500 > log$w.txt 2>&1 & done` (all six cases) or `Rscript run_sim.R 500 1,2,3,4` (Cases 1-4 only). Replications are claimed with lock directories, so workers can be added or restarted at any time; each replication is saved as `results/case<c>_rep<r>.rds` |
| `summarize_sim.R` | `Rscript summarize_sim.R results tex_sim` builds the LaTeX tables of the Appendix, the MAE tables, the coverage table, the box-plot figures and `key_numbers.txt` from whatever replications are present |
| `update_manuscript_sim.py` | inserts Section 3 (`section3_new.tex` + `case_discussion.tex` + tables) and the Appendix tables into the manuscript |

Design (see the header of `sim_functions.R`): x = (1, N(0,1), Bernoulli(0.6)); ordinal deductible with levels 500/1000/2500/5000/10000,
P(D <= d_k | x) = logistic(kappa_k + x~'beta_D), kappa = (0.1, 1.0, 1.8, 2.6), beta_D = (-0.5, -0.5); N* ~ NB(r = 1, exp(x'beta_N)), beta_N = (-0.5, 0.5, 0.2);
Y* ~ Pareto(exp(x'beta_Y), alpha = 3), beta_Y = (8.5, 0.2, 0.5); Gaussian copulas (rho_N, rho_Y) = (0,0), (-0.5,0), (0,-0.5), (-0.6,-0.4) in Cases 1-4.
Cases 5-6 violate N* _||_ Y* | D (common latent propensity; plus a direct N*-Y* correlation of 0.3 in Case 6).

Seeds: replication r of case c uses seed 100000*c + r, so the results are reproducible and extendable.

Cost per replication on one core: about 35-55 s for Cases 1-3, 5, 6 and about 120 s for Case 4, which additionally computes
the Godambe standard errors of SDM3 and FDM2 for the coverage table (the FEM2 chain is computed in all cases). With 8 cores the
full design (6 cases x 500 replications) takes roughly 7-8 hours.

Output of `summarize_sim.R` (`tex_sim/`): `tab_app_D.tex`, `tab_app_case<c>_marginal.tex`, `tab_app_case<c>_joint.tex`,
`tab_mae_case<c>.tex`, `tab_coverage_case4.tex`, `coverage.csv`, `simcase<c>.pdf`, `key_numbers.txt`.
