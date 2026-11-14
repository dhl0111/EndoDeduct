EMPIRICAL MANUSCRIPT REPORTING NOTES
Sample: In-sample descriptive evaluation
Main headline sample: 4152 / 4152 policy-years

Table 1: severity copula selection, one AIC ranking per deductible type
  and observation scenario (complete, censored, truncated).
Table 2: frequency copula selection, one ranking per deductible type
  and observation scenario (censored, truncated).
Do NOT compare AIC across different observed-data likelihoods or samples.
Only converged candidates enter delta-AIC and Akaike weights.

Table 3: all primary models are scored on the same valid policy-years.
Limited means (min(Y,C)-D)_+ with the policy's assumed ground-up C.
Mean Poisson deviance is a mean-oriented SCORE and does not require
the aggregate payments themselves to have a Poisson distribution.
Observed/predicted (O/P) near 1 is portfolio calibration, not proof
of individual predictive discrimination.

Figure 1: policy-level resampling within deductible groups. Fits are
FIXED; intervals are descriptive and do NOT propagate parameter
estimation uncertainty or account for in-sample model selection.
Bootstrap replicates: 199
Minimum independent policy clusters: 15

Figure 2: decile-like groups are based on the same XR1 premium and
are identical for every model; ties remain grouped together.

Appendix RMSRE_positive excludes policy-years with observed S=0.
RMSE remains available and is not the same as RMSRE.
Common-copula sensitivity holds FEM1 fixed, even for Scenario 3:
it isolates the severity copula and does NOT refit thinning models.
Unlimited premiums are supplementary; their finite mean is checked
conservatively in the original code.

IMPORTANT MODELING ASSUMPTIONS TO VALIDATE:
1. Coverage is a per-occurrence ground-up limit in dollars after scaling.
2. The claim data represent the observation mechanism in each likelihood.
3. The frequency-severity factorization is substantively justified.
4. Results are in-sample descriptive unless holdout_years is set.
