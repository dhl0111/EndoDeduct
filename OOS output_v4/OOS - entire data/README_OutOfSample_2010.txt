OUT-OF-SAMPLE EMPIRICAL VALIDATION: YEAR 2010
Training: policy years 2006-2009; validation: policy year 2010.
PRIMARY OOS EVALUATION: all 2010 policy-years retained.
All model and copula parameters are frozen at 2006-2009 estimates.
No selection, calibration or parameter fitting occurs on the 2010 outcomes.
Common complete-case evaluation n=1038 of 1038 policy-years (all nine headline models).
Primary loss: per-occurrence (min(Y, C_i)-D_i)_+ summed within policy-year.
C_i = Coverage_i * 1e+06 dollars; contract interpretation is assumed, not verified.
O/P = sum(observed losses) / sum(predicted pure premiums).
Poisson deviance, MAE, RMSE and relative portfolio bias use zero-loss policies.
RMSRE_positive excludes zero observed limited losses; N_positive is reported.
O/P = NA for observed=predicted=0; O/P = Inf for observed>0 and predicted=0.
Figure 1: policy-cluster bootstrap CIs condition on fixed fitted predictions.
Figure 1B: cases 1, 4c, 4d in each observation scheme.
Figure 2: test groups use only TRAINING-period XR1 risk cutpoints.
Tables 1-2 are fixed training copula specifications: NOT test-year AIC fits.
Appendix common copulas uses cached training candidates and FEM1 fixed.
If some test predictions are invalid, a common subset is used for all nine models.
Out-of-sample does not mean causally identified; frequency-severity factorization remains assumed.
