# Manuscript Revision — Version 4 (v4)

**Last updated:** November 14, 2026  
**Status:** Major revisions completed; OOS validation and final manuscript sections pending.

## 1. Summary of Revisions

The manuscript has undergone substantial revisions over the past month. The major changes are summarized below:

1. **Variable definitions (Section 2.1, Table 1)**
   - Revised the definitions and notation of \(N^*\), \(Y^*\), \(N\), and \(Y\).

2. **Terminology**
   - Replaced *reporting scenario* with *observation scenario* throughout the manuscript.
   - The revised terminology is intended to be more neutral, as "reporting" may place excessive emphasis on policyholders' reporting behavior.

3. **Standard error estimation (Section 3, Table 6)**
   - Incorporated sandwich estimators to account for estimation uncertainty.

4. **Manuscript organization**
   - Condensed and reorganized the simulation and empirical analysis sections to improve clarity and balance.

5. **Ratemaking comparison (Section 2.5)**
   - Added comparisons of pure premium estimates between exogenous and endogenous models.

### Current Manuscript Status

All sections except the Introduction and Conclusion have been revised.

- **Blue text:** Newly added material or content that may require further reorganization.
- **Red text:** Incomplete sections, primarily the out-of-sample (OOS) validation.

## 2. Methodological Considerations

### 2.1 Conditional Independence Assumption

The current framework remains subject to limitations arising from the conditional independence assumption.

Nevertheless, the manuscript focuses on the practical contribution of comparing ratemaking outcomes through conditional expectations under different observation scenarios.

The revised manuscript aims to present this contribution more concisely without overstating the implications of the dependence modeling framework.

### 2.2 Conditional Expectation and Pure Premium Calculation (Appendix E)

The mathematical formulations for conditional expectation calculations have been added to Appendix E, with AI assistance in drafting the mathematical exposition.

Further revisions will focus on retaining only the essential mathematical details.

In the empirical analysis, the GB2 distribution and certain copula families impose conditions on the existence of finite expectations. However, by restricting the integration range to the policy coverage limit, finite conditional expectations can be obtained for the fitted models.

### 2.3 Simulation Evaluation

The simulation study compares:

- **Theoretical pure premium:** Calculated using the true data-generating parameters.
- **Fitted pure premium:** Calculated using the estimated model parameters.

This comparison evaluates pure premium estimation accuracy under incomplete observation scenarios, independently of the variability caused by extreme claims in realized aggregate losses.

## 3. Empirical Analysis: LGPIF Data

### 3.1 Deductible Endogeneity

A definitive methodology for identifying deductible endogeneity has not yet been established within the current framework.

However, several empirical findings suggest potential departures from exogeneity, including:

- The coverage coefficient estimated in the severity models (SEM).
- The estimated slope of the GB2 location parameter \(\mu\) with respect to coverage.

Based on these findings, the analysis proceeds with the SDM and FDM models.

Comparisons of regression coefficients and chi-square tests are retained to examine differences between exogenous and endogenous specifications.

### 3.2 In-sample Ratemaking (Sections 4.1–4.4)

The continuous models generally improve several measures of predictive performance relative to the marginal models, although the improvements are not consistent across all metrics.

**Full and censored observation scenarios:**

The continuous models show improvements in several measures for aggregate insured losses:

- Mean Absolute Error (MAE)
- Root Mean Squared Error (RMSE)
- Observed-to-Predicted Ratio
- Mean Poisson Deviance

**Truncated observation scenario:**

The continuous models demonstrate performance broadly comparable to that of the marginal models.

These findings suggest potential benefits from incorporating dependence structures into ratemaking, particularly under full and censored observations.

## 4. Out-of-sample Validation (Section 4.5)

**Status: In progress — Tables 15–16**

The OOS analysis raises several questions that require further discussion.

### 4.1 Treatment of Extreme Claims

One exceptionally large claim occurred in the 2010 policy year, accounting for approximately **18% of total ground-up losses during 2006–2010**.

Prediction performance has been evaluated both including and excluding this observation.

**Open question:** Could excluding this extreme claim introduce bias into the evaluation of predictive performance, potentially overstating the accuracy of certain models?

### 4.2 Interpretation of Predictive Performance

Similar to the in-sample analysis, the continuous models outperform the marginal models under some evaluation measures, particularly in the full and censored observation scenarios.

However, improvements are not consistent across all metrics.

**Open question:** Are these improvements sufficiently meaningful to support a conclusion that continuous dependence models offer better predictive performance?

### 4.3 Actuarial and Practical Implications

The XR model itself demonstrates reasonably good OOS predictive performance.

Additional comparisons of prediction metrics across deductible levels, similar to Table 5 in the simulation study, have not yet revealed clear practical implications for ratemaking.

**Open question:** How should the results be interpreted when more complex dependence models offer only limited or inconsistent improvements over simpler marginal models?

## 5. Remaining Tasks

- [x] Revise variable definitions and observation scenario terminology.
- [x] Incorporate sandwich standard errors.
- [x] Reorganize simulation and empirical analysis sections.
- [x] Add exogenous/endogenous ratemaking comparisons.
- [x] Add conditional expectation formulations to Appendix E.
- [x] Complete the main in-sample ratemaking comparisons.
- [ ] Finalize the OOS validation and interpretation.
- [ ] Determine how to present results involving extreme claims.
- [ ] Clarify the practical implications of deductible-dependent ratemaking.
- [ ] Further condense Appendix E.
- [ ] Revise the Introduction and Conclusion.
- [ ] Complete the final manuscript review.

## 6. Discussion Priorities

The main issues requiring discussion are:

1. Whether the OOS analysis should present results both including and excluding the extreme claim, and how these results should be interpreted.
2. Whether the observed improvements in continuous models provide sufficient evidence of practical predictive benefits.
3. How to position the contribution of dependence-based ratemaking when the simpler XR model already performs reasonably well.
4. How to acknowledge the conditional independence limitation while maintaining a clear methodological contribution.

---

**Version:** v4  
**Updated:** November 14, 2026
