main_revised_v3.tex contains "Insurance ratemaking with endogenous deductibles" project tex file.

This project manuscript consists of three parts, Introduction & Methodology, Simulation and Empirical analysis with LGPIF dataset.

code folder inside simulation_update_code_tables provides the results, figures and tables in Section 3 simulation study. 

endogenous_deductible_pipeline.R file inside code folder gives the results in section 4.1 - 4.4.

Following works are remaining : 

1. Change of notation (N, Y, N^{*}, Y^{*})

 N^{*}, Y^{*} are chosen to define the frequency and severity of ground-up loss and N, Y are selected to describe the number and amount payment by an insurer. 
 However, some literature seems to utilize N, Y for the ground-up loss so change of notation needs to be implemented in main_revised_v3.tex.

 Then, aggregate insured loss (S) can be expressed as summation of Y^*. We aim to assess the performance of modeling with respect to S. 

 Aggregate insured loss can be found the product of conditional expectation of frequency and severity given deductible.

2. Summarize simulation section

The current simulation section spans approximately 10 pages, which seems excessively long. 
In particular, the section consists of four subsections, with figures presenting parameter estimation results and tables reporting the mean absolute errors (MAEs) of the frequency and severity models.

To streamline this section, I plan to shift the focus toward aggregate insured loss ((S)). 
Specifically, I intend to compare the MAE and mean relative error (MRE) of (S) across different models to investigate how ignoring deductible endogeneity affects the estimation of aggregate insured loss.

However, I do not intend to remove the existing results. Instead, I plan to move the detailed parameter estimation results, figures, tables, and plots to the Appendix. 
The Appendix will also provide comprehensive comparisons of the estimated (N), (Y), (N^), and (Y^) across models, using evaluation metrics such as MAE and mean squared relative error (MSRE).

Moreover, table 6 in current version (main_revised_v3.tex) illustrates the relative error of the estimated frequency and severity components of the pure premium by deductible level for only case 4. 
We need to include this framework when comparing MAE, MSRE of S with various deductible levels across all the four cases. These need to be included in the manuscript together.

Incorporating these remaining works above the new manuscript, R code need to be generated too.
