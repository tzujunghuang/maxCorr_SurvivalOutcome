# maxCorr_SurvivalOutcome
This repository contains programming codes for the numerical studies in the article "Efficient Estimation of the Maximal Association between Multiple Predictors and a Survival Outcome" that is authored by Tzu-Jung Huang, Alex Luedtke and Ian W. McKeague (https://arxiv.org/abs/2112.10996) and accepted by Annals of Statistics. Below is the layout of this repository.

1. The directory *codes for stabilized one-step* in the folder *numerical studies* contains the main script *sim_maxCorrSurv.R* that generates the rejection decision results of Stabilized One-Step Estimator ('SOSE') and its competing methods (as mentioned in Competing Methods Section of the article) in one single Monte Carlo simulation run, 
using auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_data_acquisition.R* in the same directory. 
The file *sum_maxCorrSurv.R* helps to collect the results from 1000 Monte Carlo simulation runs that are implemented via parallel computing on 1000 nodes/CPUs,
and then gives type I error and power of all the methods in the cases of independent and dependent errors.

2. The directory *standalone simulation codes* in the folder *numerical studies* contains the main script *example_sim_maxCorrSurv.R* that demonstrates how to use our codes to generate the p-values and CIs of Stabilized One-Step Estimator ('SOSE') and its competing methods (as mentioned in Competing Methods Section of the article) in a standalone-computer-manageable scenario,
where our codes include auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_data_acquisition.R* in the same directory.
  
3. The directory *codes for 10-fold stabilized one-step* in the folder *numerical studies* contains the main script *sim_maxCorrSurv_rdo.R* that generates the rejection decision results of 10-fold Stabilized One-Step Estimator (introduced in Remark 4.3 in Stabilized One-Step Estimator Section of the article) in one single Monte Carlo simulation run, 
using auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_data_acquisition.R* in the same directory.
The file *sum_maxCorrSurv_rdo.R* helps to collect the results from 1000 Monte Carlo simulation runs that are implemented via parallel computing on 1000 nodes/CPUs,
and then gives type I error and power of all the methods in the cases of independent and dependent errors.

4. The directory *codes for power analysis* in the folder *numerical studies* contains the main script *sim_poweranalysis_maxCorrSurv.R* that gives the rejection decision results of Stabilized One-Step Estimator ('SOSE') and its competing methods in one single Monte Carlo simulation run, 
using auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_poweranalysis_data_acquisition.R* in the same directory.
The file *sum_poweranalysis_maxCorrSurv.R* helps to collect the results from 1000 Monte Carlo simulation runs that are implemented via parallel computing on 1000 nodes/CPUs, and then gives the power performance of all the methods in various settings in which signal has different types and strength
along with independent or dependent errors. 

5. The folder *real data application* contains the main script *VRC01_analysis.R* that applies Bonferroni Cox Method ('BONF_COX'), 
Bonferroni One-Step Estimator ('BONF_OSE') and Stabilized One-Step Estimator ('SOSE') with one random ordering of data (that is, r = 1) to 
the VRC01-mediated neutralization data, using auxiliary functions collected in the files *code_maxCorrSurv.R* and *real_data_acquisition.R* 
in the same directory. Further values of r depend on how many random orderings will be taken to generate the results for 'SOSE', 
say X-fold stabilized one-step estimator that uses r = X. In the real data application of this article, X = 200.
The file *VRC01_summary.R* helps to collect the results from 200 Monte Carlo simulation runs implemented on 200 nodes/CPUs, where in each of them
one random ordering of the original data is conducted.

