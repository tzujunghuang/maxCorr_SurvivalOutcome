# maxCorr_SurvivalOutcome
This repository containes programming codes for the numeric studeis in the article "Efficient Estimation of the Maximal Association between Multiple Predictors and a Survival Outcome."
Below is the layout of this repository.

1. The direcotry **codes for stabilized one-step** in the folder **numerical studies** contains the main script *sim_maxCorrSurv.R* that generates the results 
of Stabilized One-Step ('SOSE') and its competing methods (mentioned in Sec. 5 of the article) in one single Monte Carlo simulation run, 
using auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_data_acquiration.R* in the same directory. 
The file *sum_maxCorrSurv.R* helps to collect the results from 1000 Monte Carlo simulation runs that are implemented via parallel computing on 1000 nodes/CPUs.
  
2. The direcotry **codes for 10-fold stabilized one-step** in the folder **numerical studies** contains the main script *sim_maxCorrSurv_rdo.R* that generates the results 
of 10-fold Stabilized One-Step (introduced in Remark 4.3 in Sec. 4 of the article) in one single Monte Carlo simulation run, 
using auxiliary functions collected in the files *code_maxCorrSurv.R*, *data_management.R* and *sim_data_acquiration.R* in the same directory.
The file *sum_maxCorrSurv_rdo.R* helps to collect the results from 1000 Monte Carlo simulation runs that are implemented via parallel computing on 1000 nodes/CPUs.

3. The direcotry **power analysis** in the folder **numerical studies**

4. The folder **real data application**

[To be continued]
