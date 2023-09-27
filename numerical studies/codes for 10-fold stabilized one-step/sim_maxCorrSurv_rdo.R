##### This script is for one single MC run for 10-fold stabilized one-step estimator ('SOSE'), say r = 1.
##### To reproduce our numeric results, it is required that r = 1,...,1000 and 
##### it is highly suggested that parallel computing techniques should be used.


##### Configuration
# Read in source codes that should be in the same working directory as is this script	
source('code_maxCorrSurv.R')
source('sim_data_acquisition.R')
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '10%'
# Number of random orderings per independent data set
num_rep = 10
# Methods to run
meths_vec = c('SOSE')
# List of (n,p) values to run simulation
np_list = list(c(500, 1e6), c(500, 1e5), c(500, 1e4), c(500, 1000), c(500, 100))
# List of correlation values between predictors
rho_val = 0.75
# Tau
quar = 0.9
# Significance level
alpha_val = 0.05
# elln
elln_part = c(2); dim_elln = length(elln_part)
# List of data generating distributions to use when running simulations
model_list = c('N.DE','A2.DE','A1.DE','N.IE','A2.IE','A1.IE')


##### In this script, we set r = 1 because it's for a single MC run.
r = 1
# slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# r <- as.numeric(slurm_arrayid)
# set.seed(r)

# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr') #'foreach', 'doParallel'
# new_packages = package_list[!(package_list %in% installed.packages()[,"Package"])]
# if (length(new_packages)) { install.packages(new_packages) } 
for (package in package_list) {
    require(package, character.only = TRUE); library(package, character.only = TRUE) }

# n_cores = detectCores() - 2 # 3 cores detected on local machine; 6 on each node of clusters
# cl = makeCluster(n_cores)
# registerDoParallel(cl) 

# Number of digits used here
num_digits = 7; options(digits = num_digits)


sim = data.frame(n = NA, p = NA, rho = NA, model = NA, method = NA, 
                 quar = NA, elln = NA, 
                 Rej1 = NA, Rej2 = NA, Rej3 = NA, Rej4 = NA,
                 test_stat = NA, mean_pval = NA, med_pval = NA, min_pval = NA)

for (model in model_list) {
	for (np in np_list) {
	  
	  n = np[1]; p = np[2]; rho = rho_val
	  dat = Simulate_data(n, p, model, censoring_rate = cr, rho = rho, 
	                      censoring_dist_par = c(0.07, 0.15, 0.28, 0.43))
	  obj0 <- NumericalStudy$new(input_data = dat)
	  print(c(n, p, rho, model, cr))
	  
	  out = list()
	  for (meth in meths_vec) {
	    if (meth == 'SOSE') {
        
	      print(meth)
	      SOSE_est = do.call(rbind, lapply(seq(num_rep), function(rep_idx){
	        print(rep_idx)
	        inds = sample(1:n)
	        dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
	                    U = dat$U[1:n,1:p][inds,])
	        obj1 <- NumericalStudy$new(input_data = dat0)
	        
	        do.call(rbind, lapply(1:length(elln_part), function(d_idx){
	          d = elln_part[d_idx]; elln = ceiling(n/d)
	          chunk_size = ceiling((n - elln)/10)
	          SOSE_est = obj1$Stab_onestep_est(all_obs = 1:n, chunk_size, elln, 
	                                           est_index = est_index_val, 
	                                           alpha = alpha_val, 
	                                           num_top = 1, quar_trunc = quar)
	          t(c( SOSE_est$p_val, SOSE_est$est, SOSE_est$se, elln )) } ) )
	      } ) )
	      SOSE_est = as.data.frame(SOSE_est); 
	      names(SOSE_est) = c('p_val', 'est', 'se', 'elln')
	      out$SOSE = SOSE_est
		        
		 } else stop('Invalid method.') }
	  gc()
	  
	  
	  sim = rbind(sim, do.call(rbind,lapply(meths_vec, function(meth) {
	    if ( meth == 'SOSE' ) {
	      
	      result = out[[meth]] %>% group_by(elln) %>% 
	        summarise(test_stat = abs(mean(est)/mean(se)), 
	                  mean_pval = mean(p_val), 
	                  med_pval = quantile(p_val, prob = 0.5), 
	                  min_pval = min(p_val)) %>%
	        mutate( n = n, p = p, rho = rho, model = model, method = meth, quar = quar,
	                Rej1 = 1*(mean_pval <= alpha_val), 
	                Rej2 = 1*(med_pval <= alpha_val),
	                Rej3 = 1*(min_pval <= alpha_val/num_rep), 
	                Rej4 = 1*(pnorm(test_stat, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE) < alpha_val/2) )
	      
	      result = as.data.frame(result)[, c('n', 'p', 'rho', 'model', 'method', 
	                                         'quar', 'elln', 
	                                         'Rej1', 'Rej2', 'Rej3', 'Rej4',
	                                         'test_stat', 
	                                         'mean_pval', 'med_pval', 'min_pval')]
	    }
	    return( result ) } )))
  }
}

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
# Used to save the simulation result in the r-th MC run (here r = 1)
save(sim, file = paste('sim_rdo_maxCorrSurv_cr', gsub("\\%", "", cr), '_', r, '.Rdata', sep = ''))


##### sim contains the simulation output. Each row contains:
# n: Sample size
# p: Covariate dimension
# rho: Designed correlation of covariates
# model: Data generating distributions
# method: Method run
# quar: The quantile of outcomes used for truncation
# elln: Value of elln used by our method
# Rej1: Whether rejecting the null of no correlation, comparing the mean_pval with the alpha_val
# Rej2: Whether rejecting the null of no correlation, comparing the med_pval with the alpha_val
# Rej3: Whether rejecting the null of no correlation, comparing the min_pval with the alpha_val
# Rej4: Whether rejecting the null of no correlation, comparing the upper tail prob 
#       of the asymptotic distribution of the test_stat with the alpha_val
# test_stat: Test statistic 
# mean_pval: the mean of the pvalues obtained from 10-fold random orderings of data
# med_pval: the median of the pvalues obtained from 10-fold random orderings of data
# min_pval: the minimum of the pvalues obtained from 10-fold random orderings of data

