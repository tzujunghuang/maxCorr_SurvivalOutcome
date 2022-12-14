##### This script is for one single MC simulation run, say r = 1.
##### To reproduce our numeric results, it is required that r = 1,...,1000 and 
##### it is highly suggested that parallel computing techniques should be used.


##### Configuration
# Read in source codes that should be in the same working directory as is this script	
source('code_maxCorrSurv.R')
source('sim_poweranalysis_data_acquiration.R')
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '10%'
# Methods to run
meths_vec = c('BONF_COX', 'BONF_OSE', 'OOSE', 'SOSE')
# Consider a list of (n,p) values
n_vals = c(500)
np_vals = data.frame('n' = 500, 'p' = 10^5)
# List of correlation values between predictors
rho_val = 0.75
# Tau
quar = 0.9
# Significance level
alpha_val = 0.05
# elln 
elln_part = c(2); dim_elln = length(elln_part)
# List of data generating distributions to use when running simulations
model_list = c('A2.DE','A1.DE','A2.IE','A1.IE')
# Number of potential coefficient values used in models
n_eta_vals = 5


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

# Number of times to run simulation
num_digits = 7; options(digits = num_digits)

# Considered settings
settings = data.frame('model' = model_list) %>% 
  mutate('eta_effect' = gsub('.DE','', gsub('.IE', '', model)), 
         'join_idx' = 1) %>%
  full_join(cbind(np_vals, 'join_idx' = 1), by = 'join_idx') %>% 
  select(-join_idx) %>%
  mutate(ntimes = n_eta_vals) %>%
  group_by(model, eta_effect, n, p) %>% slice(rep(1:n(), first(ntimes))) %>% 
  select(-ntimes) %>%
  bind_cols('eta' = c(rep(seq(0.1, 0.2, 0.025), 2*length(n_vals)), 
                      rep(seq(0.01, 0.02, 0.0025), 2*length(n_vals)))) %>%
  as.data.frame()


sim = data.frame(method = NA, model = NA, n = NA, p = NA, eta = NA, rho = NA, 
                 quar = NA, elln = NA, rej = NA, est = NA, se = NA)

for (row_idx in 1:nrow(settings)) {
    
    model = settings[row_idx, 'model']; eta_effect = settings[row_idx, 'eta_effect'];
    n = settings[row_idx, 'n']; p = settings[row_idx, 'p']; eta = settings[row_idx, 'eta']
    rho = rho_val
    
    dat = Simulate_data(n, p, model, censoring_rate = cr, rho = rho, eta, 
                        censoring_dist_par = c(0.07, 0.15, 0.28, 0.43))
    obj0 <- NumericalStudy$new(input_data = dat)
    print(c(n, p, model, eta))
    
    out = list()
    for (meth in meths_vec) {
      if (meth == 'SOSE') {
        
        print(meth) 
        SOSE_est = do.call(rbind, lapply(1:length(elln_part), function(d_indx){
          d = elln_part[d_indx]; elln = ceiling(n/d)
          chunk_size = ceiling((n - elln)/10)
          SOSE_est = obj0$Stab_onestep_est(all_obs = 1:n, chunk_size, 
                                           elln, est_index = est_index_val, 
                                           alpha = alpha_val, num_top = 1, 
                                           quar_trunc = quar)
          return( t( c( SOSE_est$ci, SOSE_est$est, SOSE_est$se, elln ) ) ) 
        } ) )
        out$SOSE = SOSE_est
        
      } else if (meth == 'OOSE') {
        
        print(meth)
        OOSE_est = obj0$Oracle_onestep_est(alpha = alpha_val, 
                                           quar_trunc = quar, idx_taken = 1)
        out$OOSE = c( OOSE_est$rej, OOSE_est$est, OOSE_est$se )
        
      } else if (meth == 'Naive_OSE') {
        
        print(meth)
        Naive_OSE_est = obj0$Naive_onestep_est(alpha = alpha_val, 
                                               quar_trunc = quar, num_top = 1)
        out$Naive_OSE = c( Naive_OSE_est$rej, Naive_OSE_est$est, Naive_OSE_est$se )
        
      } else if (meth == 'BONF_COX') {
        
        print(meth)
        BONF_COX_est = obj0$Bonf_Cox_ests(alpha = alpha_val)
        out$BONF_COX = c( BONF_COX_est$rej )
        
      } else if (meth == 'BONF_OSE') {
        
        print(meth)
        BONF_OSE_est = obj0$Bonf_onestep_ests(alpha = alpha_val, 
                                              quar_trunc = quar)
        out$BONF_OSE = c( BONF_OSE_est$rej, BONF_OSE_est$est, BONF_OSE_est$se )
        
      } else stop('Invalid method.') }
    gc()
    
    sim = rbind(sim, do.call(rbind,lapply(meths_vec, function(meth) {
      if ( meth == 'SOSE' ) {
        result = data.frame( method = rep(meth, dim_elln), model = rep(model, dim_elln),
                             n = rep(n, dim_elln), p = rep(p, dim_elln), eta = rep(eta, dim_elln),
                             rho = rep(rho, dim_elln), 
                             quar = rep(quar, dim_elln), elln = out[[meth]][,5],
                             rej = 1*( out[[meth]][,1] > 0 | out[[meth]][,2] < 0 ),
                             est = out[[meth]][,3], se = out[[meth]][,4] )
        
      } else if ( meth == 'BONF_COX' ) {
        result = data.frame( method = meth, model = model, n = n, p = p, 
                             eta = eta, rho = rho, quar = quar, elln = NA, 
                             rej = 1*(out[[meth]][1] > 0), est = NA, se = NA )
      } else {
        result = data.frame( method = meth, model = model, n = n, p = p, 
                             eta = eta, rho = rho, quar = quar, elln = NA, 
                             rej = 1*(out[[meth]][1] > 0), est = out[[meth]][2], 
                             se = out[[meth]][3] )
      }
      return( result ) })))
}

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
# Used to save simulation file
save(sim, file = paste('sim_poweranalysis_maxCorrSurv_cr', gsub("\\%", "", cr), '_', r, '.Rdata', sep = ''))


## sim contains the simulation output. Each row contains:
# method: Method run
# model: Data generating distributions
# n: Sample size
# p: Covariate dimension
# eta: Signal strength (a positive value)
# rho: Designed correlation of covariates
# quar: The quantile of outcomes used for truncation 
# elln: Value of elln used by our method
# rej: Whether rejecting the null of no correlation
# est: Estimate of the parameter
# se: Standard error of the statistic

