##### This script is used to collect the results from the parallel computing implemented on 1000 nodes.

##### Configuration
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr_num = 10
# Methods to run
meths_vec = c('BONF_COX', 'BONF_OSE', 'SOSE', 'OOSE')
# Consider a list of (n,p) values to achieve the desired 90% power 
n_vals = c(500)
np_vals = data.frame('n' = 500, 'p' = 10^5)
# List of correlation values between predictors
rho_val = 0.75; rho_char = '075'
# Tau
quars_vec = c(0.9)
# Significance level
alpha_val = 0.05
# elln 
elln_part = c(2); elln_char = 'n05'; dim_elln = length(elln_part)
elln_quar_list = expand.grid(elln_part, quars_vec)
names(elln_quar_list) = c('d', 'quar')
# List of data generating distributions to use when running simulations
model_list = c('A2.DE','A1.DE','A2.IE','A1.IE')
# Number of potential coefficient values used in models
n_eta_vals = 5

# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr') #'foreach', 'doParallel'
# new_packages = package_list[!(package_list %in% installed.packages()[,"Package"])]
# if (length(new_packages)) { install.packages(new_packages) } 
for (package in package_list) {
  require(package, character.only = TRUE); library(package, character.only = TRUE) }

# Number of times to run simulation
num_digits = 7; options(digits = num_digits)

# Considered settings 
settings = data.frame('model' = model_list) %>% mutate('eta_effect' = gsub('.DE','', gsub('.IE', '', model)), 
                                                       'join_idx' = 1) %>%
  full_join(cbind(np_vals, 'join_idx' = 1), by = 'join_idx') %>% select(-join_idx) %>%
  mutate(ntimes = n_eta_vals) %>%
  group_by(model, eta_effect, n, p) %>% slice(rep(1:n(), first(ntimes))) %>% select(-ntimes) 

scenarios = list( 'spec_1' = list(seq(0.1, 0.2, 0.025), seq(0.01, 0.02, 0.0025)),
                  'spec_2' = list(seq(0.1, 0.5, 0.1), seq(0.01, 0.05, 0.01)) )

#1000 MC runs for each scenario
num_scenarios = 1; num_r = 1000

for (spec_idx in 1:num_scenarios) { print( paste('spec_idx = ', spec_idx, sep = '') )
  
  vecs = scenarios[[spec_idx]]

  settings_spec = settings %>% bind_cols('eta' = c(rep(vecs[[1]], 2*length(n_vals)), 
                                                   rep(vecs[[2]], 2*length(n_vals)))) %>% as.data.frame()

  orig_seq = 1:num_r  ## seq(1+(spec_idx-1)*num_r, spec_idx*num_r, 1); 
  
  final = do.call(rbind, lapply(orig_seq, function(r_idx){
    file = paste('sim_poweranalysis_maxCorrSurv_cr', cr_num, '_', r_idx, '.Rdata', sep = '')
    load(file)
    if (r_idx %in% c(250, 500, 750, 1000)) {print(r_idx)} 
  
    temp_final = NULL   
    for (meth in meths_vec) {
      for (row_idx in 1:nrow(settings)) {
      
        model = settings_spec[row_idx, 'model']; 
        eta_effect = settings_spec[row_idx, 'eta_effect'];
        n = settings_spec[row_idx, 'n']; p = settings_spec[row_idx, 'p']; 
        eta = settings_spec[row_idx, 'eta']
      
        if (meth == 'SOSE') {
          #print(meth)
          results = do.call(rbind, lapply(1:dim(elln_quar_list)[1], function(row_idx){
            
            d = elln_quar_list[row_idx, 'd']; quar = elln_quar_list[row_idx, 'quar']
            elln = ceiling(n/d); 
            
            condition1 = ( sim$method == meth & sim$model == model & sim$n == n & 
                           sim$p == p & sim$eta == eta & sim$rho == rho_val & 
                           sim$quar == quar & sim$elln == elln )

            est = sim[condition1, 'est']; se =  sim[condition1, 'se']
            rej =  sim[condition1, 'rej']
            ret = data.frame(model = model, meth = meth, n = n, p = p, 
                             eta = eta, rho = rho_val, quar = quar, elln = elln, 
                             rej = rej, est = est, se = se, replic_num = r_idx)
            return( ret ) }))
        } else { 
          #print(meth)
          results = do.call(rbind, lapply(1:length(quars_vec), function(idx){
            
            quar = quars_vec[idx]
            condition1 = ( sim$method == meth & sim$model == model & 
                           sim$n == n & sim$p == p & sim$eta == eta &
                           sim$rho == rho_val & sim$quar == quar )
            est = sim[condition1, 'est']; se = sim[condition1, 'se']
            rej =  sim[condition1, 'rej']
            ret = data.frame(model = model, meth = meth, n = n, p = p, 
                             eta = eta, rho = rho_val,quar = quar, elln = NA, 
                             rej = rej, est = est, se = se, replic_num = r_idx)
            return( ret ) })) }    
        temp_final = rbind(temp_final, results)
    }        
  }  
  return(temp_final) }))

  if (anyNA(final$rej)) { print(r_idx) }  
  final$est_index = est_index_val
  save(final, file = paste('results', spec_idx, '_poweranalysis_maxCorrSurv_cr', 
                           cr_num, '_', elln_char, '_', est_index_val, '.Rdata', sep = ''))

}

