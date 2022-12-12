##### This script is used to collect the results from the parallel computing implemented on 1000 nodes.

##### Configuration
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '10%'; cr_num = 10
# Methods to run
meths_vec = c('BONF_COX', 'Naive_OSE', 'BONF_OSE', 'OOSE', 'SOSE')
# List of (n,p) values to run simulation
np_list = list(c(500, 1e6), c(500, 1e5), c(500, 1e4), c(500, 1000), c(500, 100))
# List of correlation values between predictors
rho_val = 0.75; rho_char = '075'
# Tau
quars_vec = c(0.9)
# Significance level
alpha_val = 0.05
# elln 
elln_part = c(2); dim_elln = length(elln_part); elln_char = 'n05'
# List of data generating distributions to use when running simulations
model_list = c('N.IE','A1.IE','A2.IE', 'N.DE','A1.DE','A2.DE')

elln_quar_list = expand.grid(elln_part, quars_vec)
names(elln_quar_list) = c('d', 'quar')

num_r = 1000
orig_seq = 1:num_r;

final = do.call(rbind, lapply(orig_seq, function(r_idx){
  file = paste('sim_maxCorrSurv_cr', cr_num, '_', r_idx, '.Rdata', sep = '')
  load(file)
  if (r_idx %in% c(250, 500, 750, 1000)) { print(r_idx) } 
  
  temp_final = NULL   
  for (meth in meths_vec) {
    for (model in model_list) {
      for (np in np_list) {
        
        n = np[1]; p = np[2]    
        if (meth == 'SOSE') {
          #print(meth)
          results = do.call(rbind, lapply(1:dim(elln_quar_list)[1], function(row_idx){
            
            d = elln_quar_list[row_idx, 'd']; quar = elln_quar_list[row_idx, 'quar']
            elln = ceiling(n/d); 
            
            condition1 = ( sim$method == meth & sim$model == model & sim$n == n & 
                           sim$p == p & sim$rho == rho_val &
                           sim$quar == quar & sim$elln == elln )

            est = sim[condition1, 'est']; se =  sim[condition1, 'se']
            rej =  sim[condition1, 'rej']
            ret = data.frame(model = model, rho = rho_val, meth = meth, 
                             n = n, p = p, quar = quar, elln = elln, 
                             est = est, se = se, rej = rej, replic_num = r_idx)
            return( ret ) }))
        } else { 
          #print(meth)
          results = do.call(rbind, lapply(1:length(quars_vec), function(idx){
            
            quar = quars_vec[idx]
            condition1 = ( sim$method == meth & sim$model == model & 
                           sim$n == n & sim$p == p & sim$rho == rho_val &
                           sim$quar == quar )
            est = sim[condition1, 'est']; se = sim[condition1, 'se']
            rej =  sim[condition1, 'rej']
            ret = data.frame(model = model, rho = rho_val, meth = meth, 
                             n = n, p = p, quar = quar, elln = NA, 
                             est = est, se = se, rej = rej, replic_num = r_idx)
            return( ret ) })) }    
        temp_final = rbind(temp_final, results)
      }      
    }        
  }  
  return(temp_final) }))


if (anyNA(final$rej)) { print(r_idx) }  
final$est_index = est_index_val
save(final, file = paste('results_maxCorrSurv_cr', cr_num, '_', est_index_val, '.Rdata', sep = ''))

