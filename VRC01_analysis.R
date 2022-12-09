##### This script is for the real data analysis that uses 'BONF_COX', 'BONF_OSE' and
##### 'SOSE' with one random ordering of data, that is, r = 1. 
##### Further values of r depend on how many random orderings will be taken. In the article, r=1,...,200.
# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr') #'doParallel',

for (package in package_list) {
  require(package, character.only = TRUE); library(package, character.only = TRUE) }

r = 1
# slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# r <- as.numeric(slurm_arrayid)
# set.seed(r)

# registerDoParallel(cores=3) # 3 cores on local machine; 6 on each node of clusters

# Should be in same working directory	
source('real_data_acquiration.R')
source('code_maxCorrSurv.R')
csv_path = 'fulldata.csv'

# Methods to run
if (r==1){ meth_vecs = c('BONF_COX', 'BONF_OSE', 'SOSE')
} else { meth_vecs = c('SOSE') }

# Configuration
cr = '16%'
elln_part = c(2); dim_elln = length(elln_part)
quars_list = c(0.9)

# Estimation Index
est_index_val = 'whole samples'

# Significance level
alpha_val = 0.05

# Read in analysis data
### Full data
#data0 = Get_analysis_data(csv_path, binary_var=TRUE, full=TRUE, single_layer_stratification=FALSE,
#                          outcome_label='ic50.geometric.mean.raw', censoring_label='ic50.censored',
#                          strata_var_1 = NA, strata_val_1=NA,
#                          strata_var_2 = NA, strata_val_2=NA, 
#                          clinical_cols=c('subtype.reduced', 'country.of.origin', 'geographic.region.of.origin',
#                                          'infection.stage.ordinal'),
#                          start_index_predictors=43)

# vals = c("Asia", "Europe.Americas", "N.Africa", "S.Africa")
# text_vals = c("Asia", "Europe_Americas", "N_Africa", "S_Africa")

# val = "Asia"; text_val = "Asia"
#data0 = Get_analysis_data(csv_path, binary_var=TRUE, full=FALSE, single_layer_stratification=TRUE,
#                          outcome_label='ic50.geometric.mean.raw', censoring_label='ic50.censored',
#                          strata_var_1='geographic.region.of.origin', strata_val_1=val,
#                          strata_var_2=NA, strata_val_2=NA,
#                          clinical_cols=c('subtype.reduced', 'country.of.origin', 'geographic.region.of.origin',
#                                          'infection.stage.ordinal'),
#                          start_index_predictors=43)

### Double layer stratification
#val = "S.Africa"; text_val = "subtype_C_S_Africa"
#data0 = Get_analysis_data(csv_path, binary_var=TRUE, full=FALSE, single_layer_stratification=FALSE,
#                          outcome_label='ic50.geometric.mean.raw', censoring_label='ic50.censored',
#                          strata_var_1='subtype.reduced', strata_val_1='C',
#                          strata_var_2='geographic.region.of.origin', strata_val_2=val,
#                          clinical_cols=c('subtype.reduced', 'country.of.origin', 'geographic.region.of.origin',
#                                          'infection.stage.ordinal'),
#                          start_index_predictors=43)

#Explore_data(data_used=data0)

### Single layer stratification
vals = c("B", "C"); types = c("c", "b")
vals_types = expand.grid(vals, types); colnames(vals_types) = c("val", "type")
texts = sapply(1:dim(vals_types)[1], 
                   function(i){ paste0("subtype_", vals_types[i,"val"], "_", vals_types[i,"type"]) })

sim = data.frame(val = NA, type = NA, method = NA, rej = NA, p_val = NA, est = NA, se = NA, 
                 Sn = NA, k_est = NA, lb_ci = NA, ub_ci = NA, quar = NA, elln = NA)


for (idx in 1:dim(vals_types)[1]) {
  
  print(idx)
  
  val = as.character(vals_types[idx,"val"]); type = as.character(vals_types[idx,"type"])
  is.binary_vars = ifelse(type == 'b', TRUE, FALSE)
  
  data0 = Get_analysis_data(csv_path, binary_vars = is.binary_vars, full = FALSE, 
                            single_layer_stratification = TRUE,
                            outcome_label = 'ic50.geometric.mean.imputed', 
                            censoring_label = 'ic50.censored',
                            strata_var_1 = 'subtype.reduced', strata_val_1 = val,
                            strata_var_2 = NA, strata_val_2 = NA,
                            clinical_cols = c('subtype.reduced', 'country.of.origin', 
                                              'geographic.region.of.origin',
                                              'infection.stage.ordinal'),
                            start_index_predictors = 43)

  
  n = dim(data0$U)[1]; p = dim(data0$U)[2]
  obj0 <- NumericalStudy$new(input_data = data0)

  for (quar in quars_list) {
  
   out = list();
   for (meth in meth_vecs) {
      if (meth == 'SOSE') {
        
        print(meth)
        
        inds = sample(1:n)
        data1 = list(X = round(data0$X[1:n][inds], num_digits), delta = data0$delta[1:n][inds], 
                     U = data0$U[1:n,1:p][inds,])
        obj1 <- NumericalStudy$new(input_data = data1)
        
        SOSE_est = do.call(rbind, lapply(1:length(elln_part), function(d_indx){
          d = elln_part[d_indx]; elln = ceiling(n/d)
          chunk_size = ceiling((n - elln)/10)
          SOSE_est = obj1$Stab_onestep_est(all_obs = 1:n, chunk_size, elln, est_index = est_index_val, 
                                           alpha = alpha_val, num_top = 1, quar_trunc = quar)
          return( t( c( SOSE_est$rej, SOSE_est$p_val, SOSE_est$Sn, SOSE_est$est, SOSE_est$se, SOSE_est$k_est, 
                        SOSE_est$ci, elln )  ) ) 
        } ) )
        
        out$SOSE = SOSE_est
        
      } else if (meth == 'BONF_COX') {
        
        print(meth)
        BONF_COX_est = obj0$Bonf_Cox_ests(alpha = alpha_val)
        out$BONF_COX = c( BONF_COX_est$rej, BONF_COX_est$p_val, BONF_COX_est$k_est )
     
      } else if (meth == 'BONF_OSE') {
        
        print(meth)
        data2 = list(X = data0$X[1:n], delta = data0$delta[1:n], U = data0$U[1:n,1:p])
        obj2 <- NumericalStudy$new(input_data = data2)
        
        BONF_OSE_est = obj2$Bonf_onestep_ests(alpha = alpha_val, quar_trunc = quar)
        out$BONF_OSE = c( BONF_OSE_est$rej, BONF_OSE_est$p_val, BONF_OSE_est$Sn, BONF_OSE_est$est, 
                          BONF_OSE_est$se, BONF_OSE_est$k_est, BONF_OSE_est$k_ests )
        
      } else stop('Invalid method.') }
 
   
   sim = rbind(sim, do.call(rbind,lapply(meth_vecs, function(meth) {
    
    if (meth == 'SOSE') {
      
      runs = as.matrix(out[[meth]])
      result = data.frame( val = val, type = type, method = meth, rej = runs[,1], p_val = runs[,2], 
                           est = runs[,3], se = runs[,4], Sn = runs[,5], k_est = data0$U_names[runs[,6]], 
                           lb_ci = runs[,7], ub_ci = runs[,8], quar = quar, elln = runs[,9] )

    } else if (meth == 'BONF_OSE') {
      
      runs = t(as.matrix(out[[meth]]))
      result = data.frame( val = val, type = type, method = meth, rej = runs[,1], p_val = runs[,2], 
                           est = runs[,3], se = runs[,4], Sn = runs[,5], k_est = data0$U_names[runs[,6]], 
                           lb_ci = NA, ub_ci = NA, quar = quar, elln = NA )

    } else if (meth == 'BONF_COX') {
      
      runs = t(as.matrix(out[[meth]]))
      result = data.frame( val = val, type = type, method = meth, rej = runs[,1], p_val = runs[,2], 
                           est = NA, se = NA, Sn = NA, k_est = data0$U_names[runs[,3]], 
                           lb_ci = NA, ub_ci = NA, quar = NA, elln = NA )

    }
    
    return( result )  }  ) ) ) 
  }
}  

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
save(sim, file = paste0('VRC01_subdata_analysis_', r, '.Rdata'))

