num_digits = 7
source('data_management.R')

# Read in VRC01 data and prepare analysis data
Get_analysis_data = function(csv_path, binary_vars, full, single_layer_stratification, 
                             outcome_label = 'ic50.geometric.mean.imputed', 
                             censoring_label = 'ic50.censored',
                             strata_var_1 = 'subtype.reduced', 
                             strata_val_1,
                             strata_var_2 = 'geographic.region.of.origin', 
                             strata_val_2, 
                             clinical_cols = c('subtype.reduced', 
                                               'country.of.origin', 
                                               'geographic.region.of.origin',
                                               'infection.stage.ordinal'),
                             start_index_predictors = 43){
  
  data_in = read.csv(file = csv_path)
  
  if (full) {
    data_in = data_in
  } else{
    if (single_layer_stratification) {
      data_in = data_in[data_in[strata_var_1] == strata_val_1, ] 
    } else { data_in = data_in[(data_in[strata_var_1] == strata_val_1) & 
                               (data_in[strata_var_2] == strata_val_2), ] }
  }  
  
  data = data_in; n = dim(data)[1]
  
  n_p_track = list(); n_p_track$n = n  
  
  # possible outcomes: 'ic50.geometric.mean.imputed', 'ic80.geometric.mean.imputed'
  X = data[ ,outcome_label];

  # censoring status: `ic50.censored`, `ic80.censored`
  delta = as.vector(t(1 - 1*data[ ,censoring_label]))
  
  # Clinical/Demographic variables
  U0 = data[ ,clinical_cols]
  
  # Genetic predictors
  # Remove a specific predictor
  U1 = data[ ,start_index_predictors:dim(data)[2]]
  n_p_track$p = dim(U1)[2]
  
  U1_void_idx0 = which(apply(U1, 2, function(x){ length(unique(x)) < 2 } ))
  U1_binary_idx0 = which(apply(U1, 2, function(x){ length(unique(x)) == 2 } ))
  U1_count_idx0 = which(apply(U1, 2, function(x){ (length(unique(x)) > 2 & all((x == floor(x)) & (x >= 0))) } ))
  U1_conti_idx0 = which(apply(U1, 2, function(x){ (length(unique(x)) > 2 & all(x > floor(x))) } ))
  n_p_track$p_void = length(U1_void_idx0); n_p_track$p_binary = length(U1_binary_idx0)
  n_p_track$p_count = length(U1_count_idx0); n_p_track$p_conti = length(U1_conti_idx0)
  # filepath0 = 'variable_names.csv';
  # ndf0 = data.frame(names(U0)); names(ndf0)[1] = 'ClinicalVariable'
  # res_list = check_void(U1, maineffect=TRUE)
  # df_info_U1 = res_list$df_info; res_list$void_percent
  # write.csv(ndf0, file=filepath0, sheetName='ClinicalVariables', row.names=FALSE)
  # write.csv(df_info_U1, file=filepath0, sheetName='GeneticVariables', append=TRUE, row.names=FALSE)
  
  U1_keep_idx = which(apply(U1, 2, function(x){( length(unique(x)) >= 2 & mostfrequent_percent(x) < 0.95)} ))
  n_p_track$p_valid = length(U1_keep_idx)
  
  # Check data type
  U1_binary_idx = which(apply(U1[U1_keep_idx], 2, function(x){ length(unique(x)) == 2 } ))
  U1_count_idx = which(apply(U1[U1_keep_idx], 2, function(x){ (length(unique(x)) > 2 & all((x == floor(x)) & (x >= 0))) } ))
  U1_conti_idx = which(apply(U1[U1_keep_idx], 2, function(x){ (length(unique(x)) > 2 & all(x > floor(x))) } ))
  n_p_track$p_valid_binary = length(U1_binary_idx); n_p_track$p_valid_count = length(U1_count_idx)
  n_p_track$p_valid_conti = length(U1_conti_idx)
  
  if (length(U1_count_idx) > 0 | length(U1_conti_idx) > 0) {
    taken_idx = c(U1_count_idx, U1_conti_idx)
    U2 = data.matrix(U1[U1_keep_idx][taken_idx])
    # Pre-standardization for continuous variable
    U2_colSDs = apply(U2, 2, sd); nn = dim(U2)[1]
    # U2 = (U2 - matrix(rep(U2_colMeans, nn), nrow=nn, byrow=TRUE))/matrix(rep(U2_colSDs, nn), nrow=nn, byrow=TRUE)
    U2 = colStandardization(U2, U2_colSDs)
    U2 = as.data.frame(U2)
    
    if (binary_vars) { U3 = U1[U1_keep_idx][U1_binary_idx] } else { U3 = U2 } 
  } else { if (binary_vars) { U3 = U1[U1_keep_idx][U1_binary_idx] } else { stop('No Available Data.') } }
  
  # Interactions of U1
  U1_inter = do.call(cbind, combn(colnames(U3), 2,
                                  FUN = function(x){ list( setNames(data.frame(U3[,x[1]]*U3[,x[2]]), 
                                                                    paste(x, collapse = "_")) ) } ))
  # res_list = check_void(U1_inter, maineffect=FALSE)
  # df_info_U1inter = res_list$df_info; res_list$void_percent
  # write.csv(df_info_U1inter, file=filepath0, sheetName='InteractionVariables', append=TRUE, row.names=FALSE)
  
  U1inter_keep_idx = which(apply(U1_inter, 2, function(x){( length(unique(x)) >= 2 & mostfrequent_percent(x) < 0.95)} ))
  U1_inter_valid = U1_inter[U1inter_keep_idx]
  n_p_track$p_inter = dim(U1_inter)[2]; n_p_track$p_inter_valid = length(U1inter_keep_idx)
  if (binary_vars) { U1_inter_valid_0 = U1_inter_valid
  } else { U1_inter_valid_colSDs = apply(U1_inter_valid, 2, sd); 
           U1_inter_valid_0 = colStandardization(U1_inter_valid, U1_inter_valid_colSDs) } 
  
  U = data.matrix(cbind(U3, U1_inter_valid_0))
  U_names = c(names(U3), names(U1_inter_valid))
  
  return( list(X = round(X, num_digits), delta = delta, U0 = data.matrix(U0), 
               U = U, U_names = U_names, n_p_track = n_p_track) )
}

Explore_data = function(data_used){
  
  filepath = 'data_exploration.xlsx';
  
  df1 = data.frame(table(data_used$U0$subtype.reduced)); names(df1)[1] = 'Virus_Subtype'
  df2 = data.frame(table(data_used$U0$country.of.origin)); names(df2)[1] = 'Country_Of_Origin'
  df3 = data.frame(table(data_used$U0$geographic.region.of.origin)); names(df3)[1] = 'GeographicRegion_Of_Origin'
  df4 = data.frame(table(data_used$U0$infection.stage.ordinal)); names(df4)[1] = 'Infection_Stage'
  
  write.csv(df1, file = filepath, sheetName = 'Virus_Subtype', row.names = FALSE)
  write.csv(df2, file = filepath, sheetName = 'Country_Of_Origin', append = TRUE, row.names = FALSE)
  write.csv(df3, file = filepath, sheetName = 'GeographicRegion_Of_Origin', append = TRUE, row.names = FALSE)
  write.csv(df4, file = filepath, sheetName = 'Infection_Stage', append = TRUE, row.names = FALSE)
}

