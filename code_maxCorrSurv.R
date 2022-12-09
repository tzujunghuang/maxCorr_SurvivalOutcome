num_digits = 7
source('data_management.R')

NumericalStudy <- setRefClass( "NumericalStudy",
 fields = list(
    input_data = "list",
    n = "numeric",
    p = "numeric"
  ),
  
 methods = list(
  initialize = function(input_data = dat){
    input_data <<- input_data
    n <<- dim(input_data$U)[1]
    p <<- dim(input_data$U)[2]
  },
    
  KM_weight_func = function(obs){ # obs is a vector of observation/individual indices 
    data_km = data.frame(X=input_data$X[obs], delta=input_data$delta[obs])
    data_km = data_km[order(data_km$X),]
    prod = c(1, cumprod( sapply(1:(n-1), function(i){ ((n-i)/(n-i+1))^(data_km$delta[i]) }) ) )
    kmwts = sapply(1:n, function(i){ data_km$delta[i]*prod[i]/(n-i+1) } )
    return ( kmwts )
  },  
    
  KM_SurF = function(t, obs){ # t is a scalar; obs is a vector of observation/individual indices 
    data_km = data.frame(X=input_data$X[obs], delta=input_data$delta[obs])
    km = survfit(Surv(X, 1-delta)~1, data=data_km)
    rm(data_km)
    
    survest = cbind(km$time, km$surv)
    if ( length(which(survest[, 1] <= t)) > 0 ) {
        return( survest[max(which(survest[, 1] <= t)), 2] )  
    } else { return( 1 ) }
  },  
  
  KM_SurF_self = function(obs, quar_trunc){ # obs is a vector of observation/individual indices 
    data_km = data.frame(X=input_data$X[obs], delta=input_data$delta[obs])
    tau = quantile(data_km$X, probs = quar_trunc)
    tau_surv = KM_SurF(tau, obs)
    
    km = survfit(Surv(X, 1-delta)~1, data=data_km)
    rm(data_km)
    
    km$surv[km$time > tau] = tau_surv
    survest = cbind(km$time, km$surv)
    return ( survest )
  },
  
  Inverse_weight_func = function(x0, delta0, obs0, quar_trunc, err_msg){
    tau = quantile(x0, probs = quar_trunc)
    tau_surv = KM_SurF(tau, obs=obs0)

    if (length(x0) == 1) {
      inverse_weight = ifelse(x0 < tau, KM_SurF(x0, obs=obs0), tau_surv)
    } else {    
      if (length(unique(x0))==length(x0)) { KM_table = KM_SurF_self(obs=obs0, quar_trunc) 
      } else {
        KM_table0 = data.frame(KM_SurF_self(obs=obs0, quar_trunc))
        names(KM_table0) = c("time", "surv_prob")
        n_occur = data.frame(table(x0)); names(n_occur) = c("time", "freq")
        n_freq = n_occur$freq; n1 = n_freq[n_freq > 1]; n1[is.nan(n1)|is.na(n1)|is.null(n1)] = 0
        
        tryCatch({ temp = data.frame(cbind(rep(KM_table0[KM_table0$time %in% n_occur$time[n_freq > 1],]$time, n1),
                                           rep(KM_table0[KM_table0$time %in% n_occur$time[n_freq > 1],]$surv_prob, n1))) },
          error = function(err_msg) {
            message("Original error message:"); message(paste0(err_msg,"\n"))
            save( list(x0=x0, delta0=delta0), file=paste0('errordata_', cr, '_', r, '.Rdata') )
            stop(err_msg) })
        names(temp) = names(KM_table0)
        # times_once = as.numeric( as.character(n_occur$time[n_freq == 1]) )
        KM_table = as.matrix(rbind(temp, KM_table0[KM_table0$time %in% n_occur$time[n_freq == 1],])) }
        
      diff = setdiff(x0, KM_table[,1])  
      if (length(diff) > 0) {
        for (item in diff) { KM_table = rbind(KM_table, c(item, KM_table[max(which(item >= KM_table[,1])), 2])) } }
      KM_table[KM_table[,1] > tau, 2] = tau_surv
      inverse_weight = as.vector(KM_table[order(KM_table[,1]),][,2]) }  
    
    if (length(inverse_weight) != length(x0)) {
      save( list(x0=x0, delta0=delta0), file=paste('errordata_', cr, '_', r, '.Rdata', sep=''))
      stop('Dimension Problem!')  
    } else { return( inverse_weight ) }
  },
  
  Est_Psi0_d = function(U_index, obs, quar_trunc){
    # U_index is a scalar or a vector of predictor indices; 
    # obs is a vector of observation/individual indices  
    data0 = cbind(input_data$X[obs], input_data$delta[obs], input_data$U[obs, U_index])
    data0 = data0[order(data0[,1]),]
    
    x_0 = data0[,1]; delta_0 = data0[,2]
    selectU_0 = data0[,3:dim(data0)[2]]
    rm(data0)
    
    inverse_weight = Inverse_weight_func(x0=x_0, delta0=delta_0, obs0=obs, quar_trunc, 
                                         err_msg='Error in Est_Psi0_d Inverse_weight KM_table')
    Y_Ghat = (x_0*delta_0) / inverse_weight; Y_Ghat[is.nan(Y_Ghat)|is.na(Y_Ghat)] = 0
    cov_Y_selectU = colCovs(selectU_0, Y_Ghat); var_selectU = colVars(selectU_0, sd_use=FALSE)
    est = cov_Y_selectU / var_selectU
    return( est ) 
  },
  
  IF_star_self = function(m, U_index, obs, quar_trunc){
    # m is a scalar; U_index is a scalar or a vector of predictor indices; 
    # obs is the vector of old/whole individuals indices
    data0 = cbind(input_data$X[obs], input_data$delta[obs], input_data$U[obs, U_index])
    data0 = data0[order(data0[,1]),]
    X = data0[,1]; delta = data0[,2]; selectU_0 = data0[,3:dim(data0)[2]]
    rm(data0)
    tau = quantile(X, probs = quar_trunc)
    inverse_weight = Inverse_weight_func(x0=X, delta0=delta, obs0=obs, quar_trunc, 
                                         err_msg='Error in IF_star_self Inverse_weight KM_table')

    if ( is.matrix(selectU_0) ) {
       mu_selectU_0 = matrix(rep(my_colMeans(selectU_0), length(obs)), 
                             byrow=TRUE, nrow=length(obs))
       selectU_0_colVars = colVars(selectU_0, sd_use=FALSE)
    } else { mu_selectU_0 = mean(selectU_0) 
             selectU_0_colVars = var(selectU_0) }
    
    var_selectU_0 = matrix(rep(selectU_0_colVars, length(obs)), byrow=TRUE, nrow=length(obs))
    gc()
        
    Y_Ghat = (X*delta) / inverse_weight; Y_Ghat[is.nan(Y_Ghat)|is.na(Y_Ghat)] = 0
    mean_YGhat = mean(Y_Ghat); 
    
    if ( is.matrix(selectU_0) ) { cov_mat = colCovs(selectU_0, Y_Ghat)
    } else { cov_mat = cov(selectU_0, Y_Ghat) }
    
    a = ((selectU_0-mu_selectU_0)*(Y_Ghat-mean_YGhat)) / var_selectU_0
    b = ((selectU_0-mu_selectU_0)^2)*matrix(rep(cov_mat, length(obs)), nrow=length(obs), 
                                            byrow=TRUE) / (var_selectU_0)^2
    IF_ipw = matrix(rep(m,length(obs)), nrow=length(obs), byrow=TRUE)*(a-b)
    gc(); rm(a); rm(b)
    
    time_comparison1 = outer(X, X, '=='); time_comparison2 = outer(X, X, '>=')
    event_nums = colSums(time_comparison1*(1-delta)); risk_set = colSums(time_comparison2)
    Y_Ghat_mat = ( time_comparison2 * (X*delta) / inverse_weight ) 
    Y_Ghat_mat[is.nan(Y_Ghat_mat)|is.na(Y_Ghat_mat)] = 0
    mean_YGhat_mat = colMeans(Y_Ghat_mat) 
    
    ### Slow when obs is large
    if ( is.matrix(selectU_0) ) { 
    
     Ehat_pointwise_val = sapply( 1:length(obs), function(i){
       a = mean_YGhat_mat[i]; n_s = length(obs)-(i-1); p_s = length(U_index)
       if (i==length(obs)) { 
           selectU = t(as.matrix(selectU_0[length(obs),]))
       } else{ selectU = selectU_0[i:length(obs),] }
       
       tryCatch({ mean_selectU = my_colMeans(selectU); var_selectU = colVars(selectU, sd_use=FALSE) },
           error = function(msg='Error in IF_star_self Ehat_pointwise_val') {
           message("Original error message:"); message(paste0(msg,"\n"))
           save( input_data, file=paste0('errordata_in_selectU_for', i, '_in', r, '.Rdata') )
           stop(msg) })
       
       cov_Y_selectU = my_colMeans( (selectU - matrix(rep(mean_selectU, n_s), byrow=TRUE, ncol=p_s))
                                    *matrix(rep(Y_Ghat_mat[i:length(obs),i], p_s), ncol=p_s) )
       b = cov_Y_selectU / var_selectU; b[is.nan(b)|is.na(b)] = 0
       c = selectU_0[i,]-mean_selectU
       rm(selectU); rm(mean_selectU); rm(var_selectU); rm(cov_Y_selectU)
       return ( a+c*b ) } )
     Ehat_pointwise_val[is.nan(Ehat_pointwise_val)|is.na(Ehat_pointwise_val)] = 0
    
    } else {
      
     Ehat_pointwise_val = sapply( 1:length(obs), function(i){
        a = mean_YGhat_mat[i]; n_s = length(obs)-(i-1); p_s = length(U_index)
        if (i==length(obs)) { 
          selectU = selectU_0[length(obs)]
        } else{ selectU = selectU_0[i:length(obs)] }
        
        mean_selectU = mean(selectU);
        if (length(selectU) > 1) { var_selectU = var(selectU); 
          cov_Y_selectU = mean( (selectU - matrix(rep(mean_selectU, n_s), byrow=TRUE, ncol=p_s))
                                *matrix(rep(Y_Ghat_mat[i:length(obs),i], p_s), ncol=p_s) )                           
          b = cov_Y_selectU / var_selectU
          b[is.nan(b)|is.na(b)] = 0
        } else { var_selectU = 'NA'; cov_Y_selectU = 0; b = 0 } 
        
        c = selectU_0[i]-mean_selectU
        rm(selectU); rm(mean_selectU); rm(var_selectU); rm(cov_Y_selectU)
        return ( a+c*b ) } )
     Ehat_pointwise_val[is.nan(Ehat_pointwise_val)|is.na(Ehat_pointwise_val)] = 0 
      
    } 
    
    hazard = event_nums/risk_set; hazard[is.nan(hazard)|is.na(hazard)]=0
    mart_X = (X <= tau & delta==0) - (X <= tau)*hazard
    val_1 = ( matrix(t(Ehat_pointwise_val), nrow=length(obs))*mart_X )
    IF_CAR = ( matrix(rep(m,length(obs)), nrow=length(obs), byrow=TRUE)
               * (selectU_0-mu_selectU_0) * val_1 / var_selectU_0 )
    return( round((IF_ipw - IF_CAR), 3) )  
  },
  
  IF_star = function(m, U_index, obs_all, obs0, obs1, quar_trunc){
    # m is a scalar; U_index is a scalar or a vector of predictor indices; 
    # (obs_all, obs0, obs1) are the vectors of whole, old and new individual indices
    data_all = cbind(input_data$X[obs_all], input_data$delta[obs_all])
    data_all = data_all[order(data_all[,1]),]
    x_all = data_all[,1]; delta_all = data_all[,2]
    rm(data_all)
    
    data0 = cbind(input_data$X[obs0], input_data$delta[obs0], input_data$U[obs0, U_index])
    data0 = data0[order(data0[,1]),]
    x_0 = data0[,1]; delta_0 = data0[,2]
    selectU_0 = matrix(data0[,3:dim(data0)[2]], nrow=length(obs0), ncol=length(U_index))
    rm(data0)
    
    dup_vals_xall = x_all[duplicated(x_all)]; dup_vals_x0 = x_0[duplicated(x_0)]
    extra_dup_vals = setdiff(dup_vals_xall, dup_vals_x0)
    if (length(dup_vals_xall) > length(dup_vals_x0) & length(dup_vals_x0) == 0) {
      indices = order(x_all)[!duplicated(x_all) & x_all %in% x_0]
    } else if (length(dup_vals_xall) > length(dup_vals_x0) & length(dup_vals_x0) > 0) {
      indices_to_drop = order(x_all)[duplicated(x_all) & x_all %in% extra_dup_vals]
      indices = setdiff(order(x_all)[x_all %in% x_0], indices_to_drop)
    } else if (length(dup_vals_xall) == length(dup_vals_x0)) { indices = order(x_all)[x_all %in% x_0] }
    
    tau = quantile(x_all, probs = quar_trunc)
    tau_surv = KM_SurF(tau, obs=obs_all)
    
    data1 = list(X=input_data$X[obs1], delta=input_data$delta[obs1], U=input_data$U[obs1, U_index])
    X = data1$X; delta = data1$delta; selectU_1 = matrix(data1$U, nrow=length(obs1), ncol=length(U_index))
    rm(data1)
    
    if (length(obs1) == 1){
        KM_SurF_terms = ifelse(X < tau, KM_SurF(X, obs=obs_all), tau_surv)
    } else {    
        KM_SurF_terms = sapply(X, KM_SurF, obs=obs_all)
        KM_SurF_terms[X >= tau] = tau_surv }  
    
    mu_selectU_0 = matrix(rep(my_colMeans(selectU_0), length(obs1)), nrow=length(obs1), byrow=TRUE)
    selectU_0_colVars = colVars(selectU_0, sd_use=FALSE)
    var_selectU_0 = matrix(rep(selectU_0_colVars, length(obs1)), nrow=length(obs1), byrow=TRUE)
    
    inverse_weight = Inverse_weight_func(x0=x_all, delta0=delta_all, obs0=obs_all, quar_trunc, 
                                         err_msg='Error in IF_star Inverse_weight KM_table')
    Y_Ghat_0 = (x_all*delta_all) / inverse_weight; Y_Ghat_0[is.nan(Y_Ghat_0)|is.na(Y_Ghat_0)] = 0
    Y_Ghat_0 = Y_Ghat_0[indices]
    mean_YGhat_0 = mean(Y_Ghat_0)
    
    time_comparison1 = outer(x_all, X, '==')
    delta_matrix = 1-delta_all
    time_comparison2 = outer(x_all, X, '>=')
    tryCatch({
      event_nums = colSums(time_comparison1*delta_matrix); risk_set = colSums(time_comparison2)
      inv_weight_mat = inverse_weight
      Y_Ghat_mat = time_comparison2*x_all*delta_all / inv_weight_mat
      Y_Ghat_mat[is.nan(Y_Ghat_mat)|is.na(Y_Ghat_mat)] = 0
      Y_Ghat_mat = Y_Ghat_mat[indices,]
      mean_YGhat_mat = colMeans(Y_Ghat_mat)
    
      cov_mat_0 = colCovs(selectU_0, Y_Ghat_0)
      Y_Ghat_1 = delta*X/KM_SurF_terms; Y_Ghat_1[is.nan(Y_Ghat_1)|is.na(Y_Ghat_1)] = 0
    },
    error = function(msg='Error in IF_star time_comparison dimension') {
      message("Original error message:"); message(paste0(msg,"\n"))
      save( input_data, file=paste0('errdata_', r, '.Rdata') )
      stop(msg) })
    
    a = ((selectU_1-mu_selectU_0) * (Y_Ghat_1-mean_YGhat_0)) / var_selectU_0
    b = ((selectU_1-mu_selectU_0)^2) * matrix(rep(cov_mat_0, length(obs1)), nrow=length(obs1), byrow=TRUE) / (var_selectU_0)^2
    IF_ipw = matrix(rep(m,length(obs1)), nrow=length(obs1), byrow=TRUE)*(a-b)
    rm(a); rm(b)
    
    ### Slow when obs1 is large;
    Ehat_pointwise_val = sapply(1:length(obs1), function(i){
       a = mean_YGhat_mat[i]
       selectU = matrix(selectU_0[time_comparison2[indices,i],], ncol=length(U_index))
       n_s = dim(selectU)[1]; p_s = dim(selectU)[2]

       if (n_s==0) { NaN
       } else {
         mean_selectU = my_colMeans(selectU); var_selectU = colVars(selectU, sd_use=FALSE)
         cov_Y_selectU = my_colMeans((selectU - matrix(rep(mean_selectU, n_s), byrow=TRUE, ncol=p_s))
                                      *matrix(rep(Y_Ghat_mat[time_comparison2[indices,i],i], p_s), ncol=p_s))
         b = cov_Y_selectU / var_selectU; b[is.nan(b)|is.na(b)] = 0
         u = matrix(selectU_1[i,], byrow=TRUE, ncol=length(U_index))
         c = u-matrix(rep(mean_selectU, dim(u)[1]), byrow=TRUE, ncol=dim(u)[2])
         rm(u)
         return ( a+c*b ) }})
    
    Ehat_pointwise_val[is.nan(Ehat_pointwise_val)|is.na(Ehat_pointwise_val)] = 0
    rm(selectU_0)
    
    hazard = event_nums/risk_set; hazard[is.nan(hazard)|is.na(hazard)]=0
    mart_X = (X <= tau & delta==0) - (X <= tau)*hazard
    val_1 = ( matrix(t(Ehat_pointwise_val), nrow=length(obs1), ncol=length(U_index)) 
              * matrix(rep(mart_X, length(U_index)), nrow=length(obs1), ncol=length(U_index)) )
    IF_CAR = ( matrix(rep(m,length(obs1)), nrow=length(obs1), byrow=TRUE)
               * (selectU_1-mu_selectU_0) * val_1 / var_selectU_0 ) 
    rm(selectU_1)
    return( round((IF_ipw - IF_CAR), 3) )  
  },
    
  Stab_onestep_est = function(all_obs, chunk_size, elln, est_index, alpha, num_top=1, quar_trunc){ # chunk_size, elln, alpha are scalars.
    mt = rowMeans( sapply(0:(ceiling((n-elln)/chunk_size)-1), function(i){
                  # if (i==(ceiling((n-elln)/chunk_size)-1)) {print(i)}
                  # print(i)
                  old_obs = all_obs[1:(elln + i*chunk_size)]
                  if(i<ceiling((n-elln)/chunk_size)-1){
                     new_obs = all_obs[(elln + i*chunk_size+1):(elln + (i+1)*chunk_size)]
                  } else {
                     new_obs = all_obs[(elln + i*chunk_size+1):n]  }
                  
                  cors0 = Est_Psi0_d(U_index=1:p, obs=old_obs, quar_trunc); sgn0 = rep(1, p)
                     
                  ### EXPENSIVE!!!
                  if (p <= 1e4){
                    IF_star_mat = IF_star_self(m=sgn0, U_index=1:p, obs=old_obs, quar_trunc)
                    utility = (cors0 + my_colMeans(IF_star_mat))
                  } else { utility = cors0 }   
                     
                  if (num_top < p) {
                    k0 = order(abs(utility), decreasing=TRUE)[seq(num_top)]
                  } else { k0 = seq(p) }
        
                  cor_k0 = cors0[k0]
                  sgn_k0 = 2*(utility[k0]>=0) - 1
        
                  Uk0 = input_data$U[old_obs, k0, drop=FALSE]
        
                  mu_Uk0 = colMeans(Uk0)
                  sd_Uk0 = apply(Uk0, 2, sd)
            
                  if (est_index == 'subsamples') { obs_all_val=old_obs; obs0_val=old_obs
                  } else if (est_index == 'partial subsamples') { obs_all_val=all_obs; obs0_val=old_obs
                  } else if (est_index == 'whole samples') { obs_all_val=all_obs; obs0_val=all_obs }
      
         curr_sigma_inv0 = 1/sd(rowMeans(IF_star(m=sgn_k0, U_index=k0, obs_all=obs_all_val, obs0=obs0_val, obs1=old_obs, quar_trunc)))
         est0 = mean( (mean(sgn_k0*cor_k0) 
                       + rowMeans(IF_star(m=sgn_k0, U_index=k0, obs_all=obs_all_val, obs0=obs0_val, obs1=new_obs, quar_trunc)))
                       *curr_sigma_inv0 )
    # main_term = mean(sgn_k0*cor_k0) 
    # mean_pred_IF_star = mean(rowMeans(IF_star(m=sgn_k0, U_index=k0, obs_all=obs_all_val, obs0=obs0_val, obs1=new_obs, quar_trunc)))
    # mean_train_IF_star = mean(rowMeans(IF_star(m=sgn_k0, U_index=k0, obs_all=obs_all_val, obs0=obs0_val, obs1=old_obs, quar_trunc)))
         return( c(est0, curr_sigma_inv0) )  } ) )
       est = mt[1]/mt[2]
       se = 1/(mt[2]*sqrt(n-elln))
       ci = c( est-qnorm(1-alpha/2)/(mt[2]*sqrt(n-elln)), est+qnorm(1-alpha/2)/(mt[2]*sqrt(n-elln)) )
       rej = 1*( ci[1]>0 | ci[2]<0 )
       p_val = 2*pnorm(abs(est/se), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
       
       
    ###### To be checked   
       k_est = sapply(n:n, function(i){
         old_obs = all_obs[1:n]
         
         cors0 = .self$Est_Psi0_d(U_index=1:p, obs=old_obs, quar_trunc)
         #cors0 = Est_Psi0_d(U_index=1:p, obs=old_obs, quar_trunc)
         sgn0 = rep(1, p)
         
         ### EXPENSIVE!!!
         if (p <= 1e6){
           # IF_star_mat = IF_star_self(m=sgn0, U_index=1:p, obs=old_obs, quar_trunc)
           IF_star_mat = .self$IF_star_self(m=sgn0, U_index=1:p, obs=old_obs, quar_trunc)
           utility = (cors0 + my_colMeans(IF_star_mat))
         } else { utility = cors0 }   
         
         if (num_top < p) {
           # k0 = order(abs(cors0), decreasing=TRUE)[seq(num_top)]
           k0 = order(abs(utility), decreasing=TRUE)[seq(num_top)]
         } else { k0 = seq(p) }
         
         return( k0 ) } )

       if (is.na(mt[2])) {
         save(input_data, file=paste('NAdata_stabOSE_model', model, '_n', n, '_elln', elln, '_r', r, '.Rdata', sep='') )
       }
           
       return( list(rej = rej, Sn = abs(est/se), est = est, se = se, ci = ci, p_val = p_val, k_est = k_est) )  
  },
  
  Oracle_onestep_est = function(alpha, quar_trunc, idx_taken=1){
    # Using 1:2 here is to compile the vectorized syntax in codes, even if we only need the results of the first predictor 
    cors0 = Est_Psi0_d(U_index=1:2, obs=1:n, quar_trunc)
    IF_star_mat = IF_star_self(m=rep(1,2), U_index=1:2, obs=1:n, quar_trunc)
    sigma_IF_star = colVars(IF_star_mat, sd_use=TRUE)[idx_taken]
    est = (cors0 + my_colMeans(IF_star_mat))[idx_taken]
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star) > qnorm(1-alpha/2, 0, 1, lower.tail=TRUE, log.p=FALSE) )
    p_val = 2*pnorm( sqrt(n)*abs(est/sigma_IF_star), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE )

    if (is.nan(est) | is.na(sigma_IF_star)) {
      save(input_data, file=paste('NAdata_singleOSE_model', model, '_n', n, '_r', r, '.Rdata', sep='') )
    }
    return( list( rej = rej, Sn = sqrt(n)*abs(est/sigma_IF_star), est = est, se = sigma_IF_star, p_val = p_val,
                  k_est = 1 ) )
  },
  
  Naive_onestep_est = function(alpha, quar_trunc, num_top=1){
    
    cors0 = Est_Psi0_d(U_index=1:p, obs=1:n, quar_trunc); sgn0 = rep(1, p)
    
    ### EXPENSIVE!!!
    if (p <= 1e4){
      IF_star_mat0 = IF_star_self(m=sgn0, U_index=1:p, obs=1:n, quar_trunc)
      utility = (cors0 + my_colMeans(IF_star_mat0))
    } else { utility = cors0 }   
    
    if (num_top < p) {
      k0 = order(abs(utility), decreasing=TRUE)[seq(num_top)]
    } else { k0 = seq(p) }
    
    cor_k0 = cors0[k0]
    sgn_k0 = 2*(utility[k0]>=0) - 1
    
    IF_star_mat = IF_star_self(m=sgn_k0, U_index=k0, obs=1:n, quar_trunc)
    est = (cor_k0 + mean(IF_star_mat))
    sigma_IF_star = sd(IF_star_mat)
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star) 
              > qnorm(1-alpha/2, 0, 1, lower.tail=TRUE, log.p=FALSE) )
    abs_test_stat_val = sqrt(n)*abs(est/sigma_IF_star)
    p_val = min(1, 2*pnorm(abs_test_stat_val, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE))
    
    return( list( rej = rej, Sn = abs_test_stat_val,
                  est = est, se = sigma_IF_star, p_val = p_val, k_est = k0 ) )
  },
  
  Bonf_onestep_ests = function(alpha, quar_trunc){
    cors0 = Est_Psi0_d(U_index=1:p, obs=1:n, quar_trunc)
    IF_star_mat = IF_star_self(m=rep(1,p), U_index=1:p, obs=1:n, quar_trunc)
    est = (cors0 + my_colMeans(IF_star_mat))
    sigma_IF_star = colVars(IF_star_mat, sd_use=TRUE)
    rej = 1*( sqrt(n)*max(abs(est/sigma_IF_star)) 
                            > qnorm(1-alpha/(2*p), 0, 1, lower.tail=TRUE, log.p=FALSE) )
    abs_test_stat_vals = sqrt(n)*abs(est/sigma_IF_star)
    p_vals = 2*pnorm(abs_test_stat_vals, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
    p_val = min(1, p*2*pnorm(max(abs_test_stat_vals), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE))
    
    return( list( rej = rej, Sn = sqrt(n)*max(abs(est/sigma_IF_star)),
                  est = est[which.max(abs_test_stat_vals)], se = sigma_IF_star[which.max(abs_test_stat_vals)],
                  p_val = p_val, k_est = which.max(abs_test_stat_vals), k_ests = which(p_vals <= alpha/p) ) )
  },
  
  Bonf_Cox_ests = function(alpha){ # alpha is a scalar
    data_est = cbind(input_data$X, input_data$delta, input_data$U)
    data_est = data_est[order(data_est[,1]),] %>% as.data.frame()
    colnames(data_est) = c('X', 'delta', paste0('V', 1:p))
    
    result = do.call(rbind, lapply(1:p, function(idx){
      cox_fit = coxph(Surv(X, delta) ~ data_est[,paste0('V',idx)], data_est, 
                      control = coxph.control(eps=1e-05))
      pval = summary(cox_fit)$coefficients[,5]
      return( c(idx, pval) ) }))
    
    return( list( rej = 1*(min(result[,2]) <= alpha/p),
                  p_val = min(1, p*min(result[,2])),
                  k_est = result[which.min(result[,2]), 1], k_ests = result[result[,2]<= alpha/p, 1] ) )
  }
 )
)