num_digits = 7
source('data_management.R')

# Simulate data, according to distributions described in Simulation Section of Sinica paper
Simulate_data = function(n, p, model, censoring_rate, rho, 
                         censoring_dist_par=c(0.07, 0.15, 0.28, 0.43), 
                         just.first100=FALSE, mc_true_par=FALSE, block_size=1e4){
  
  # Simulate data of predictor structure
  # Only working for p <= 1000
  # U1= rho*matrix(1,p,p) - diag(rho*rep(1,p));
  # U2= diag(rep(1,p)) + U1;
  # U = rmvnorm(n,rep(0,p),diag(1,p))%*%chol(U2);
  
  # For any p, especially p > 1000
  if (rho==0) {
    U = matrix(rnorm(n*p), ncol=p)
  } else {
    a = uniroot( function(x){ 1 - rho - (x-sqrt((1-x^2)/(p-1)))^2 }, c(0,1), 
                 tol=.Machine$double.eps )$root
    b = sqrt((1-a^2)/(p-1))
    p0 = ifelse(just.first100, 100, p)
    tmp = matrix(rnorm(n*p0), ncol=p0)
    U = matrix( rep(rowSums(tmp), p0), ncol=p0 )*b + tmp[,1:p0]*(a-b)
    rm(tmp) }
  
  # Updating p by the column number of U to accommodate the case of just.first100 = TRUE
  p0 = dim(U)[2]
  
  # Simulate data of error structure
  if (grepl('IE', model, fixed=TRUE)) {
    error_IE = rnorm(n, 0, 1)
  } else {  
    error_DE = rnorm(n, 0, 0.7*(abs(U[,1])+0.7)) } 
  
  # Simulate data of survival time based on the specified AFT model
  if(model=='N.IE'){
    log_T = error_IE
  } else if (model=='A1.IE'){
    log_T = 0.25*U[,1] + error_IE; 
  } else if (model=='A2.IE'){
    log_T = c( as.matrix(U[,1:10]) %*% matrix(c(rep(0.15,5),rep(-0.1,5)), ncol=1) ) + error_IE
  } else if(model=='N.DE'){
    log_T = error_DE
  } else if (model=='A1.DE'){
    log_T = 0.25*U[,1] + error_DE
  } else if (model=='A2.DE'){
    log_T = c( as.matrix(U[,1:10]) %*% matrix(c(rep(0.15,5),rep(-0.1,5)), ncol=1) ) + error_DE
  } else stop('Invalid model choice.')
  
  if (!(mc_true_par)){
    # Simulate data of censoring time
    if (censoring_rate=='0%'){
      C = 1e5
    }
    else{
      for (i in 1:length(censoring_dist_par)){ assign( paste('par', i, sep=''), censoring_dist_par[i] )}
      cr_par = par1*(censoring_rate=='10%') + par2*(censoring_rate=='20%') + par3*(censoring_rate=='30%') + par4*(censoring_rate=='40%')
      C = log(rexp(n, cr_par)) }
    
    X = round(pmin(log_T, C), num_digits); 
    delta = 1*(log_T <= C);
    
    # Pre-standardization and saved for further analyses
    U_colSDs = colVars(U, sd_use=TRUE, na.rm=TRUE)
    U = colStandardization(U, U_colSDs)
    # U = (U - matrix(rep(colMeans(U), n), nrow=n, byrow=TRUE))/matrix(rep(U_colSDs, n), nrow=n, byrow=TRUE)
    input_data = list(X = X, delta = delta, U = U)
  } else {
    U_colSDs = colVars(U, sd_use=TRUE, na.rm=TRUE)
    U = colStandardization(U, U_colSDs)
    # U = (U - matrix(rep(colMeans(U), n), nrow=n, byrow=TRUE))/matrix(rep(U_colSDs, n), nrow=n, byrow=TRUE)
    log_T = round(log_T, num_digits)
    input_data = list(timeT = log_T, U = U)
  }
  return( input_data )
}    
