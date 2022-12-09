# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '10%'; cr_num = 10
# elln 
elln_part = c(2); elln_char = 'n05'

num_r = 1000
orig_seq = 1:num_r;

final = do.call(rbind, lapply(orig_seq, function(r_idx){
  file=paste('sim_rdo_maxCorrSurv_cr', gsub("\\%", "", cr), '_', r_idx, '.Rdata', sep='')
  if (r_idx %in% c(250, 500, 750, 1000)) { print(r_idx) } 
  d1 = get(load(file)); d1$rep = r_idx
  return( d1 )
} ) )

final$est_index = est_index_val
save(final, file=paste0('rdoresults_maxCorrSurv_cr', cr_num, '_', est_index_val, '.Rdata'))
