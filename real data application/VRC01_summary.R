##### This script is used to collect the results from the parallel computing implemented on 200 nodes.

num_r = 200
meth_vecs = c('BONF_COX', 'BONF_OSE', 'SOSE')

elln_part = c(2)
quar = 0.9

final = NULL
orig_seq = 1:num_r

for (r in orig_seq) {
        
    file = paste('VRC01_subdata_analysis_', r, '.Rdata', sep = '')
    d1 = get(load(file));
    
    if (r == 1) {
    d_SOSE = d1[(d1$quar == quar & d1$method == 'SOSE'),] 
    d_notSOSE = d1[d1$method %in% c('BONF_COX', 'BONF_OSE'),]
    
    save(d_notSOSE, file = 'results_VRC01_subdata_analysis_nonSOSE.Rdata')
    } else { d_SOSE = d1 }  
    
    d_SOSE$rep = r
    if (r %in% c(250, 500, 750, 1000)) { print(r) }
    final = rbind(final, d_SOSE)
}

save(final, file = paste('results_VRC01_subdata_analysis_SOSE', '_rep', num_r, '.Rdata', sep = ''))

