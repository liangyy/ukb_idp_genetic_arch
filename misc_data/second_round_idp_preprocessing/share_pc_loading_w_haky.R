# setwd('misc_data/second_round_idp_preprocessing/')


save_pc_info = function(tag, output) {
  if(file.exists(output)) {
    return(NULL)
  }
  ff = arrow::read_parquet(paste0('output/', tag, '.parquet'))
  ee = readRDS(paste0('output/', tag, '.w_pc.pca_results.rds'))
  pc_loadings = as.data.frame(ee$pc_loadings)
  colnames(pc_loadings) = paste0('PC-', 1 : ncol(pc_loadings))
  pc_loadings$IDP = colnames(ff)[-1]
  pve = c(ee$pve, 0) - c(0, ee$pve)
  pve = data.frame(
    PC = colnames(pc_loadings)[- ncol(pc_loadings)], 
    pve = pve[-length(pve)]
  )
  saveRDS(list(pc_loading = pc_loadings, pve = pve), output)
}


save_pc_info('t1.scaled.all_covar', 't1.scaled.all_covar.pc_info.rds')
save_pc_info('dmri.original.all_covar', 'dmri.original.all_covar.pc_info.rds')
