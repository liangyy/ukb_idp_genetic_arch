# setwd('misc_data/third_round_idp_preprocessing/')

second_dir = '../second_round_idp_preprocessing'

library(dplyr)
library(ggplot2)
library(patchwork)
source(paste0(second_dir, '/scripts/rlib.R'))
source('rlib.R')

force_run = T
outdir = 'output'
dir.create(outdir)

# fig size
ww = 12
hh = 10
# END

tags = c('ICVF', 'ISOVF', 'OD', 'FA')

# load dmri
dmri_mat = arrow::read_parquet(paste0(second_dir, '/output/dmri.original.all_covar.parquet'))

# load idp annot
idps = read.delim2('../supplementary_materials/supp_table_1.tsv')

# out_mat
out_mat = list()

for(tag in tags) {
  outputs = c(paste0(outdir, '/', tag, '_probt_tbss.png'), paste0(outdir, '/', tag, '_probt.png'), paste0(outdir, '/', tag, '_tbss.png'))
  if(sum(file.exists(outputs)) != 3) {
    doit = T
  } else {
    doit = F
  }
  
  message('--- Working on ', tag, ' ----')
  
  if(doit | force_run) {
    # ICVF as example!
    
    # extract ICVF
    idps_icvf = idps %>% filter(t1_or_dmri == 'dMRI', dmri_measure == tag)
    
    # ICVF all
    icvf_all = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf$ukb_field)])
    kk = do_all(icvf_all)
    message('ProbTrack+TBSS PVE = ', kk$res$pve)
    
    
    # ICVF ProbTrack
    idps_icvf_prob = idps_icvf %>% filter(measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)')
    icvf_prob = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_prob$ukb_field)])
    kk2 = do_all(icvf_prob)
    message('ProbTrack PVE = ', kk2$res$pve)
    
    
    # ICVF TBSS
    idps_icvf_tbss = idps_icvf %>% filter(measurement_type == 'dMRI skeleton (TBSS-style measurement)')
    icvf_tbss = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_tbss$ukb_field)])
    kk3 = do_all(icvf_tbss)
    message('TBSS PVE = ', kk3$res$pve)
    
    ggsave(outputs[1], plot_all(kk$ps, paste0(tag, ': TBSS + ProbTrack')), width = ww, height = hh)
    ggsave(outputs[2], plot_all(kk2$ps, paste0(tag, ': ProbTrack only')), width = ww, height = hh)
    ggsave(outputs[3], plot_all(kk3$ps, paste0(tag, ': TBSS only')), width = ww, height = hh)
    
    # collect results
    mat_res2 = kk2$res$residual
    mat_res2 = as.data.frame(mat_res2)
    pc_all2 = kk2$res$pc
    pc_all2 = as.data.frame(pc_all2)
    colnames(mat_res2) = colnames(icvf_prob)
    colnames(pc_all2) = paste0('PC-', tag, '-ProbTrack-', 1 : ncol(pc_all2))
    mat_res2 = cbind(mat_res2, pc_all2)
    
    mat_res3 = kk3$res$residual
    mat_res3 = as.data.frame(mat_res3)
    pc_all3 = kk3$res$pc
    pc_all3 = as.data.frame(pc_all3)
    colnames(mat_res3) = colnames(icvf_tbss)
    colnames(pc_all3) = paste0('PC-', tag, '-TBSS-', 1 : ncol(pc_all3))
    mat_res3 = cbind(mat_res3, pc_all3)
    
    out_mat[[length(out_mat) + 1]] = cbind(mat_res2, mat_res3)
    
    saveRDS(
      list(
        TBSS = list(pc_loadings = kk3$res$pc_loadings, pve = kk3$res$pve),
        ProbTrack = list(pc_loadings = kk2$res$pc_loadings, pve = kk2$res$pve),
        Combined = list(pc_loadings = kk$res$pc_loadings, pve = kk$res$pve)
      ), 
      paste0(outdir, '/', tag, '.pca_results.rds'))
  
  }
}

out_mat = do.call(cbind, out_mat)
df_out = data.frame(individual = dmri_mat$individual)
df_out = cbind(df_out, out_mat)

message('Performing inverse normalization on residuals.')
df_out[, -1] = apply(df_out[, -1], 2, inv_norm)
message('Saving data.frame: shape = ', nrow(df_out), ' x ', ncol(df_out))
arrow::write_parquet(df_out, paste0(outdir, '/third_round_dmri.parquet'))


