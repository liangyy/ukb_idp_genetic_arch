# setwd('misc_data/fourth_round_idp_preprocessing/')

second_dir = '~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/second_round_idp_preprocessing/'
third_dir = '~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/third_round_idp_preprocessing/'
library(dplyr)
library(ggplot2)
library(patchwork)
source(paste0(second_dir, '/scripts/rlib.R'))
# source(paste0(third_dir, '/rlib.R'))
source('rlib.R')
source('rlib_imac.R')

force_run = T
save_yaml = T
save_df = F
save_small = T
outdir = 'output'
dir.create(outdir)

# fig size
ww = 12
hh = 10
# END

tags = c('ICVF', 'ISOVF', 'OD', 'FA')

# load dmri
dmri_mat = read_parquet(paste0(second_dir, '/output/dmri.original.all_covar.parquet'))

# load idp annot
idps = read.delim2(paste0(second_dir, '../supplementary_materials/supp_table_1.tsv'))

# out_mat
out_mat = list()

# yaml
idp_dict = list()

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
    
    
    # ICVF ProbTrack
    idps_icvf_prob = idps_icvf %>% filter(measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)')
    icvf_prob = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_prob$ukb_field)])
    kk2 = do_all(icvf_prob)
    
    # ICVF TBSS
    idps_icvf_tbss = idps_icvf %>% filter(measurement_type == 'dMRI skeleton (TBSS-style measurement)')
    icvf_tbss = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_tbss$ukb_field)])
    kk3 = do_all(icvf_tbss)
    
    if(isTRUE(save_df)) {
      ggsave(outputs[1], plot_all(kk$ps, paste0(tag, ': TBSS + ProbTrack')), width = ww, height = hh)
      ggsave(outputs[2], plot_all(kk2$ps, paste0(tag, ': ProbTrack only')), width = ww, height = hh)
      ggsave(outputs[3], plot_all(kk3$ps, paste0(tag, ': TBSS only')), width = ww, height = hh)
    }
    
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
    
    idp_dict[[paste0(tag, '-ProbTrack')]] = list(
      covariates = colnames(pc_all2),
      x = colnames(icvf_prob)
    )
    idp_dict[[paste0(tag, '-TBSS')]] = list(
      covariates = colnames(pc_all3),
      x = colnames(icvf_tbss)
    )
    
    if(isTRUE(save_df) | isTRUE(save_small)) {
      saveRDS(
        list(
          TBSS = list(
            pc_loadings = data.frame(pc_loading = kk3$res$pc_loadings, IDP = colnames(icvf_tbss)), 
            pve = kk3$res$pve
          ),
          ProbTrack = list(
            pc_loadings = data.frame(kk2$res$pc_loadings, IDP = colnames(icvf_prob)), 
            pve = kk2$res$pve
          ),
          Combined = list(
            pc_loadings = kk$res$pc_loadings, 
            pve = kk$res$pve
          )
        ), 
        paste0(outdir, '/', tag, '.pca_results.rds'))
    }
    
  }
}

if(isTRUE(save_yaml)) {
  yaml::write_yaml(idp_dict, paste0(outdir, '/dmri_covar.yaml'))
  kk = idp_dict
  for(n in names(kk)) {
    tmp_x = kk[[n]][['x']]
    kk[[n]][['x']] = kk[[n]][['covariates']]
    kk[[n]][['covariates']] = tmp_x
  }
  yaml::write_yaml(kk, paste0(outdir, '/dmri_covar_for_pc.yaml'))
}

if(doit | force_run) {
  out_mat = do.call(cbind, out_mat)
  df_out = data.frame(individual = dmri_mat$individual)
  df_out = cbind(df_out, out_mat)

  message('Performing inverse normalization on residuals.')
  df_out[, -1] = apply(df_out[, -1], 2, inv_norm)
  message('Saving data.frame: shape = ', nrow(df_out), ' x ', ncol(df_out))
  if(isTRUE(save_df)) {
    write_parquet(df_out, paste0(outdir, '/fourth_round.dmri_w_pc.parquet'))
  }
  
  message('Performing inverse normalization on originals.')
  idps_sub = idps %>% filter(t1_or_dmri == 'dMRI', dmri_measure %in% tags)
  dmri_mat_sub = dmri_mat[, c('individual', paste0('IDP-', idps_sub$ukb_field))]
  dmri_mat_sub[, -1] = apply(dmri_mat_sub[, -1], 2, inv_norm)
  message('Saving data.frame: shape = ', nrow(dmri_mat_sub), ' x ', ncol(dmri_mat_sub))
  if(isTRUE(save_df)) {
    write_parquet(dmri_mat_sub, paste0(outdir, '/fourth_round.dmri_no_pc.parquet'))
  }
}

