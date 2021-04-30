# setwd('misc_data/fourth_round_idp_preprocessing/')

second_dir = '../second_round_idp_preprocessing'
third_dir = '../third_round_idp_preprocessing'
library(dplyr)
library(ggplot2)
library(patchwork)
source(paste0(second_dir, '/scripts/rlib.R'))
# source(paste0(third_dir, '/rlib.R'))
source('rlib.R')


force_run = T
save_yaml = T
save_df = F
outdir = 'output'
dir.create(outdir)

# fig size
ww = 12
hh = 10
# END

do_nothing = c('Total', 'Brainstem')
exclude_idps = c(25901)  # it has skewed distribution

# load t1
t1_mat = arrow::read_parquet(paste0(second_dir, '/output/t1.scaled.all_covar.parquet'))

# load idp annot
idps = read.delim2('../supplementary_materials/supp_table_1.tsv')


df_tag = idps %>% 
  filter(t1_or_dmri == 'T1') %>% 
  select(t1_anatomy_group, measurement_type) %>% 
  distinct() %>% rename(tag = t1_anatomy_group) 



# out_mat
out_mat = list()


# yaml
idp_dict = list()

for(i in 1 : nrow(df_tag)) {
  
  tag = df_tag$tag[i]
  if(tag == 'Brainstem' & df_tag$measurement_type[i] == 'Regional grey matter volumes (FAST)') {
    next
  }
  if(tag == 'Subcortical') {
    if(df_tag$measurement_type[i] == 'Subcortical volumes (FIRST)') {
      subtag = 'vol'
    } else if (df_tag$measurement_type[i] == 'Regional grey matter volumes (FAST)') {
      subtag = 'GMvol'
    } else {
      subtag = NA
    }
    tag = paste0(tag, '_', subtag)
  }
  
  outputs = paste0(outdir, '/', tag, '_t1.png')
  if(sum(file.exists(outputs)) != 1) {
    doit = T
  } else {
    doit = F
  }
  if(doit | force_run) {
    message('--- Working on ', tag, ' ----')
    
    # cortical!
    
    # extract ICVF
    idps_cort = idps %>% filter(t1_or_dmri == 'T1', t1_anatomy_group == df_tag$tag[i])
    if(df_tag$tag[i] == 'Subcortical') {
      idps_cort = idps_cort %>% filter(measurement_type == df_tag$measurement_type[i])
    }
    
    # ICVF all
    cols = paste0('IDP-', idps_cort$ukb_field)
    cols = cols[ ! cols %in% paste0('IDP-', exclude_idps) ]
    cort_all = as.matrix(t1_mat[, cols])
    kk = do_all(cort_all, skip_pca = tag %in% do_nothing)
    
    if(isTRUE(save_df)) {
      ggsave(outputs[1], plot_all(kk$ps, tag), width = ww, height = hh)
    }
    
    # collect results
    mat_res = kk$res$residual
    mat_res = as.data.frame(mat_res)
    
    colnames(mat_res) = colnames(cort_all)
    pc_all = kk$res$pc
    if(!is.null(pc_all)) {
      pc_all = as.data.frame(pc_all)
      colnames(pc_all) = paste0('PC-', tag, '-', 1 : ncol(pc_all))
      mat_res = cbind(mat_res, pc_all)
      idp_dict[[as.character(tag)]] = list(
        covariates = colnames(pc_all),
        x = colnames(cort_all)
      )
    }
    
    out_mat[[length(out_mat) + 1]] = mat_res
    if(isTRUE(save_df)) {
      saveRDS(
        list(pc_loadings = kk$res$pc_loadings, pve = kk$res$pve), 
        paste0(outdir, '/', tag, '.pca_results.rds')
      )
    }
  }
}

if(isTRUE(save_yaml)) {
  yaml::write_yaml(idp_dict, paste0(outdir, '/t1_covar.yaml'))
  kk = idp_dict
  for(n in names(kk)) {
    tmp_x = kk[[n]][['x']]
    kk[[n]][['x']] = kk[[n]][['covariates']]
    kk[[n]][['covariates']] = tmp_x
  }
  yaml::write_yaml(kk, paste0(outdir, '/t1_covar_for_pc.yaml'))
}

if(doit | force_run) {
  # add cols being excluded
  dd = as.matrix(t1_mat[, paste0('IDP-', exclude_idps) ])
  dd = standardize(dd)
  mat_res = as.data.frame(dd)
  colnames(mat_res) = colnames(dd)
  out_mat = cbind(out_mat, mat_res)

  out_mat = do.call(cbind, out_mat)
  df_out = data.frame(individual = t1_mat$individual)
  df_out = cbind(df_out, out_mat)

  message('Performing inverse normalization on residuals.')
  df_out[, -1] = apply(df_out[, -1], 2, inv_norm)
  message('Saving data.frame: shape = ', nrow(df_out), ' x ', ncol(df_out))
  if(isTRUE(save_df)) {
    arrow::write_parquet(df_out, paste0(outdir, '/fourth_round.t1_w_pc.parquet'))
  }

  message('Performing inverse normalization on originals.')
  t1_mat[, -1] = apply(t1_mat[, -1], 2, inv_norm)
  message('Saving data.frame: shape = ', nrow(t1_mat), ' x ', ncol(t1_mat))
  if(isTRUE(save_df)) {
    arrow::write_parquet(t1_mat, paste0(outdir, '/fourth_round.t1_no_pc.parquet'))
  }
}