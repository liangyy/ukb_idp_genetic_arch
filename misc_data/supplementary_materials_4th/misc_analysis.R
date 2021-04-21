# setwd('misc_data/supplementary_materials_4th/')
library(ggplot2)
library(dplyr)
source('misc_rlib.R')

# correlation structure in IDPs
{
  second_dir = '../second_round_idp_preprocessing'
  idps = read.delim2('../supplementary_materials/supp_table_1.tsv')
  df_cor = list()
  
  # dMRI
  tags = c('ICVF', 'ISOVF', 'OD', 'FA')
  dmri = arrow::read_parquet(paste0(second_dir, '/output/dmri.original.all_covar.parquet'))
  for(tag in tags) {
    idps_sub = idps %>% filter(t1_or_dmri == 'dMRI', dmri_measure == tag)
    mat_sub = as.matrix(dmri[, paste0('IDP-', idps_sub$ukb_field)])
    cc0 = get_corr(mat_sub)
    df_cor[[length(df_cor) + 1]] = cc0 %>% mutate(tag = tag)
  }
  
  # T1
  do_nothing = c('Total', 'Brainstem')
  exclude_idps = c(25901)  # it has skewed distribution
  df_tag = idps %>% 
    filter(t1_or_dmri == 'T1') %>% 
    select(t1_anatomy_group, measurement_type) %>% 
    distinct() %>% rename(tag = t1_anatomy_group) 
  t1 = arrow::read_parquet(paste0(second_dir, '/output/t1.scaled.all_covar.parquet'))
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
    if(tag %in% do_nothing) {
      next
    }
    idps_sub = idps %>% filter(t1_or_dmri == 'T1', t1_anatomy_group == df_tag$tag[i])
    if(df_tag$tag[i] == 'Subcortical') {
      idps_sub = idps_sub %>% filter(measurement_type == df_tag$measurement_type[i])
    }
    
    # ICVF all
    cols = paste0('IDP-', idps_sub$ukb_field)
    cols = cols[ ! cols %in% paste0('IDP-', exclude_idps) ]
    mat_sub = as.matrix(t1[, cols])
    cc0 = get_corr(mat_sub)
    df_cor[[length(df_cor) + 1]] = cc0 %>% mutate(tag = tag)
  }
  df_cor = do.call(rbind, df_cor)
  message('------------- Within group correlation ------------')
  df_cor %>% group_by(tag) %>% summarize(median(corr))
}