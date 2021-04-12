# prepare results into one RDS
# setwd('misc_data/supplementary_materials_4th/')

load_sbxcan = function(folder, trait_list, idp_type) {
  # idp_type = list(dMRI = 'fourth_round.dmri_no_pc', T1 = 'third_round_t1')
  models = list(ridge = 'gw_ridge', EN = 'gw_elastic_net')
  df1 = list()
  for(t in trait_list) {
    for(i in names(idp_type)) {
      for(m in names(models)) {
        df1[[length(df1) + 1]] = read.csv(paste0(folder, '/', idp_type[[i]], '.', models[[m]], '_beta_x_', t, '_x_simagexcan.csv'), header = T) %>% mutate(idp_type = i, model = m, phenotype = t)
      }
    }
  }
  df1 = do.call(rbind, df1)
  df1
}

library(dplyr)
options(stringsAsFactors = F)
source('../../rmd/rlib.R')

sbxcan = F#T
ibxcan = F#T

foldern = 'bxcan2rds'
dir.create(foldern)

if(isTRUE(sbxcan)) {
  df_gwas = read.delim2('../supplementary_materials_3rd/supp_table_4.tsv', header = T)
  folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_4th', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_4th')
  df_gwas$folder = 'gtex_gwas'
  df_gwas$folder[25:35] = 'psychiatric'
  idp_type1 = list(dMRI = 'fourth_round.dmri_no_pc', T1 = 'fourth_round.t1_no_pc')
  idp_type2 = list(dMRI = 'fourth_round.dmri_w_pc', T1 = 'fourth_round.t1_w_pc')
  df1 = list()
  df2 = list()
  for(cc in c('gtex_gwas', 'psychiatric')) {
    traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
    df1[[length(df1) + 1]] = load_sbxcan(folders[[cc]], traits, idp_type1) %>% mutate(source = cc, idp_scale = 'original')
    df2[[length(df2) + 1]] = load_sbxcan(folders[[cc]], traits, idp_type2) %>% mutate(source = cc, idp_scale = 'residual')
  }
  df1 = do.call(rbind, df1)
  df2 = do.call(rbind, df2)
  df1$model[df1$model == 'EN'] = 'elastic net'
  df2$model[df2$model == 'EN'] = 'elastic net'
  df1 = df1 %>% mutate(zscore = p2z(pval, bhat))
  df2 = df2 %>% mutate(zscore = p2z(pval, bhat))
  df = rbind(df1, df2)
  saveRDS(df, paste0(foldern, '/dataframe_full.sbxcan.rds'))
}

if(isTRUE(ibxcan)) {
  pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes', 'bmi', 'height')
  models = list(ridge = 'ridge', EN = 'elastic_net')
  idps = list(dMRI = 'dmri_no_pc', T1 = 't1_no_pc') #  
  idps2 = list(dMRI = 'dmri_w_pc', T1 = 't1_w_pc')
  # types = c('linear', 'susie')
  df = list()
  for(m in names(models)) {
    for(t in names(idps)) {
      fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp1 = read.csv(fn1)
      fn2 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_susie.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp2 = read.csv(fn2)
      tmp = inner_join(tmp1, tmp2, by = c('IDP', 'phenotype'))
      df[[length(df) + 1]] = tmp %>% mutate(model = m, idp_type = t, idp_scale = 'original')
      
      fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear.fourth_round.', idps2[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp1 = read.csv(fn1)
      fn2 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_susie.fourth_round.', idps2[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp2 = read.csv(fn2)
      tmp = inner_join(tmp1, tmp2, by = c('IDP', 'phenotype'))
      df[[length(df) + 1]] = tmp %>% mutate(model = m, idp_type = t, idp_scale = 'original')
    }
  }
  df = do.call(rbind, df)
  df$bhat = - df$bhat
  df_all = df
  # df_all$model[df_all$model == 'EN'] = 'elastic net'
  df = df_all %>% filter(phenotype %in% pheno_interest) %>% mutate(zscore = p2z(pval, bhat))
  saveRDS(df, paste0(foldern, '/dataframe_full.indiv_bxcan.rds'))
  
}