# prepare results into RDSs
# 1. S-BrainXcan
# 2. Indiv-BrainXcan
# 3. Heritability
# 4. CV performance

# setwd('misc_data/supplementary_materials_4th/')

load_sbxcan = function(folder_prefix, trait_list, idp_type, dir_suffix = NULL) {
  # idp_type = list(dMRI = 'fourth_round.dmri_no_pc', T1 = 'third_round_t1')
  models = list(ridge = 'ridge', en = 'elastic_net')
  df1 = list()
  for(t in trait_list) {
    for(i in names(idp_type)) {
      for(m in names(models)) {
        df1[[length(df1) + 1]] = read.csv(paste0(folder_prefix, '_', m, dir_suffix, '/', idp_type[[i]], '_x_', t, '_x_simagexcan.csv'), header = T) %>% mutate(idp_type = i, model = m, phenotype = t)
      }
    }
  }
  df1 = do.call(rbind, df1)
  df1
}

load_cv_perf_w_pc = function(idp_tag, model_tag, prefix = '~/Desktop/tmp/ukb_idp/idp_models_4th/fourth_round') {
  df1 = read.table(paste0(prefix, '.', idp_tag, '_no_pc', '.gw_', model_tag, '_beta.perf.tsv.gz'), header = T)
  df2 = read.table(paste0(prefix, '.', idp_tag, '_w_pc', '.gw_', model_tag, '_beta.perf.tsv.gz'), header = T)
  df = rbind(df1, df2 %>% filter(substr(phenotype, 1, 2) == 'PC'))
  df
}

trim_test = function(ss) {
  ss = as.character(ss)
  unlist(lapply(strsplit(ss, ':'), function(x) {x[1]}))
}

library(dplyr)
options(stringsAsFactors = F)
source('../../rmd/rlib.R')

sbxcan = F
ibxcan = F
herit = F
cv_perf = F
sbxcan_wPCadj = F
ibxcan_wPCadj = T

foldern = 'result2rds'
dir.create(foldern)

if(isTRUE(sbxcan)) {
  df_gwas = read.delim2('../supplementary_materials_3rd/supp_table_4.tsv', header = T)
  folder_prefixes = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_4th', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_4th')
  df_gwas$folder = 'gtex_gwas'
  df_gwas$folder[25:35] = 'psychiatric'
  idp_type = list(dMRI = 'dmri', T1 = 't1')
  df = list()
  for(cc in c('gtex_gwas', 'psychiatric')) {
    traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
    df[[length(df) + 1]] = load_sbxcan(folder_prefixes[[cc]], traits, idp_type)
  }
  df = do.call(rbind, df)
  df = df %>% mutate(zscore = p2z(pval, bhat))
  df$test = trim_test(df$test)
  df$pip = NULL
  df$cs95 = NULL
  df = left_join(
    df %>% filter(test == 'univariate') %>% select(-test),
    df %>% filter(test == 'adj_covar') %>% select(-test),
    by = c('IDP', 'idp_type', 'model', 'phenotype'),
    suffix = c('.univariate', '.adj_covar')
  )
  saveRDS(df, paste0(foldern, '/sbxcan.rds'))
}

if(isTRUE(sbxcan_wPCadj)) {
  df_gwas = read.delim2('../supplementary_materials_3rd/supp_table_4.tsv', header = T)
  folder_prefixes = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_4th', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_4th')
  df_gwas$folder = 'gtex_gwas'
  df_gwas$folder[25:35] = 'psychiatric'
  idp_type = list(dMRI = 'dmri', T1 = 't1')
  df = list()
  for(cc in c('gtex_gwas', 'psychiatric')) {
    traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
    df[[length(df) + 1]] = load_sbxcan(folder_prefixes[[cc]], traits, idp_type, dir_suffix = '_wPCadj')
  }
  df = do.call(rbind, df)
  df = df %>% mutate(zscore = p2z(pval, bhat))
  df$test = trim_test(df$test)
  df$pip = NULL
  df$cs95 = NULL
  df = left_join(
    df %>% filter(test == 'univariate') %>% select(-test),
    df %>% filter(test == 'adj_covar') %>% select(-test),
    by = c('IDP', 'idp_type', 'model', 'phenotype'),
    suffix = c('.univariate', '.adj_covar')
  )
  saveRDS(df, paste0(foldern, '/sbxcan_wPCadj.rds'))
}

if(isTRUE(ibxcan)) {
  pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes', 'bmi', 'height')
  models = list(ridge = 'ridge', en = 'elastic_net')
  idps = list(dMRI = 'dmri', T1 = 't1') 
  # types = c('linear', 'susie')
  df = list()
  for(m in names(models)) {
    for(t in names(idps)) {
      fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp1 = read.csv(fn1)
      fn2 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear_adj.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp2 = read.csv(fn2)
      tmp = rbind(tmp1, tmp2)
      df[[length(df) + 1]] = tmp %>% mutate(model = m, idp_type = t)
    }
  }
  df = do.call(rbind, df)
  df = df %>% filter(phenotype %in% pheno_interest) %>% mutate(zscore = p2z(pval, bhat))
  df$test = trim_test(df$test)
  df = left_join(
    df %>% filter(test == 'univariate') %>% select(-test),
    df %>% filter(test == 'adj_covar') %>% select(-test),
    by = c('IDP', 'idp_type', 'model', 'phenotype'),
    suffix = c('.univariate', '.adj_covar')
  )
  saveRDS(df, paste0(foldern, '/indiv_bxcan.rds'))
}

if(isTRUE(ibxcan_wPCadj)) {
  pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes', 'bmi', 'height')
  models = list(ridge = 'ridge', en = 'elastic_net')
  idps = list(dMRI = 'dmri', T1 = 't1') 
  # types = c('linear', 'susie')
  df = list()
  for(m in names(models)) {
    for(t in names(idps)) {
      fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear_residual.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp1 = read.csv(fn1)
      fn2 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear_wPCadj.fourth_round.', idps[[t]], '.gw_', models[[m]], '_beta.csv')
      tmp2 = read.csv(fn2)
      tmp = rbind(tmp1, tmp2)
      df[[length(df) + 1]] = tmp %>% mutate(model = m, idp_type = t)
    }
  }
  df = do.call(rbind, df)
  df = df %>% filter(phenotype %in% pheno_interest) %>% mutate(zscore = p2z(pval, bhat))
  df$test = trim_test(df$test)
  df = left_join(
    df %>% filter(test == 'univariate') %>% select(-test),
    df %>% filter(test == 'adj_covar') %>% select(-test),
    by = c('IDP', 'idp_type', 'model', 'phenotype'),
    suffix = c('.univariate', '.adj_covar')
  )
  saveRDS(df, paste0(foldern, '/indiv_bxcan_wPCadj.rds'))
}

if(isTRUE(herit)) {
  df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_no_pc.tsv.gz', header = T)
  df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_no_pc.tsv.gz', header = T)
  df3 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
  df4 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
  df = rbind(
    df1 %>% mutate(idp_type = 'T1'),
    df2 %>% mutate(idp_type = 'dMRI'),
    df3 %>% filter(substr(phenotype, 1, 2) == 'PC') %>% mutate(idp_type = 'T1'),
    df4 %>% filter(substr(phenotype, 1, 2) == 'PC') %>% mutate(idp_type = 'dMRI')
  )
  df = df %>% select(phenotype, idp_type, h2, h2_SE, Chisq_pval) %>% rename(IDP = phenotype, pval = Chisq_pval)
  saveRDS(df, paste0(foldern, '/heritability.rds'))
}

if(isTRUE(cv_perf)) {
  models = list(ridge = 'ridge', en = 'elastic_net')
  idps = list(dMRI = 'dmri', T1 = 't1') 
  df = list()
  for(idp in names(idps)) {
    for(model in names(models)) {
      idp_tag = idps[[idp]]
      model_tag = models[[model]]
      df[[length(df) + 1]] = load_cv_perf_w_pc(idp_tag, model_tag) %>% mutate(model = model_tag, idp_type = idp)
    }
  }
  df = do.call(rbind, df)
  df = df %>% rename(IDP = phenotype)
  saveRDS(df, paste0(foldern, '/cv_performance.rds'))
}
