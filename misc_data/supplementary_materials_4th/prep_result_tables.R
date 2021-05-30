# Prepare result tables
# 1. BrainXcan/S-BrainXcan table
# 2. Genetic correlation table
# 3. Heritability and polygenicity table

# setwd('misc_data/supplementary_materials_4th/')

library(dplyr)
options(stringsAsFactors = F)
source('rlib.R')

forceall = T

outdir = 'result_tables'
dir.create(outdir)
f1 = paste0(outdir, '/', 'supp_brainxcan.csv')
f2 = paste0(outdir, '/', 'supp_gencor.csv')
f3 = paste0(outdir, '/', 'supp_genetic_arch.csv')

# IDP annotation
idp_annot = load_idp_annot()

# load prediction performance
df_perf = read.delim2('supp_table_2.tsv') %>% select(IDP, model_name, Spearman) %>% rename(model = model_name, CV_Spearman = Spearman)

# BrainXcan/S-BrainXcan
if(!file.exists(f1) | forceall) {
  k1 = readRDS('indiv_bxcan/dataframe_full.indiv_bxcan.rds')
  k2 = readRDS('s_bxcan/dataframe_full.sbxcan.rds')
  # k1 = filter_pred_perf(k1)
  # k2 = filter_pred_perf(k2)
  k = rbind(
    k1 %>% select(phenotype, IDP, model, zscore, pval), 
    k2 %>% select(phenotype, IDP, model, zscore, pval)
  )
  k = left_join(
    k, idp_annot %>% select(IDP, subtype, region, left_or_right), 
    by = 'IDP'
  ) %>% rename(side = left_or_right)
  k = left_join(k, df_perf, by = c('IDP', 'model'))
  write.csv(k, f1, row.names = F)
}
k %>% group_by(phenotype, model) %>% summarize(n = n()) %>% print(n = 100)

# Genetic correlation
if(!file.exists(f2) | forceall) {
  k3 = readRDS('s_bxcan/dataframe_full.gencor.rds')
  k3 = k3 %>% select(phenotype, IDP, rg, pval) 
  k3 = left_join(
    k3, idp_annot %>% select(IDP, subtype, region, left_or_right), 
    by = 'IDP'
  ) %>% rename(side = left_or_right)
  write.csv(k3, f2, row.names = F)
}
k3 %>% group_by(phenotype) %>% summarize(n = n()) %>% print(n = 100)

# Heritability and polygenicity
if(!file.exists(f3) | forceall) {
  # loading h2
  { 
    df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
    df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
    df1$is_pc = substr(df1$phenotype, 1, 2) == 'PC'
    df2$is_pc = substr(df2$phenotype, 1, 2) == 'PC'
    df1$pc = rep('Region-Specific', nrow(df1))
    df1$pc[df1$is_pc] = 'Common Factor'
    df2$pc = rep('Region-Specific', nrow(df2))
    df2$pc[df2$is_pc] = 'Common Factor'
    df_h2 = rbind(
      df1 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), 
      df2 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))
    )
    df_h2 = df_h2 %>% rename(IDP = phenotype)
  }
  
  pcs = df_h2 %>% filter(is_pc) %>% select(IDP, idp_type) 
  
  # loading polygenicity
  {
    idps = read.delim2('../supplementary_materials_3rd/supp_table_1.tsv', header = T) %>% mutate(IDP = paste0('IDP-', ukb_field), idp_type = t1_or_dmri) %>% select(IDP, idp_type)
    type_list = list(dMRI = 'fourth_round_dmri', T1 = 'fourth_round_t1')
    idps = rbind(
      idps, 
      pcs
    )
    dd2 = list()
    for(i in 1 : nrow(idps)) {
      fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/', type_list[[idps$idp_type[i]]], '/', idps$IDP[i], '.sld4m_all.csv')
      if(file.exists(fn)) {
        tmp = read.csv(fn)
        dd2[[length(dd2) + 1]] = tmp %>% mutate(IDP = idps$IDP[i], idp_type = idps$idp_type[i])
      }
    }
    dd2 = do.call(rbind, dd2)
    df_poly = dd2 %>% filter(Var1 == 'Manual_aggregated')
    df_poly = df_poly %>% mutate(stable = Ma_est > 0 & Ma_est > 1.96 * Ma_err)
  }
  df_both = inner_join(
    df_h2 %>% select(IDP, h2, h2_SE), 
    df_poly %>% select(IDP, Ma_est, Ma_err), 
    by = 'IDP'
  ) %>% rename(Me = Ma_est, Me_SE = Ma_err) 
  df_both = left_join(
    df_both, idp_annot %>% select(IDP, subtype, region, left_or_right), 
    by = 'IDP'
  ) %>% rename(side = left_or_right)
  write.csv(df_both, f3, row.names = F)
}