# setwd('misc_data/supplementary_materials_4th/')
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)
df_map = read.table('~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/supplementary_materials_3rd/Elliot_2018_to_UKB_field.tsv', header = T)
df_iwas = rbind(
  read.csv('../iwas_s1_structural.csv'),
  read.csv('../iwas_s1_diffusion.csv')
) %>% inner_join(df_map, by = c('Phenotype' = 'n2018_id')) %>% mutate(IDP = paste0('IDP-', ukb_field))
df_idp = read.delim2('supp_table_1.tsv') %>% mutate(IDP = paste0('IDP-', ukb_field))
df_iwas = df_iwas %>% filter(IDP %in% df_idp$IDP)
df_bxcan = rbind(
  readRDS('s_bxcan/dataframe_full.sbxcan.rds') %>% select(IDP, phenotype, bhat, pval, zscore, model),
  readRDS('indiv_bxcan/dataframe_full.indiv_bxcan.rds') %>% select(IDP, phenotype, bhat, pval, zscore, model)
) 
ads = c('IGAP_Alzheimer', 'AD_Jansen_2019', 'parent_AD')
df_bxcan = df_bxcan %>% filter(phenotype %in% ads)
df_merge = inner_join(df_bxcan, df_iwas, by = 'IDP')

df_merge %>% ggplot() + geom_point(aes(x = zscore, y = beta_MVIWAS / se_MVIWAS, color =  Type)) + facet_grid(phenotype~model)
df_merge %>% ggplot() + geom_boxplot(aes(x = MVIWAS.Sig., y = abs(zscore))) + facet_grid(phenotype~model)
