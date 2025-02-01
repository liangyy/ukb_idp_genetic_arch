# setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/')
# update table s8
old <- read.csv('~/Downloads/Table_S8.xlsx - ..csv')
new <- readRDS("s_bxcan_twasinf_vc/dataframe_full.sbxcan_twasinf.rds")

library(dplyr)
cols.to.update <- c('pval_corrected', 'z_corrected')
cols.join <- c('IDP', 'model', 'phenotype')
old <- left_join(old, new[, c(cols.join, cols.to.update)], by = cols.join)
write.csv(old %>% select(-z_adj_perm_null, -pval_adj_perm_null), 'Table_S8-twasinf_vc.csv')