# setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/')
# update table s8
old <- read.csv('~/Downloads/Table_S8.xlsx - ..csv')
new <- readRDS("s_bxcan_permz/dataframe_full.sbxcan_permz.rds")

library(dplyr)
cols.to.update <- c('z_adj_perm_null', 'pval_adj_perm_null')
cols.join <- c('IDP', 'model', 'phenotype')
old <- left_join(old, new[, c(cols.join, cols.to.update)], by = cols.join, suffix = c('.old', ''))
write.csv(old %>% select(-z_adj_perm_null.old, -pval_adj_perm_null.old), 'Table_S8-w-factor.csv')