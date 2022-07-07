# table_s6
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- readRDS('/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/indiv_bxcan/dataframe_full.indiv_bxcan.rds')
df <- df %>% mutate(idp_id = paste(idp_type, IDP, model))


idp_sig <- read.table('/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv', header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
df <- df[df$idp_id %in% idp_sig$idp_id, ]
df %>% group_by(phenotype, model) %>% summarize(n())
write.csv(
  df %>% select(model, IDP, phenotype, bhat, pval, zscore), 
  paste0(outdir, '/table_s6.csv'), row.names = FALSE)