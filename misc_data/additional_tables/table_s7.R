# table_s7
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- readRDS('/Users/yanyuluchicago/Documents/repo/GitHub/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/s_bxcan_permz/dataframe_full.sbxcan_permz.rds')
df <- df %>% mutate(idp_id = paste(modality, IDP, model))


idp_sig <- read.table('/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv', header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
df <- df[df$idp_id %in% idp_sig$idp_id, ]
df %>% group_by(phenotype, model) %>% summarize(n())
write.csv(
  df %>% select(model, IDP, phenotype, bhat, pval, zscore, nsnp_used, nsnp_total, z_adj_perm_null, pval_adj_perm_null),
  paste0(outdir, '/table_s7.csv'), row.names = FALSE)

# limit to ridge and no probtrack
df <- df %>% filter(model == 'ridge') 
is_pt <- substr(df$subtype, 1, 2) == 'w-' | !is.na(stringr::str_match(df$IDP, 'ProbTrack')[, 1] )
df <- df %>% filter(!is_pt)
sig <- df %>% group_by(phenotype) %>% filter(pval_adj_perm_null < 0.05 / n())
ns <- length(unique(sig$IDP))
nt <- length(unique(df$IDP))
message(ns, ' out of ', nt, ' are Bonferroni significant alpha = 0.05', ' fraction = ', ns / nt)
