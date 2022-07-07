# table_s3
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
input_r2 <- '/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv'
df_r2 <- read.delim2(input_r2) 
write.csv(
  df_r2 %>% select(IDP, model_name, Spearman, Pearson, R2, is_kept), 
  paste0(outdir, '/table_s3.csv'), row.names = FALSE)
df_r2 %>% group_by(model_name) %>% summarize(
  min(Spearman), max(Spearman), median(Spearman), mean(Spearman > 0),
  median(R2))

mm <- read.csv(paste0(outdir, '/table_s1.csv'))
sec <- mm[substr(mm$subtype, 1, 1) == 'w', ]
sec <- unique(c(sec$IDP, sec$pc1_name))
df_r2 %>%
  mutate(is_pc = substr(IDP, 1, 2) == 'PC') %>%
  mutate(is_sec = IDP %in% sec) %>%
  group_by(model_name, IDP_type, is_pc, is_sec) %>% summarize(
  sum(Spearman > 0.1), n())
df_r2 %>%
  mutate(is_pc = substr(IDP, 1, 2) == 'PC') %>%
  mutate(is_sec = IDP %in% sec) %>%
  group_by(model_name, IDP_type, is_pc) %>% summarize(
    sum(Spearman > 0.1), n())
df_r2 %>%
  filter(substr(IDP, 1, 2) == 'PC') %>%
  filter(Spearman < 0.1)
