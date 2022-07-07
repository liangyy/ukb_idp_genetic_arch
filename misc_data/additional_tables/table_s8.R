# table_s8
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- readRDS('/Users/yanyuluchicago/Documents/repo/GitHub/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/s_bxcan_permz/dataframe_full.gencor.rds')


df %>% group_by(phenotype) %>% summarize(n())
# write.csv(
#   df %>% select(IDP, phenotype, rg, pval, zscore),
#   paste0(outdir, '/table_s8.csv'), row.names = FALSE)

# misc
df2 <- read.csv(paste0(outdir, '/table_s7.csv'))
kk <- read.csv(paste0(outdir, '/table_s1.csv')) %>% select(subtype, IDP)
df2 <- left_join(df2, kk, by = 'IDP')
df <- left_join(df, kk, by = 'IDP')
df2 <- df2 %>% filter(model == 'ridge') 
is_pt <- (!is.na(df2$subtype) & substr(df2$subtype, 1, 2) == 'w-') | !is.na(stringr::str_match(df2$IDP, 'ProbTrack')[, 1] )
is_pt1 <- (!is.na(df$subtype) & substr(df$subtype, 1, 2) == 'w-') | !is.na(stringr::str_match(df$IDP, 'ProbTrack')[, 1] )
df2 <- df2 %>% filter(!is_pt)
df <- df %>% filter(!is_pt1)
mm <- inner_join(df, df2, by = c('IDP', 'phenotype'))
mm %>% group_by(phenotype) %>% 
  # filter(!is.infinite(zscore.x), !is.infinite(z_adj_perm_null)) %>%
  summarize(cor = cor(zscore.x, z_adj_perm_null)) %>% 
  summarize(min(cor), max(cor), median(cor))
sig2 <- df2 %>% group_by(phenotype) %>% filter(pval_adj_perm_null < 0.05 / n())
sig1 <- df %>% group_by(phenotype) %>% filter(pval < 0.05 / n())

# misc 2
scz <- df2 %>% filter(phenotype == 'SCZ_PGC_2020', model == 'ridge')
table(scz$subtype)
scz %>% filter(pval_adj_perm_null < 0.05 / n())
