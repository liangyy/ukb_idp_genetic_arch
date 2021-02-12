# setwd('misc_data/supplementary_materials/')

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/dmri.original.all_covar.w_pc.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/t1.scaled.all_covar.w_pc.tsv.gz', header = T)
df1$is_pc = substr(df1$phenotype, 1, 2) == 'PC'
df2$is_pc = substr(df2$phenotype, 1, 2) == 'PC'


p = rbind(df1 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), df2 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))) %>%
  mutate(pheno = factor(paste(idp_type, phenotype), levels = paste(idp_type, phenotype)[order(idp_type, h2)]), PC = as.character(is_pc)) %>% 
  ggplot() + 
  geom_errorbar(aes(x = pheno, ymax = h2 + 1.96 * h2_SE, ymin = h2 - 1.96 * h2_SE), color = 'gray') + 
  geom_point(aes(x = pheno, y = h2, color = PC)) +
  theme(axis.text.x = element_blank(), legend.position = c(0.9, 0.18)) +
  scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  th2 + 
  facet_grid(.~idp_type, scales = 'free_x', space = "free_x")
ggsave('heritability.png', p, width = 7, height = 4)
               