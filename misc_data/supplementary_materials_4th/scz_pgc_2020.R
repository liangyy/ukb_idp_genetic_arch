# setwd('misc_data/supplementary_materials_4th/')

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
options(stringsAsFactors = F)

outdir = 'scz_pgc_2020/'
dir.create(foldern)

kk = readRDS('result2rds/sbxcan_wPCadj.rds')
mm = kk %>% filter(phenotype == 'SCZ_PGC_2020')
mm$pval = mm$pval.univariate
is_pc = substr(mm$IDP, 1, 2) == 'PC'
mm$pval[ is_pc ] = mm$pval.adj_covar[ is_pc ]
mm %>% mutate(p_exp = rank(pval) / (n() + 1)) %>% ggplot() + 
  geom_point(aes(x = -log10(p_exp), y = -log10(pval))) + geom_abline(slope = 1)
