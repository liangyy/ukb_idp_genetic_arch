source('../fourth_round_idp_preprocessing/rlib_imac.R')
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
library(patchwork)

outdir = 'plot_preprocess'
dir.create(outdir)

ff = readRDS('../fourth_round_idp_preprocessing/output/ICVF.pca_results.rds')
pc_loadings = ff$TBSS$pc_loadings
pc_loadings$IDP = as.character(pc_loadings$IDP)

d1 = read_parquet('~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/fourth_round_idp_preprocessing/output/fourth_round.dmri_w_pc.parquet')[, pc_loadings$IDP]
d0 = read_parquet('~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/fourth_round_idp_preprocessing/output/fourth_round.dmri_no_pc.parquet')[, pc_loadings$IDP]
cor0 = reshape2::melt(cor(d0)) 
cor1 = reshape2::melt(cor(d1))
tmp = rbind(
  cor0 %>% mutate(type = 'Original'), 
  cor1 %>% mutate(type = 'After removing PC1')
) %>% mutate(type = factor(type, levels = c('Original', 'After removing PC1')))
tmp$Var1 = as.character(tmp$Var1)
tmp$Var2 = as.character(tmp$Var2)

p1 = tmp %>%
  ggplot() + 
  geom_raster(aes(x = Var1, y = Var2, fill = value)) + 
  scale_fill_gradient2(low = 'blue', high = 'red') + 
  facet_wrap(~type) + 
  coord_equal() + th2 + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + theme(legend.title = element_blank()); p1
p2 = tmp %>% filter(Var1 > Var2) %>% ggplot() + geom_histogram(aes(x = value, fill = type), position = 'dodge', bins = 30) + theme(legend.position = c(0.85, 0.9), legend.title = element_blank()) + th + xlab('Correlation between IDPs') + ylab('Occurrence'); p2
p3 = pc_loadings %>% ggplot() + geom_point(aes(x = IDP, y = pc_loading)) + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab('ICVF for brain regions') + 
  ylab('IDP loadings on PC 1'); p3
ggsave(paste0(outdir, '/ICVF_1.png'), p1, width = 6, height = 3)
ggsave(paste0(outdir, '/ICVF_2.png'), p3 / p2, width = 6, height = 4)
