source('../fourth_round_idp_preprocessing/rlib_imac.R')
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
library(patchwork)

run_dmri = T
run_t1 = T

outdir = 'plot_preprocess'
dir.create(outdir)
dmri_modalities = c('FA', 'ICVF', 'ISOVF', 'OD')
t1_modalities = c('Subcortical_vol', 'Subcortical_GMvol', 'Cortical', 'Cerebellum')

if(isTRUE(run_dmri)) {
  for(dmri in dmri_modalities) {
    ff = readRDS(paste0('../fourth_round_idp_preprocessing/output/', dmri, '.pca_results.rds'))
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
    p2 = tmp %>% filter(Var1 > Var2) %>% ggplot() + geom_histogram(aes(x = value, fill = type), position = 'dodge', bins = 30) + theme(legend.position = 'right', legend.title = element_blank()) + th + xlab('Correlation between IDPs') + ylab('Occurrence'); p2
    p3 = pc_loadings %>% ggplot() + geom_point(aes(x = IDP, y = pc_loading)) + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      xlab(paste0(dmri, ' for brain regions')) + 
      ylab('IDP loadings on PC 1'); p3
    ggsave(paste0(outdir, '/', dmri, '_1.png'), p1 / p2, width = 6, height = 4.5)
    ggsave(paste0(outdir, '/', dmri, '_2.png'), p3, width = 4, height = 2)
  }
}

if(isTRUE(run_t1)) {
  for(t1 in t1_modalities) {
    ff = readRDS(paste0('../fourth_round_idp_preprocessing/output/', t1, '.pca_results.rds'))
    pc_loadings = ff$pc_loadings
    pc_loadings$IDP = as.character(pc_loadings$IDP)
    
    d1 = read_parquet('~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/fourth_round_idp_preprocessing/output/fourth_round.t1_w_pc.parquet')[, pc_loadings$IDP]
    d0 = read_parquet('~/Desktop/tmp/backup_from_macbook/ukb_idp_genetic_arch/misc_data/fourth_round_idp_preprocessing/output/fourth_round.t1_no_pc.parquet')[, pc_loadings$IDP]
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
    p2 = tmp %>% filter(Var1 > Var2) %>% ggplot() + geom_histogram(aes(x = value, fill = type), position = 'dodge', bins = 30) + theme(legend.position = 'right', legend.title = element_blank()) + th + xlab('Correlation between IDPs') + ylab('Occurrence'); p2
    p3 = pc_loadings %>% ggplot() + geom_point(aes(x = IDP, y = pc_loading)) + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      xlab(paste0(pc_map(t1), ' for brain regions')) + 
      ylab('IDP loadings on PC 1'); p3
    ggsave(paste0(outdir, '/', t1, '_1.png'), p1 / p2, width = 6, height = 4.5)
    ggsave(paste0(outdir, '/', t1, '_2.png'), p3, width = 4, height = 2)
  }
}
