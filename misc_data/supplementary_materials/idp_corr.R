# plot the observed phenotypic correlation before and after PC adjustment

dir.create('idp_corr')

plot_dmri_sep = function(mytmp, annot_sub, myorder) {
  # mytmp = df_cor0
  mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% left_join(myorder, by = c('Var1' = 'tag')) %>% left_join(annot_sub %>% select(idp, dmri_measure), by = c('Var1' = 'idp')) %>% group_by(dmri_measure) %>% summarize(pos = mean(order))
  mytmp %>% rename(obs_cor = value) %>% ggplot() + 
    geom_raster(aes(x = Var1, y = Var2, fill = obs_cor)) +
    annotate(geom = "text", x = mytmp_g$pos, y = -7, label = mytmp_g$dmri_measure, size = 4, vjust = 1) + 
    scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
    theme(
      plot.margin = margin(0.1, 0.1, 0.5, 0.1, "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}

plot_dmri = function(mytmp, annot_sub, myorder) {
  # mytmp = df_cor0
  mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% left_join(myorder, by = c('Var1' = 'tag')) %>% left_join(annot_sub %>% select(idp, dmri_measure, measurement_type), by = c('Var1' = 'idp')) %>% group_by(measurement_type, dmri_measure) %>% summarize(pos = mean(order))
  mytmp_g2 = mytmp[!duplicated(mytmp$Var1), ] %>% left_join(myorder, by = c('Var1' = 'tag')) %>% left_join(annot_sub %>% select(idp, measurement_type), by = c('Var1' = 'idp')) %>% group_by(measurement_type) %>% summarize(pos = max(order), posm = mean(order)) 
  mytmp %>% rename(obs_cor = value) %>% ggplot() + 
    geom_raster(aes(x = Var1, y = Var2, fill = obs_cor)) +
    geom_vline(data = mytmp_g2[-nrow(mytmp_g2), ], aes(xintercept = pos + 0.5), size = 0.5) +
    geom_hline(data = mytmp_g2[-nrow(mytmp_g2), ], aes(yintercept = pos + 0.5), size = 0.5) +
    annotate(geom = "text", x = mytmp_g$pos, y = -8, label = mytmp_g$dmri_measure, size = 4, hjust = 1, angle = 90) + 
    annotate(geom = "text", y = mytmp_g2$posm, x = -8, label = mytmp_g2$measurement_type, size = 4, hjust = 1) +
    scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1))), xlim = c(0, sum(!duplicated(mytmp$Var2)))) +
    theme(
      plot.margin = margin(0.1, 0.1, 0.5, 2.5, "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}

plot_t1 = function(mytmp, annot_sub, myorder) {
  # mytmp = df_cor0
  mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% left_join(myorder, by = c('Var1' = 'tag')) %>% left_join(annot_sub %>% select(idp, t1_anatomy_group), by = c('Var1' = 'idp')) %>% group_by(t1_anatomy_group) %>% summarize(pos = max(order), posm = mean(order)) 
  # mytmp_g = mytmp_g[-nrow(mytmp_g), ]
  mytmp %>% rename(obs_cor = value) %>% ggplot() + 
    geom_raster(aes(x = Var1, y = Var2, fill = obs_cor)) +
    annotate(geom = "text", y = mytmp_g$posm, x = -3, label = mytmp_g$t1_anatomy_group, size = 4, hjust = 1) +
    geom_vline(data = mytmp_g[-nrow(mytmp_g), ], aes(xintercept = pos + 0.5), size = 0.3) +
    geom_hline(data = mytmp_g[-nrow(mytmp_g), ], aes(yintercept = pos + 0.5), size = 0.3) +
    scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    coord_equal(expand = FALSE, clip = "off", xlim = c(0, sum(!duplicated(mytmp$Var2)))) +
    theme(
      plot.margin = margin(0.1, 0.1, 0.1, 2.5, "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}

library(dplyr)
library(ggplot2)

# meta data
annot = read.delim2('supp_table_1.tsv', sep = '\t')
annot = annot %>% mutate(idp = paste0('IDP-', ukb_field))

skip_dmri_sep = T
# dmri separate
if(skip_dmri_sep == F) {
  df0 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.no_pc.parquet')
  df1 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.w_pc.parquet')
  
  # tbss-based
  {
    annot_sub = annot %>% filter(measurement_type == 'dMRI skeleton (TBSS-style measurement)')
    annot_sub$dmri_measure = factor(annot_sub$dmri_measure, levels = c('FA', 'MD', 'MO', 'L1', 'L2', 'L3', 'OD', 'ICVF', 'ISOVF'))
    sub0 = df0[, annot_sub$idp]
    sub1 = df1[, annot_sub$idp]
    cor0 = cor(sub0)
    cor1 = cor(sub1)
    my_order = annot_sub %>% arrange(dmri_measure, anatomy, left_or_right) %>% pull(idp)
    df_order = data.frame(tag = my_order, order = 1 : length(my_order))
    
    df_cor0 = cor0 %>% reshape2::melt() 
    df_cor0$Var1 = factor(df_cor0$Var1, levels = my_order)
    df_cor0$Var2 = factor(df_cor0$Var2, levels = my_order)
    
    df_cor1 = cor1 %>% reshape2::melt() 
    df_cor1$Var1 = factor(df_cor1$Var1, levels = my_order)
    df_cor1$Var2 = factor(df_cor1$Var2, levels = my_order)
    
    p1 = plot_dmri_sep(df_cor0, annot_sub, df_order)
    p2 = plot_dmri_sep(df_cor1, annot_sub, df_order)
    
    ggsave('idp_corr/tbss_dmri_before_adj.png', p1, width = 5.5, height = 5)
    ggsave('idp_corr/tbss_dmri_after_adj.png', p2, width = 5.5, height = 5)
  }
  
  # probabilistic tractograph based
  {
    annot_sub = annot %>% filter(measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)')
    annot_sub$dmri_measure = factor(annot_sub$dmri_measure, levels = c('FA', 'MD', 'MO', 'L1', 'L2', 'L3', 'OD', 'ICVF', 'ISOVF'))
    sub0 = df0[, annot_sub$idp]
    sub1 = df1[, annot_sub$idp]
    cor0 = cor(sub0)
    cor1 = cor(sub1)
    my_order = annot_sub %>% arrange(dmri_measure, anatomy, left_or_right) %>% pull(idp)
    df_order = data.frame(tag = my_order, order = 1 : length(my_order))
    
    df_cor0 = cor0 %>% reshape2::melt() 
    df_cor0$Var1 = factor(df_cor0$Var1, levels = my_order)
    df_cor0$Var2 = factor(df_cor0$Var2, levels = my_order)
    
    df_cor1 = cor1 %>% reshape2::melt() 
    df_cor1$Var1 = factor(df_cor1$Var1, levels = my_order)
    df_cor1$Var2 = factor(df_cor1$Var2, levels = my_order)
    
    p1 = plot_dmri_sep(df_cor0, annot_sub, df_order)
    p2 = plot_dmri_sep(df_cor1, annot_sub, df_order)
    
    ggsave('idp_corr/probtract_dmri_before_adj.png', p1, width = 5.5, height = 5)
    ggsave('idp_corr/probtract_dmri_after_adj.png', p2, width = 5.5, height = 5)
  }
}

# dmri
{
  df0 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.no_pc.parquet')
  df1 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.w_pc.parquet')
  
  
  annot_sub = annot %>% filter(measurement_type %in% c('dMRI skeleton (TBSS-style measurement)', 'dMRI weighted means (probabilistic-tractography-based measurement)'))
  tmp = as.character(annot_sub$measurement_type)
  tmp[ tmp == 'dMRI skeleton (TBSS-style measurement)'] = 'TBSS'
  tmp[ tmp == 'dMRI weighted means (probabilistic-tractography-based measurement)'] = 'ProbTrack'
  annot_sub$measurement_type = factor(tmp, levels = c('TBSS', 'ProbTrack'))
  annot_sub$dmri_measure = factor(annot_sub$dmri_measure, levels = c('FA', 'MD', 'MO', 'L1', 'L2', 'L3', 'OD', 'ICVF', 'ISOVF'))
  sub0 = df0[, annot_sub$idp]
  sub1 = df1[, annot_sub$idp]
  cor0 = cor(sub0)
  cor1 = cor(sub1)
  my_order = annot_sub %>% arrange(measurement_type, dmri_measure, anatomy, left_or_right) %>% pull(idp)
  df_order = data.frame(tag = my_order, order = 1 : length(my_order))
  
  df_cor0 = cor0 %>% reshape2::melt() 
  df_cor0$Var1 = factor(df_cor0$Var1, levels = my_order)
  df_cor0$Var2 = factor(df_cor0$Var2, levels = my_order)
  
  df_cor1 = cor1 %>% reshape2::melt() 
  df_cor1$Var1 = factor(df_cor1$Var1, levels = my_order)
  df_cor1$Var2 = factor(df_cor1$Var2, levels = my_order)
  
  p1 = plot_dmri(df_cor0, annot_sub, df_order)
  p2 = plot_dmri(df_cor1, annot_sub, df_order)
  
  ggsave('idp_corr/dmri_before_adj.png', p1, width = 6, height = 5)
  ggsave('idp_corr/dmri_after_adj.png', p2, width = 6, height = 5)
  
}

# t1
{
  df0 = arrow::read_parquet('../second_round_idp_preprocessing/output/t1.scaled.all_covar.no_pc.parquet')
  df1 = arrow::read_parquet('../second_round_idp_preprocessing/output/t1.scaled.all_covar.w_pc.parquet')
  
  annot_sub = annot %>% filter(t1_or_dmri == 'T1')
  annot_sub$t1_anatomy_group = factor(annot_sub$t1_anatomy_group, levels = c('Cortical', 'Cerebellum', 'Subcortical', 'Total'))
  
  sub0 = df0[, annot_sub$idp]
  sub1 = df1[, annot_sub$idp]
  cor0 = cor(sub0)
  cor1 = cor(sub1)
  my_order = annot_sub %>% arrange(t1_anatomy_group, anatomy, left_or_right) %>% pull(idp)
  df_order = data.frame(tag = my_order, order = 1 : length(my_order))
  
  df_cor0 = cor0 %>% reshape2::melt() 
  df_cor0$Var1 = factor(df_cor0$Var1, levels = my_order)
  df_cor0$Var2 = factor(df_cor0$Var2, levels = my_order)
  
  df_cor1 = cor1 %>% reshape2::melt() 
  df_cor1$Var1 = factor(df_cor1$Var1, levels = my_order)
  df_cor1$Var2 = factor(df_cor1$Var2, levels = my_order)
  
  p1 = plot_t1(df_cor0, annot_sub, df_order)
  p2 = plot_t1(df_cor1, annot_sub, df_order)
  
  ggsave('idp_corr/t1_before_adj.png', p1, width = 5.5, height = 4)
  ggsave('idp_corr/t1_after_adj.png', p2, width = 5.5, height = 4)
  
}
