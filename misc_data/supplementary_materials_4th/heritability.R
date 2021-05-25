# setwd('misc_data/supplementary_materials_4th/')

gen_panel = function(tmp, tag, idp_type = 'T1') {
  tmp_ = tmp %>% filter(t1_or_dmri == idp_type, subtype %in% tag)
  tmp_$subtype = factor(tmp_$subtype, levels = tags)
  p1 = tmp_ %>% 
    mutate(pheno = factor(paste(idp_type, phenotype), levels = paste(idp_type, phenotype)[order(idp_type, h2)]), PC = as.character(is_pc)) %>% 
    ggplot() + 
    geom_errorbar(aes(x = pheno, ymax = h2 + 1.96 * h2_SE, ymin = h2 - 1.96 * h2_SE), color = 'gray') + 
    geom_point(aes(x = pheno, y = h2, color = pc)) +
    theme(axis.text.x = element_blank(), legend.position = 'none', legend.title = element_blank()) +
    scale_color_manual(values = c('Common Factor' = 'red', 'Region-Specific' = 'black')) +
    th2 + 
    xlab('Brain IDPs ordered by heritability') +
    ylab('Heritability') +
    facet_grid(.~subtype, scales = 'free_x', space = "free_x") +
    theme(axis.ticks.x = element_blank()) + 
    coord_cartesian(ylim = c(0, 0.45))
  p1
}

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
library(patchwork)
source('rlib.R')

outdir = 'heritability'
dir.create(outdir)

df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
df1$is_pc = substr(df1$phenotype, 1, 2) == 'PC'
df2$is_pc = substr(df2$phenotype, 1, 2) == 'PC'
df1$pc = rep('Region-Specific', nrow(df1))
df1$pc[df1$is_pc] = 'Common Factor'
df2$pc = rep('Region-Specific', nrow(df2))
df2$pc[df2$is_pc] = 'Common Factor'

idp_annot = load_idp_annot()

tmp = rbind(
  df1 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), 
  df2 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))
) %>% left_join(idp_annot, by = c('phenotype' = 'IDP')) 
tmp$subtype = as.character(tmp$subtype)
other_subtypes = c('Brainstem', 'Total', 'Gray-Brainstem')
tmp$subtype[ tmp$subtype %in% other_subtypes ] = 'Other'


{
  tags = c('Gray-Subcortical', 'Subcortical', 'Gray-Cerebellum', 'Other')
  p1 = gen_panel(tmp, tags) + theme(axis.title.x = element_blank())
  p1
  
  tags = 'Gray-Cortical'
  p2 = gen_panel(tmp, tags) + theme(legend.position = c(0.2, 0.65)) + theme(axis.title.x = element_blank()) + theme(legend.background=element_rect(fill = alpha("white", 0)))
  
  # ggsave(paste0(outdir, '/heritability', '.t1', '.png'), (p1 / p2), width = 7, height = 5)

  tags = c('FA', 'ICVF', 'ISOVF', 'OD')
  p3 = gen_panel(tmp, tags, idp_type = 'dMRI')
  p3
  
  
  ggsave(paste0(outdir, '/heritability.png'), (p1 / p2 / p3), width = 7, height = 5)
  
  # another layout 
  p4 = p3 + facet_wrap(.~subtype, ncol = 2, scales = 'free_x')
  ggsave(paste0(outdir, '/heritability2.png'), ((p1 / p2) | p4)  + plot_layout(ncol = 2, widths = c(1, 1)), width = 11, height = 5)
  
}

df = remove_probtrack_idp(tmp %>% mutate(IDP = phenotype))
df %>% summarize(min = min(h2), max = max(h2))
