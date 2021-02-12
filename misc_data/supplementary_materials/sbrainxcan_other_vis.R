# setwd('misc_data/supplementary_materials/')



load_sbxcan = function(folder, trait_list) {
  df1 = list()
  for(t in trait_list) {
    for(i in names(idp_type)) {
      for(m in names(models)) {
        df1[[length(df1) + 1]] = read.csv(paste0(folder, '/', idp_type[[i]], '.', models[[m]], '_beta_x_', t, '_x_simagexcan.csv'), header = T) %>% mutate(idp_type = i, model = m, phenotype = t)
      }
    }
  }
  df1 = do.call(rbind, df1)
  df1
}


plot_mr = function(idp, pheno, model) {
  mr_res = readRDS(paste0('~/Desktop/tmp/ukb_idp/mr_psychiatric_2nd/MR_local.psychiatric_2nd_', model, '.', idp, '_x_', pheno, '.rds'))
  p1 = mr_res$idp2pheno$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
    geom_vline(xintercept = 0, color = 'gray') + 
    geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
    geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
    geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
    th + xlab('SNP estimated effect in IDP') + ylab('SNP estimated effect in phenotype') + ggtitle('MR: IDP -> Phenotype')
  p2 = mr_res$pheno2idp$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
    geom_vline(xintercept = 0, color = 'gray') + 
    geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
    geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
    geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
    th + xlab('SNP estimated effect in phenotype') + ylab('SNP estimated effect in IDP') + ggtitle('MR: Phenotype -> IDP')
  p1 + p2
  mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
  df_mr = rbind(
    mr_res$idp2pheno$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype'),
    mr_res$pheno2idp$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP') 
  ) %>% select(direction, method, nsnp, b, pval)
  tmp = inner_join(
    data.frame(model = c('ridge', 'EN')), 
    df %>% filter(phenotype == pheno, IDP == idp),
    by = 'model'
  )
  df_mr = rbind(
    df_mr, data.frame(direction = NA, method = c('S-BrainXcan ridge', 'S-BrainXcan EN'), nsnp = NA, b = tmp$bhat, pval = tmp$pval)
  )
  print(
    xtable::xtable(df_mr, display=c("s", "s", "s", "g", "g", 'g'), digits = 3), 
    math.style.exponents = TRUE, 
    include.rownames = FALSE, 
    file =paste0('sbrainxcan_other_vis/mr_', idp, '_x_', pheno, '.tex'), 
    sanitize.colnames.function=bold,
    booktabs = TRUE
  )
  ggsave(paste0('sbrainxcan_other_vis/mr_', idp, '_x_', pheno, '.png'), p1 + p2, width = 8, height = 4)
}


myplot = function(annot_sub, mytmp, this_idp) {
  my_order = annot_sub %>% select(anatomy, tag) %>% arrange(tag, anatomy) %>% distinct() 
  my_order = my_order %>% mutate(order = 1 : n()) # %>% select(-tag)
  mytmp_g = mytmp %>% left_join(my_order, by = c('anatomy', 'tag')) %>% group_by(tag) %>% summarize(pos = max(order) + 0.5, posm = mean(order)) %>% ungroup()
  mytmp = mytmp %>% left_join(my_order, by = c('anatomy', 'tag'))
  
  pp = mytmp %>% filter(IDP %in% this_idp)
  
  kk = unique(mytmp$left_or_right)
  kk = kk[which.max(order(kk))]
  p = mytmp %>% 
    # mutate(idp_ = factor(IDP, levels = unique(IDP[order(tag, left_or_right)]))) %>% 
    ggplot() +
    geom_point(data = pp, aes(x = order, y = -log10(pval)), shape = 1, size = 4) # +
  if('pip' %in% colnames(mytmp)) {
    p = p + geom_point(aes(x = order, y = -log10(pval), color = pip)) +
      scale_color_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0.5) 
  } else {
    p = p + geom_point(aes(x = order, y = -log10(pval)))
  }
  p = p +
    theme(axis.text.x = element_blank()) +
    facet_grid(left_or_right~.) + th2 + 
    theme(axis.ticks.x = element_blank()) +
    geom_vline(data = mytmp_g, aes(xintercept = pos)) + 
    geom_text(data = mytmp_g %>% mutate(left_or_right = kk), aes(x = posm, y = -2, label = tag), size = 3) +
    theme(legend.position = 'bottom') +
    # theme(axis.title.x = element_blank()) +
    # annotate(geom = "text", x = mytmp_g$posm, y = -2, label = mytmp_g$tag, size = 4, hjust = 1) +
    theme(
      # plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      # axis.title.y = element_blank(),
      # axis.text.y = element_blank()
    ) # + 
    # coord_equal(expand = FALSE, clip = "off", ylim = c(0, max(-log10(mytmp$pval) + 0.5)), xlim = c(0, max(mytmp$order)))
  p
}



library(ggplot2)
library(dplyr)
library(patchwork)




theme_set(theme_bw(base_size = 12))
dir.create('sbrainxcan_other_vis')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

# gtex_gwas = read.table('../../submission/simagexcan/gtex_gwas_list.txt')$V1
psychiatric = c('SCZ_PGC_2020', 'pgc.scz2')  # read.table('../preprocess_psychiatric_traits/trait_list.txt')$V1

idp_type = list(dMRI = 'dmri.original.all_covar.w_pc', T1 = 't1.scaled.all_covar.w_pc')
models = list(ridge = 'gw_ridge', EN = 'gw_elastic_net')

pheno = 'SCZ_PGC_2020'
pheno2 = 'pgc.scz2'

df1 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd', psychiatric[2])
df2 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/', psychiatric[1])
df = rbind(df1, df2)
df$model = factor(df$model, levels = c('ridge', 'EN'))

# load ldsc gen cor
df_cor = list()
tags = list(dmri = 'dmri.original.all_covar.w_pc', t1 = 't1.scaled.all_covar.w_pc')
source('../../rmd/rlib.R')
source('../../rmd/rlib_calc.R')
for(rr in psychiatric[1]) {
  for(nn in names(tags)) {
    tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_2nd/', nn, '_2nd_x_psychiatric_x_', rr, '.ldsc_rg.log')
    tmp = load_ldsc_rg(tmp)
    tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
    df_cor[[length(df_cor) + 1]] = tmp %>% 
      mutate(idp_type = nn, phenotype = rr)
  }
}
for(rr in psychiatric[2]) {
  for(nn in names(tags)) {
    tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_2nd/', nn, '_2nd_x_gtex-gwas_x_', rr, '.ldsc_rg.log')
    tmp = load_ldsc_rg(tmp)
    tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
    df_cor[[length(df_cor) + 1]] = tmp %>% 
      mutate(idp_type = nn, phenotype = rr)
  }
}
df_cor = do.call(rbind, df_cor)


annot = read.delim2('supp_table_1.tsv', header = T) %>% mutate(IDP = paste0('IDP-', ukb_field))



df = df %>% inner_join(annot, by = 'IDP')
df_cor = df_cor %>% inner_join(annot, by = 'IDP')


annot_sub_t1 = annot %>% filter(t1_or_dmri == 'T1', t1_anatomy_group != 'Total', !is.na(left_or_right))
annot_sub_t1$t1_anatomy_group = factor(annot_sub_t1$t1_anatomy_group, levels = c('Cortical', 'Cerebellum', 'Subcortical', 'Total'))
annot_sub_t1$tag = annot_sub_t1$t1_anatomy_group


mytmp_t1 = df %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno, model == 'ridge', tag != 'Total', !is.na(left_or_right))
mytmp_t1_cor = df_cor %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno, tag != 'Total', !is.na(left_or_right))
mytmp2_t1 = df %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno2, model == 'ridge', tag != 'Total', !is.na(left_or_right))
mytmp2_t1_cor = df_cor %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno2, tag != 'Total', !is.na(left_or_right))


annot_sub_dmri = annot %>% filter(t1_or_dmri == 'dMRI', measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
annot_sub_dmri$dmri_measure = factor(annot_sub_dmri$dmri_measure, levels = c('FA', 'MD', 'MO', 'L1', 'L2', 'L3', 'OD', 'ICVF', 'ISOVF'))
annot_sub_dmri$tag = annot_sub_dmri$dmri_measure


mytmp_dmri = df %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno, model == 'ridge', measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
mytmp_dmri_cor = df_cor %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno, measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
mytmp2_dmri = df %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno2, model == 'ridge', measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
mytmp2_dmri_cor = df_cor %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno2, measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))

p_t1 = myplot(annot_sub_t1, mytmp_t1, this_idp = 'IDP-25020'); p_t1
p_t1_c = myplot(annot_sub_t1, mytmp_t1_cor, this_idp = 'IDP-25020'); p_t1_c
p_dmri = myplot(annot_sub_dmri, mytmp_dmri, this_idp = c('IDP-25361', 'IDP-25430')); p_dmri
p_dmri_c = myplot(annot_sub_dmri, mytmp_dmri_cor, this_idp = c('IDP-25361', 'IDP-25430')); p_dmri_c

p2_t1 = myplot(annot_sub_t1, mytmp2_t1, this_idp = 'IDP-25020'); p2_t1
p2_t1_c = myplot(annot_sub_t1, mytmp2_t1_cor, this_idp = 'IDP-25020'); p2_t1_c
p2_dmri = myplot(annot_sub_dmri, mytmp2_dmri, this_idp = c('IDP-25361', 'IDP-25430')); p2_dmri
p2_dmri_c = myplot(annot_sub_dmri, mytmp2_dmri_cor, this_idp = c('IDP-25361', 'IDP-25430')); p2_dmri_c

mytmp_t1_e = df %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno, model == 'EN', tag != 'Total', !is.na(left_or_right))
mytmp2_t1_e = df %>% mutate(tag = t1_anatomy_group) %>% filter(t1_or_dmri == 'T1', phenotype == pheno2, model == 'EN', tag != 'Total', !is.na(left_or_right))
p_t1_e = myplot(annot_sub_t1, mytmp_t1_e, this_idp = 'IDP-25020'); p_t1_e
p2_t1_e = myplot(annot_sub_t1, mytmp2_t1_e, this_idp = 'IDP-25020'); p2_t1_e

mytmp_dmri_e = df %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno, model == 'EN', measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
mytmp2_dmri_e = df %>% mutate(tag = dmri_measure) %>% filter(t1_or_dmri == 'dMRI', phenotype == pheno2, model == 'EN', measurement_type != 'dMRI weighted means (probabilistic-tractography-based measurement)') %>% mutate(left_or_right = as.character(left_or_right))
p_dmri_e = myplot(annot_sub_dmri, mytmp_dmri_e, this_idp = c('IDP-25361', 'IDP-25430')); p_dmri_e
p2_dmri_e = myplot(annot_sub_dmri, mytmp2_dmri_e, this_idp = c('IDP-25361', 'IDP-25430')); p2_dmri_e


ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_ridge.png'), p_t1, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_en.png'), p_t1_e, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_ridge.png'), p_dmri, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_en.png'), p_dmri_e, height = 4, width = 8)

ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_gencor.png'), p_t1_c, height = 3.5, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_gencor.png'), p_dmri_c, height = 3.5, width = 8)

ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_ridge_2.png'), p2_t1, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_en_2.png'), p2_t1_e, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_ridge_2.png'), p2_dmri, height = 4, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_en_2.png'), p2_dmri_e, height = 4, width = 8)

ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_t1_gencor_2.png'), p2_t1_c, height = 3.5, width = 8)
ggsave(paste0('sbrainxcan_other_vis/', pheno, '_x_dmri_gencor_2.png'), p2_dmri_c, height = 3.5, width = 8)

model = 't1'
idp1 = 'IDP-25020'
plot_mr(idp = idp1, pheno = pheno, model = model)


model = 'dmri'
idp2 = 'IDP-25430'
plot_mr(idp = idp2, pheno = pheno, model = model)

model = 'dmri'
idp = 'IDP-25361'
plot_mr(idp = idp, pheno = pheno, model = model)


annot %>% filter(IDP %in% c(idp1, idp2))


mm = df %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'pip') %>% rename(pgc2 = pgc.scz2, pgc3 = SCZ_PGC_2020)
p1 = mm %>% ggplot() +
  geom_point(aes(x = pgc2, y = pgc3), alpha = 0.5) + th2 + 
  facet_wrap(~model) + coord_equal() + geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) + 
  geom_point(data = mm %>% filter(IDP %in% c('IDP-25020', 'IDP-25361', 'IDP-25430')), aes(x = pgc2, y = pgc3), shape = 1, size = 4)
mm = df %>% mutate(zscore = p2z(p = pval, b = bhat)) %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'zscore') %>% rename(pgc2 = pgc.scz2, pgc3 = SCZ_PGC_2020) 
p2 = mm %>% ggplot() +
  geom_point(aes(x = pgc2, y = pgc3), alpha = 0.5) + th2 + 
  facet_wrap(~model) + coord_equal() + geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # theme(legend.position = c(0.6, 0.8), legend.title = element_blank()) + 
  geom_point(data = mm %>% filter(IDP %in% c('IDP-25020', 'IDP-25361', 'IDP-25430')), aes(x = pgc2, y = pgc3), shape = 1, size = 4)
ggsave('sbrainxcan_other_vis/pgc2_vs_3_pip.png', p1, height = 3, width = 5)
ggsave('sbrainxcan_other_vis/pgc2_vs_3_zscore.png', p2, height = 3.5, width = 5)
