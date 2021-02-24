# setwd('misc_data/supplementary_materials/')

idp_type = list(dMRI = 'dmri.original.all_covar.w_pc', T1 = 't1.scaled.all_covar.w_pc')
models = list(ridge = 'gw_ridge', EN = 'gw_elastic_net')
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

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(dplyr)
library(patchwork)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib.R')
source('../../rmd/rlib_calc.R')

foldern = 'scz_vignette'
dir.create(foldern)

psychiatric = c('SCZ_PGC_2020', 'pgc.scz2')  


df1 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd/', psychiatric[2])
df2 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/', psychiatric[1])
df = rbind(df1, df2)
df$model[df$model == 'EN'] = 'elastic net'
df$model = factor(df$model, levels = c('ridge', 'elastic net'))
df = df %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))

idp_sig = read.table('supp_table_2.tsv', header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name), id = paste(IDP, tolower(IDP_type)))
df = df[df$idp_id %in% idp_sig$idp_id, ]
alpha = 0.05
min_pval = 1e-30

# load gen cor
{
  df_cor = list()
  tags = list(dmri = 'dmri.original.all_covar.w_pc', t1 = 't1.scaled.all_covar.w_pc')
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
  df_cor = df_cor %>% mutate(id = paste(IDP, idp_type))
  df_cor = df_cor[df_cor$id %in% idp_sig$id, ]
  df_cor$idp_type[df_cor$idp_type == 'dmri'] = 'dMRI'
  df_cor$idp_type[df_cor$idp_type == 't1'] = 'T1'
  df = rbind(
    df %>% select(IDP, phenotype, idp_type, model, zscore, pval, pip), 
    df_cor %>% mutate(model = 'genetic correlation', pip = NA) %>% select(IDP, phenotype, idp_type, model, zscore, pval, pip)
  )
  
  df = df %>% group_by(phenotype, model) %>% mutate(p_adj = pval * n()) %>% ungroup()
  
}

# load annot
{
  annot = read.delim2('supp_table_1.tsv', header = T) %>% mutate(IDP = paste0('IDP-', ukb_field))
  annot = rbind(
    annot %>% select(IDP, t1_or_dmri, anatomy, left_or_right, measurement_type, dmri_measure, t1_anatomy_group),
    data.frame(
      IDP = c(paste0('PC-', 1:9), paste0('PC-', 1:5)), 
      t1_or_dmri = c(rep('dMRI', 9), rep('T1', 5)),
      anatomy = c(paste0('PC-', 1:9), paste0('PC-', 1:5)), 
      left_or_right = NA,
      measurement_type = c(rep('dMRI PCA', 9), rep('T1 PCA', 5)),
      dmri_measure = 'PCA',
      t1_anatomy_group = 'PCA'
    )
  )
  annot[is.na(annot)] = 'NA'
  annot$measurement_type[annot$measurement_type == 'dMRI skeleton (TBSS-style measurement)'] = 'dMRI TBSS'
  annot$measurement_type[annot$measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)'] = 'dMRI ProbTrack'
  annot$measurement_type[annot$measurement_type == 'Other T1 measurements'] = 'T1 Global volumes'
  annot$measurement_type = factor(
    annot$measurement_type, 
    levels = c(
      "T1 PCA", "T1 Global volumes", "Regional grey matter volumes (FAST)", "Subcortical volumes (FIRST)",
      "dMRI PCA", "dMRI TBSS", "dMRI ProbTrack"
    )
  )
}

df = df %>% left_join(annot, by = c('IDP', 'idp_type' = 't1_or_dmri'))
df = df %>% mutate(idp_f = paste(IDP, idp_type))

# get idp order
{
  idp_order = df %>% arrange(idp_type, measurement_type, dmri_measure, t1_anatomy_group, anatomy, left_or_right) %>% select(idp_f) %>% distinct()
  order_idp = function(x) {
    factor(x, levels = idp_order$idp_f)
  }
  
}

df = df %>% mutate(idp_f = order_idp(idp_f))

not_run = F
if(isTRUE(not_run)) {
  p = df %>% filter(phenotype == psychiatric[1]) %>% ggplot() + 
    geom_point(
      data = df %>% filter(phenotype == psychiatric[1], pip > 0.5),
      aes(x = idp_f, y = -log10(pval)), size = 4, shape = 1
    ) +
    geom_point(aes(x = idp_f, y = -log10(pval), color = measurement_type)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    theme(legend.position = 'right', legend.title = element_blank()) +
    th2 +
    scale_color_manual(values = cbPalette) +
    facet_wrap(~model, ncol = 1) +
    geom_hline(
      data = df %>% filter(phenotype == psychiatric[1]) %>% 
        group_by(model) %>% summarize(cutoff = alpha / n()) %>% ungroup(),
      aes(yintercept = -log10(cutoff)), linetype = 2
    ) + 
    xlab('Brain IDPs') + 
    ylab(expression(paste(-log[10], p[S-BrainXcan]))) 
  ggsave(paste0(foldern, '/scz_overview.png'), p, width = 8, height = 5)
}

# mr preprocess
mr_prep = F
if(isTRUE(mr_prep)) {
  df_sig = df %>% filter(phenotype == psychiatric[1]) %>% filter(pip > 0.5 | (model == 'genetic correlation' & p_adj < alpha))
  # kk = df_sig %>% mutate(p2 = phenotype) %>% select(phenotype, IDP, p2)
  for(i in c('T1', 'dMRI')) {
    tmp = df_sig %>% filter(idp_type == i) %>% mutate(pp = phenotype) %>% select(phenotype, IDP, pp) %>% rename(pheno = phenotype, idp = IDP, pheno_code = pp) %>% distinct()
    write.table(tmp, paste0(foldern, '/scz2020.', i, '.signif.tsv'), quo = F, col = T, row = F, sep = '\t')
  }
}



# try visualization
no_run = T
if(!isTRUE(no_run)) {
  df %>% filter(idp_type == 'T1', phenotype == psychiatric[1]) %>% 
    ggplot() + geom_point(aes(x = idp_f, y = -log10(pval), color = model)) + facet_grid(~measurement_type, scales = 'free_x', space = 'free_x') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + th2
  p = df %>% filter(idp_type == 'dMRI', phenotype == psychiatric[1]) %>% 
    ggplot() + geom_point(aes(x = idp_f, y = -log10(pval), color = model)) + facet_wrap(~measurement_type, scales = 'free_x', ncol = 1) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + th2 
}

# en vs ridge
not_run = T
if(!isTRUE(not_run)) {
  tmp = full_join(
    df %>% filter(model == 'ridge'),
    df %>% filter(model == 'elastic net') %>% select(idp_f, phenotype, zscore, pval, pip, p_adj),
    by = c('idp_f', 'phenotype'),
    suffix = c('.ridge', '.en')
  )
  tmp2 = df %>% filter(model == 'rg') %>% select(idp_f, phenotype, zscore, pval, pip, p_adj)
  colnames(tmp2)[3:6] = paste0(colnames(tmp2)[3:6], '.rg')
  tmp = full_join(
    tmp,
    tmp2,
    by = c('idp_f', 'phenotype')
  )
  tmp$zscore.en[ is.na(tmp$zscore.en) ] = 0
  tmp$pip.en[ is.na(tmp$pip.en) ] = 0
  tmp$pip.ridge[ is.na(tmp$pip.ridge) ] = 0
  tmp$p_adj.en[ is.na(tmp$p_adj.en) ] = 1
  
  label_pip_pair = function(p1, p2, l, threshold = 0.5) {
    ll = rep('none', length(p1))
    ll[(p1 > threshold) & (p2 > threshold)] = 'both'
    ll[(p1 > threshold) & (p2 <= threshold)] = l[1]
    ll[(p1 <= threshold) & (p2 > threshold)] = l[2]
    factor(ll, levels = c('both', l[1], l[2], 'none'))
  }
  
  tmp %>% mutate(pip_label = label_pip_pair(-p_adj.en, -p_adj.ridge, c('en', 'ridge'), threshold = -0.05)) %>% 
    ggplot() + geom_point(aes(x = idp_f, y = -log10(pval.ridge), color = pip_label)) + 
    facet_wrap(~idp_type, scales = 'free_x') + theme(axis.text.x = element_blank()) + th2
  
}
