# setwd('misc_data/supplementary_materials/')

load_sbxcan = function(folder, trait_list) {
  idp_type = list(dMRI = 'dmri.original.all_covar.w_pc', T1 = 't1.scaled.all_covar.w_pc')
  models = list(ridge = 'gw_ridge', EN = 'gw_elastic_net')
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

library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib.R')

foldern = 's_bxcan'
dir.create(foldern)

plot_s_bxcan = F 
plot_s_qq = F
prep_s_mr = T

color_code = c('ridge' = 'blue', 'elastic net' = 'orange')
color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
t1_col = 'orange'  
dmri_col = 'blue'  
factor_idp = function(cc) {
  factor(cc, levels = c('T1', 'dMRI'))
}

df_gwas = read.delim2('supp_table_4.tsv', header = T)
folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd')
df_gwas$folder = 'gtex_gwas'
df_gwas$folder[31:44] = 'psychiatric'
df = list()
for(cc in c('gtex_gwas', 'psychiatric')) {
  traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
  df[[length(df) + 1]] = load_sbxcan(folders[[cc]], traits) %>% mutate(source = cc)
}
df = do.call(rbind, df)
df$model[df$model == 'EN'] = 'elastic net'
df = df %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))

idp_sig = read.table('supp_table_2.tsv', header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
df = df[df$idp_id %in% idp_sig$idp_id, ]
alpha = 0.05
min_pval = 1e-30
df = df %>% group_by(model, idp_type) %>% mutate(pval_cap = pmax(min_pval, pval)) %>% mutate(p_adj = pval * n()) %>% ungroup()

if(isTRUE(plot_s_bxcan)) {
  df_sig0 = df %>% filter(p_adj < alpha)
  df_sig1 = df %>% filter(p_adj < alpha, pip > 0.5)
  tmp0 = df %>% mutate(Total = T, Bonferroni = p_adj < alpha, Bonferroni_and_PIP = p_adj < alpha & pip > 0.5) %>% select(IDP, phenotype, model, idp_type, Bonferroni, Bonferroni_and_PIP) %>% mutate(idp_id = paste(IDP, phenotype))
  # tmp = inner_join(
  #   tmp0 %>% filter(model == 'ridge') %>% select(-model),
  #   tmp0 %>% filter(model == 'elastic net') %>% select(-model),
  #   by = c('IDP', 'phenotype', 'idp_type'),
  #   suffix = c('.ridge', '.en')
  # )
  # tmp[is.na(tmp)] = F
  # tmp_ = as.matrix(tmp[, (ncol(tmp) - 3) : ncol(tmp)])
  # tmp__ = tmp_
  # tmp_[tmp__] = 1
  # tmp_[!tmp__] = 0
  # tmp[, (ncol(tmp) - 3) : ncol(tmp)] = tmp_
  # tmp = as.data.frame(tmp)
  # upset(tmp %>% filter(idp_type == 'T1'), main.bar.color = t1_col, point.size = 5, line.size = 0.5, text.scale = 2)
  # upset(tmp %>% filter(idp_type == 'dMRI'), main.bar.color = dmri_col, point.size = 5, line.size = 0.5, text.scale = 2)
  tmp0 %>% group_by(idp_type) %>% summarize(length(unique(idp_id[Bonferroni])), length(unique(idp_id[Bonferroni_and_PIP])))
  tmp = inner_join(
    df %>% filter(model == 'ridge') %>% select(IDP, phenotype, model, idp_type, zscore, pip),
    df %>% filter(model == 'elastic net') %>% select(IDP, phenotype, model, idp_type, zscore, pip),
    by = c('IDP', 'phenotype', 'idp_type'),
    suffix = c('.ridge', '.en')
  )
  p1 = tmp %>% mutate(idp_type = factor_idp(idp_type)) %>% ggplot() + geom_point(aes(x = zscore.ridge, y = zscore.en), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 + facet_wrap(~idp_type) + 
    xlab('z-score from ridge predictors') + 
    ylab('z-score from elastic net predictors')
  p2 = tmp %>% mutate(idp_type = factor_idp(idp_type)) %>% 
    mutate(pip.ridge = pmax(1e-3, pip.ridge), pip.en = pmax(1e-3, pip.en)) %>% 
    ggplot() + geom_point(aes(x = pip.ridge, y = pip.en), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 + facet_wrap(~idp_type) + scale_x_log10() + scale_y_log10() + 
    xlab('PIP from ridge predictors') + 
    ylab('PIP from elastic net predictors')
  
  ggsave(paste0(foldern, '/s_bxcan_ridge_vs_en_zscore.png'), p1, width = 5.5, height = 3)
  ggsave(paste0(foldern, '/s_bxcan_ridge_vs_en_pip.png'), p2, width = 5.5, height = 3)
}

if(isTRUE(plot_s_qq)) {
  tmp = df %>% group_by(model) %>% mutate(pexp = rank(pval) / (n() + 1)) %>% arrange(pval) 
  p = tmp %>% 
    ggplot() + 
    geom_path(aes(x = -log10(pexp), y = -log10(pval_cap), color = model)) + 
    geom_point(data = tmp %>% filter(p_adj < 0.05), aes(x = -log10(pexp), y = -log10(pval_cap), color = model)) +
    facet_wrap(~idp_type) + 
    geom_abline(slope = 1, intercept = 0) + th2 +
    scale_color_manual(values = color_code) + 
    xlab(expression(paste(-log[10], p[expected]))) + 
    ylab(expression(paste(-log[10], p[observed]))) +
    theme(legend.position = c(0.875, 0.4), legend.title = element_blank())
  ggsave(paste0(foldern, '/s_bxcan_qqplot.png'), p, width = 5.5, height = 3)
}

if(isTRUE(prep_s_mr)) {
  df_sig = df %>% filter(p_adj < alpha, pip > 0.5)
  for(ss in unique(df_sig$source)) {
    for(i in c('T1', 'dMRI')) {
      tmp = df_sig %>% filter(idp_type == i, source == ss) %>% mutate(pp = phenotype) %>% select(phenotype, IDP, pp) %>% rename(pheno = phenotype, idp = IDP, pheno_code = pp) %>% distinct()
      write.table(tmp, paste0(foldern, '/s_bxcan_mr.', ss, '.', i, '.signif.tsv'), quo = F, col = T, row = F, sep = '\t')
    }
  }
}
