# setwd('misc_data/supplementary_materials/')

# set to TRUE if want to replot the qq-plot and pip
overview_plot = F
example_plot = F
overview_plot2 = T
gencor_plot = T

library(ggplot2)
library(dplyr)
library(patchwork)
theme_set(theme_bw(base_size = 15))
dir.create('sbrainxcan_run')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib_calc.R')

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

gtex_gwas = read.table('../../submission/simagexcan/gtex_gwas_list.txt')$V1
psychiatric = read.table('../preprocess_psychiatric_traits/trait_list.txt')$V1

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


df1 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd', gtex_gwas)
df2 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/', psychiatric)
df = rbind(df1, df2)
df$model = factor(df$model, levels = c('ridge', 'EN'))
df = df %>% mutate(zscore = p2z(pval, bhat))

if(isTRUE(overview_plot)) {
  p = df %>% group_by(idp_type, model) %>% mutate(p_exp = rank(pval) / (n() + 1)) %>% ungroup() %>% 
    mutate(pval_ = pmax(pval, 1e-20)) %>% 
    ggplot() + geom_point(aes(x = -log10(p_exp), y = -log10(pval_), color = idp_type)) + 
    xlab(expression(paste(-log[10], p[expected]))) + 
    ylab(expression(paste(-log[10], p[observed]))) +
    # scale_shape_manual(values = c('ridge' = 1, 'EN' = 3)) +
    scale_color_manual(values = c('dMRI' = 'orange', 'T1' = 'blue')) +
    geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 + 
    theme(legend.position = c(0.63, 0.7)) +
    facet_wrap(~model)
  ggsave('sbrainxcan_run/qqplot.png', p, width = 6, height = 3)
  
  p = df %>% ggplot() + geom_jitter(aes(x = phenotype, y = pip, color = idp_type), height = 0, alpha = 0.5) + 
    theme(axis.text.x = element_blank(), legend.position = 'bottom') + facet_wrap(~model) +
    scale_color_manual(values = c('dMRI' = 'orange', 'T1' = 'blue')) + th2 
  ggsave('sbrainxcan_run/pip.png', p, width = 6, height = 3)
  
}

if(isTRUE(overview_plot2)) {
  df_sig = df %>% group_by(idp_type, model) %>% summarize(p_adj_cutoff = 0.05 / n()) %>% ungroup() %>% mutate(z = p2z(p_adj_cutoff, 1))
  tmp = df %>% reshape2::dcast(IDP + phenotype + idp_type ~ model, value.var = 'zscore')
  tmp = tmp %>% left_join(df_sig %>% filter(model == 'ridge') %>% select(idp_type, z), by = 'idp_type') 
  tmp$ridge_sig = tmp$z < abs(tmp$ridge) 
  tmp = tmp %>% select(-z)
  tmp = tmp %>% left_join(df_sig %>% filter(model == 'EN') %>% select(idp_type, z), by = 'idp_type') 
  tmp$en_sig = tmp$z < abs(tmp$EN) 
  tmp = tmp %>% select(-z)
  tmp$label = apply(tmp %>% select(ridge_sig, en_sig), 1, function(x) {
    if(sum(is.na(x)) > 0) {
      return(NA)
    }
    if(x[1] == T & x[2] == T) {
      return('both sig')
    } else if(x[1] == T & x[2] == F) {
      return('ridge sig')
    } else if(x[1] == F & x[2] == T) {
      return('EN sig')
    } else {
      return('not sig')
    }
  })
  p1 = tmp %>% filter(!is.na(EN)) %>% ggplot() + geom_point(aes(x = ridge, y = EN, color = label), alpha = 0.5) + coord_equal() + 
    th2 + geom_abline(slope = 1, intercept = 0, linetype = 2) + facet_wrap(~idp_type) +
    theme(legend.position = 'bottom')
  
  tmp = df %>% reshape2::dcast(IDP + phenotype + idp_type ~ model, value.var = 'pip')
  p2 = tmp %>% ggplot() + geom_point(aes(x = ridge, y = EN), alpha = 0.5) + coord_equal() + 
    th2 + geom_abline(slope = 1, intercept = 0, linetype = 2) + facet_wrap(~idp_type)
  ggsave('sbrainxcan_run/pairwise_zscore.png', p1, width = 5.5, height = 4)
  ggsave('sbrainxcan_run/pairwise_pip.png', p2, width = 5.5, height = 3.5)
  
}

if(isTRUE(gencor_plot)) {
  # load ldsc gen cor
  gtex_gwas = read.table('../../submission/simagexcan/gtex_gwas_list.txt')$V1
  psychiatric = read.table('../preprocess_psychiatric_traits/trait_list.txt')$V1
  df_cor = list()
  tags = list(dmri = 'dmri.original.all_covar.w_pc', t1 = 't1.scaled.all_covar.w_pc')
  source('../../rmd/rlib.R')
  for(rr in psychiatric) {
    for(nn in names(tags)) {
      tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_2nd/', nn, '_2nd_x_psychiatric_x_', rr, '.ldsc_rg.log')
      tmp = load_ldsc_rg(tmp)
      tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
      df_cor[[length(df_cor) + 1]] = tmp %>% 
        mutate(idp_type = nn, phenotype = rr)
    }
  }
  for(rr in gtex_gwas) {
    for(nn in names(tags)) {
      tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_2nd/', nn, '_2nd_x_gtex-gwas_x_', rr, '.ldsc_rg.log')
      tmp = load_ldsc_rg(tmp)
      tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
      df_cor[[length(df_cor) + 1]] = tmp %>% 
        mutate(idp_type = nn, phenotype = rr)
    }
  }
  df_cor = do.call(rbind, df_cor)
  label = rep('T1', nrow(df_cor))
  label[df_cor$idp_type == 'dmri'] = 'dMRI'
  df_cor$idp_type = label
  dfm = df_cor %>% select(IDP, phenotype, idp_type)
  dfm = dfm %>% left_join(
    df_cor %>% select(IDP, phenotype, idp_type, zscore) %>% rename(zscore.rg = zscore),
    by = c('IDP', 'phenotype', 'idp_type')
  )
  dfm = dfm %>% left_join(
    df %>% filter(model == 'EN') %>% select(IDP, phenotype, idp_type, zscore) %>% rename(zscore.EN = zscore),
    by = c('IDP', 'phenotype', 'idp_type')
  )
  dfm = dfm %>% left_join(
    df %>% filter(model == 'ridge') %>% select(IDP, phenotype, idp_type, zscore) %>% rename(zscore.ridge = zscore),
    by = c('IDP', 'phenotype', 'idp_type')
  )
  p1 = dfm %>% ggplot() + geom_point(aes(x = zscore.rg, y = zscore.EN), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'gray') + th
  p2 = dfm %>% ggplot() + geom_point(aes(x = zscore.rg, y = zscore.ridge), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'gray') + th
  ggsave('sbrainxcan_run/vs_gen_cor_en.png', p1, width = 4, height = 4)
  ggsave('sbrainxcan_run/vs_gen_cor_ridge.png', p2, width = 4, height = 4)
}

if(isTRUE(example_plot)) {
  # plot a specific IDP/phenotype pair with MR results
  idp_annot = read.delim2('supp_table_1.tsv', header = T)
  pheno = 'EA_Lee_2018'
  idp = idp_annot %>% filter(
    dmri_measure == 'ICVF', 
    measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)',
    anatomy == 'inferior fronto-occipital fasciculus',
    left_or_right == 'right'
  ) %>% pull(ukb_field)
  idp = paste0('IDP-', idp)
  mr_res = readRDS(paste0('~/Desktop/tmp/ukb_idp/mr_psychiatric_2nd/MR_local.psychiatric_2nd_dmri.', idp, '_x_', pheno, '.rds'))
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
  
  ggsave(paste0('sbrainxcan_run/mr_', idp, '_x_', pheno, '.png'), p1 + p2, width = 8, height = 4)
  
  mr_methods = c('Inverse variance weighted', 'Weighted mode', 'Weighted median')
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
    file =paste0('sbrainxcan_run/mr_', idp, '_x_', pheno, '.tex'), size="\\scriptsize", 
    sanitize.colnames.function=bold,
    booktabs = TRUE
  )
}

