# setwd('misc_data/supplementary_materials_3rd/')

load_mr = function(dd, df, prefix = '~/Desktop/tmp/ukb_idp/mr_scz2020_3rd/MR_local.scz2020_') {
  df_mr = list()
  for(i in 1 : nrow(dd)) {
    tmp = readRDS(paste0(prefix, dd$idp_type[i], '.', dd$idp[i], '_x_', dd$pheno[i], '.rds'))
    df_mr0 = rbind(
      tmp$idp2pheno$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype'),
      tmp$pheno2idp$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP') 
    ) %>% select(direction, method, nsnp, b, pval)
    kk = df_mr0 %>% group_by(direction) %>% summarize(nsig = sum(pval < 0.05), sign = max(sum(b > 0), sum(b <= 0))) %>% ungroup()
    if(max(kk$nsig) < 2 | sum(kk$sign[kk$nsig >= 2] == 3) == 0) {
      next
    }
    tmp2 = inner_join(
      data.frame(model = c('ridge', 'elastic net', 'genetic correlation'), method = c('BrainXcan ridge', 'BrainXcan EN', 'Genetic Correlation')), 
      df %>% filter(phenotype == dd$pheno[i], IDP == dd$idp[i], tolower(idp_type) == dd$idp_type[i]),
      by = 'model'
    )
    df_mr0 = rbind(
      df_mr0 %>% mutate(pip = NA), data.frame(direction = NA, method = tmp2$method, nsnp = NA, b = tmp2$bhat, pval = tmp2$pval, pip = tmp2$pip)
    )
    tmp = df_mr0$b[7:nrow(df_mr0)]
    tmp = tmp[df_mr0$pval[7:nrow(df_mr0)] < 5e-2]
    if(
      (!check_sign(df_mr0$b[1:3], tmp)) &
      (!check_sign(df_mr0$b[4:6], tmp))
    ) {
      next
    }
    df_mr[[length(df_mr) + 1]] = df_mr0 %>% mutate(phenotype = dd$pheno[i], IDP = dd$idp[i], idp_type = dd$idp_type[i])
  }
  df_mr = do.call(rbind, df_mr)
  df_mr
}

plot_mr = function(model, idp, pheno, prefix = '~/Desktop/tmp/ukb_idp/mr_scz2020_3rd/MR_local.scz2020_') {
  mr_res = readRDS(paste0(prefix, model, '.', idp, '_x_', pheno, '.rds'))
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
  list(p1 = p1, p2 = p2)
}

save_mr_table = function(df_mr, filename) {
  print(
    xtable::xtable(
      df_mr %>% select(IDP, model, direction, method, nsnp, b, pval, pip), 
      display=c("s", "s", "s", "s", "s", "g", "g", 'g', 'g'), digits = 3
    ), 
    math.style.exponents = TRUE, 
    include.rownames = FALSE, 
    file = filename, 
    sanitize.colnames.function=bold,
    booktabs = TRUE
  )
}

idp_type = list(dMRI = 'third_round_dmri', T1 = 'third_round_t1')
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
source('rlib.R')

plot_overview = T
mr_prep = T
mr_check = F
compare_scz2 = T
not_run = T  # playground, skip if you'd like

pcs = get_pcs_for_3rd()

foldern = 'scz_vignette'
dir.create(foldern)

psychiatric = c('SCZ_PGC_2020', 'pgc.scz2')  


df1 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_3rd/', psychiatric[2])
df2 = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_3rd/', psychiatric[1])
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
  tags = list(dmri = 'third_round_dmri', t1 = 'third_round_t1')
  for(rr in psychiatric[1]) {
    for(nn in names(tags)) {
      tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_3rd/', nn, '_3rd_x_psychiatric_x_', rr, '.ldsc_rg.log')
      tmp = load_ldsc_rg(tmp)
      tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
      df_cor[[length(df_cor) + 1]] = tmp %>% 
        mutate(idp_type = nn, phenotype = rr)
    }
  }
  for(rr in psychiatric[2]) {
    for(nn in names(tags)) {
      tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_3rd/', nn, '_3rd_x_gtex-gwas_x_', rr, '.ldsc_rg.log')
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
    df %>% select(IDP, phenotype, idp_type, model, zscore, pval, pip, bhat), 
    df_cor %>% mutate(model = 'genetic correlation', pip = NA) %>% select(IDP, phenotype, idp_type, model, zscore, pval, pip, rg) %>% rename(bhat = rg)
  )
  
  df = df %>% group_by(phenotype, model) %>% mutate(p_adj = pval * n()) %>% ungroup()
  
}

# load annot
{
  annot = read.delim2('supp_table_1.tsv', header = T) %>% mutate(IDP = paste0('IDP-', ukb_field))
  annot = rbind(
    annot %>% select(IDP, t1_or_dmri, anatomy, left_or_right, measurement_type, dmri_measure, t1_anatomy_group),
    data.frame(
      IDP = pcs$IDP,  
      t1_or_dmri = pcs$idp_type,
      anatomy = pcs$IDP, 
      left_or_right = NA,
      measurement_type = paste0(pcs$idp_type, ' PCA'),
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

if(isTRUE(plot_overview)) {
  p = df %>% filter(phenotype == psychiatric[1]) %>% ggplot() + 
    geom_point(
      data = df %>% filter(phenotype == psychiatric[1]) %>% 
        filter((pip > 0.5 & p_adj < alpha) | (model == 'genetic correlation' & p_adj < alpha)),
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
if(isTRUE(mr_prep)) {
  df_sig = df %>% filter(phenotype == psychiatric[1]) %>% filter((pip > 0.5 & p_adj < alpha) | (model == 'genetic correlation' & p_adj < alpha))
  # kk = df_sig %>% mutate(p2 = phenotype) %>% select(phenotype, IDP, p2)
  for(i in c('T1', 'dMRI')) {
    tmp = df_sig %>% filter(idp_type == i) %>% mutate(pp = phenotype) %>% select(phenotype, IDP, pp) %>% rename(pheno = phenotype, idp = IDP, pheno_code = pp) %>% distinct()
    write.table(tmp, paste0(foldern, '/scz2020.', i, '.signif.tsv'), quo = F, col = T, row = F, sep = '\t')
  }
}

if(isTRUE(mr_check)) {
  idp_meta = read.delim2('supp_table_1.tsv') %>% mutate(IDP = paste0('IDP-', ukb_field))
  mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
  dd = rbind(
    read.table(paste0(foldern, '/scz2020.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
    read.table(paste0(foldern, '/scz2020.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
  )
  
  df_mr = load_mr(dd, df)
  
  df_mr$idp_type[df_mr$idp_type == 'dmri'] = 'dMRI'
  df_mr$idp_type[df_mr$idp_type == 't1'] = 'T1'
  
  df_mr = df_mr %>% left_join(
    df_sig %>% select(IDP, phenotype, idp_type, model, anatomy, left_or_right, measurement_type, dmri_measure, t1_anatomy_group), 
    by = c('IDP', 'phenotype', 'idp_type'))
  
  df_mr_entries = df_mr %>% select(phenotype, IDP, idp_type) %>% distinct()
  df_mr_entries = left_join(df_mr_entries, idp_meta %>% select(IDP, notes), by = 'IDP')
  message('distinct pairs passing the criteria: ', nrow(df_mr_entries))
  message('distinct GWASs passing the criteria: ', length(unique(df_mr_entries$phenotype)))
  message('distinct GWASs (psychiatric) passing the criteria: ', sum(!unique(df_mr_entries$phenotype) %in% not_psych))
  df_sig %>% left_join(df_mr_entries %>% select(-notes) %>% mutate(pass_mr = T), by = c('IDP', 'phenotype', 'idp_type')) %>% group_by(model) %>% summarize(sum(pass_mr, na.rm = T))
  # idx = 1
  tt = list()
  for(idx in 1 : nrow(df_mr_entries)) {
    pp = plot_mr(df_mr_entries$idp_type[idx], df_mr_entries$IDP[idx], df_mr_entries$phenotype[idx])
    tmp = df_mr_entries[idx, ] %>% inner_join(df_mr, by = c('idp_type', 'IDP', 'phenotype'))
    tt[[length(tt) + 1]] = tmp
    # save_mr_table(tmp, paste0(foldern, '/scz2020_x_', df_mr_entries$IDP[idx], '_', df_mr_entries$idp_type[idx], '.tex'))
    ggsave(paste0(foldern, '/scz2020_x_', df_mr_entries$IDP[idx], '_', df_mr_entries$idp_type[idx], '.png'), pp[[1]] + pp[[2]], height = 4, width = 8)
  }
  tt = do.call(rbind, tt)
  write.table(tt %>% select(IDP, model, direction, method, nsnp, b, pval, pip, notes), 'supp_table_5.tsv', quo = F, col = T, row = F, sep = '\t')
  # save_mr_table(tt, paste0(foldern, '/scz2020_mr.tex'))
}

if(isTRUE(compare_scz2)) {
  # p = df %>% filter(phenotype == psychiatric[2]) %>% ggplot() + 
  #   geom_point(
  #     data = df %>% filter(phenotype == psychiatric[2], pip > 0.5),
  #     aes(x = idp_f, y = -log10(pval)), size = 4, shape = 1
  #   ) +
  #   geom_point(aes(x = idp_f, y = -log10(pval), color = measurement_type)) +
  #   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  #   theme(legend.position = 'right', legend.title = element_blank()) +
  #   th2 +
  #   scale_color_manual(values = cbPalette) +
  #   facet_wrap(~model, ncol = 1) +
  #   geom_hline(
  #     data = df %>% filter(phenotype == psychiatric[1]) %>% 
  #       group_by(model) %>% summarize(cutoff = alpha / n()) %>% ungroup(),
  #     aes(yintercept = -log10(cutoff)), linetype = 2
  #   ) + 
  #   xlab('Brain IDPs') + 
  #   ylab(expression(paste(-log[10], p[S-BrainXcan]))) 
  # ggsave(paste0(foldern, '/scz_overview_scz2.png'), p, width = 8, height = 5)
  p1 = df %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'zscore') %>%
    ggplot() + geom_point(aes(x = pgc.scz2, y = SCZ_PGC_2020), alpha = 0.5) + facet_wrap(~model) + 
    geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 +
    ylab('S-BrainXcan PIP from \n Ripke et al. (2020)') +
    xlab('S-BrainXcan PIP from Ripke et al. (2014)')
  p2 = df %>% mutate(pip = pmax(pip, 1e-3)) %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'pip') %>%
    filter(model != 'genetic correlation') %>% 
    ggplot() + geom_point(aes(x = pgc.scz2, y = SCZ_PGC_2020), alpha = 0.5) + 
    facet_wrap(~model) + 
    geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 + scale_x_log10() + scale_y_log10() +
    ylab('S-BrainXcan PIP from \n Ripke et al. (2020)') +
    xlab('S-BrainXcan PIP from Ripke et al. (2014)')
  ggsave(paste0(foldern, '/scz2020_vs_scz2_pip.png'), p2, height = 3, width = 6)
  ggsave(paste0(foldern, '/scz2020_vs_scz2_zscore.png'), p1, height = 3, width = 6)
}

# -------------- playgroud -------------------

# try visualization
no_run = T
if(!isTRUE(no_run)) {
  df %>% filter(idp_type == 'T1', phenotype == psychiatric[1]) %>% 
    ggplot() + geom_point(aes(x = idp_f, y = -log10(pval), color = model)) + facet_grid(~measurement_type, scales = 'free_x', space = 'free_x') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + th2
  p = df %>% filter(idp_type == 'dMRI', phenotype == psychiatric[1]) %>% 
    ggplot() + geom_point(aes(x = idp_f, y = -log10(pval), color = model)) + facet_wrap(~measurement_type, scales = 'free_x', ncol = 1) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + th2 
}

# en vs ridge
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
