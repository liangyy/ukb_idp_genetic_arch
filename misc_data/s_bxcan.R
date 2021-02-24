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

load_mr = function(dd, df, prefix = '~/Desktop/tmp/ukb_idp/mr_indiv_bxcan/MR.indiv_bxcan_') {
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
      data.frame(model = c('ridge', 'elastic net'), method = c('BrainXcan ridge', 'BrainXcan EN')), 
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

plot_mr = function(model, idp, pheno, prefix = '~/Desktop/tmp/ukb_idp/mr_indiv_bxcan/MR.indiv_bxcan_') {
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

check_sign = function(target, btarget) {
  if(length(unique(sign(btarget))) > 1 | length(unique(sign(target))) > 1) {
    return(F)
  } else {
    bsign = sign(btarget[1])
    sig_ = sign(target[1])
    return(sig_ == bsign)
  }
}

factor_methods = function(x) {
  factor(x, levels = c('ridge', 'elastic net', 'rg'))
}

library(dplyr)
library(patchwork)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib.R')
source('../../rmd/rlib_calc.R')

foldern = 's_bxcan'
dir.create(foldern)

plot_s_bxcan = F 
plot_s_qq = F
prep_s_mr = F
check_s_mr = F
common_pheno_s = F
gen_cor_s = T

color_code = c('ridge' = 'blue', 'elastic net' = 'orange', 'rg' = 'pink')
color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
t1_col = 'orange'  
dmri_col = 'blue'  
factor_idp = function(cc) {
  factor(cc, levels = c('T1', 'dMRI'))
}

df_gwas = read.delim2('supp_table_4.tsv', header = T)
folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd')
df_gwas$folder = 'gtex_gwas'
df_gwas$folder[25:35] = 'psychiatric'
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
    ylab('PIP from elastic net predictors') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
  
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

if(isTRUE(check_s_mr)) {
  idp_meta = read.delim2('supp_table_1.tsv') %>% mutate(IDP = paste0('IDP-', ukb_field))
  mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
  dd_psych = rbind(
    read.table(paste0(foldern, '/s_bxcan_mr.psychiatric.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
    read.table(paste0(foldern, '/s_bxcan_mr.psychiatric.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
  )
  dd_gg = rbind(
    read.table(paste0(foldern, '/s_bxcan_mr.gtex_gwas.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
    read.table(paste0(foldern, '/s_bxcan_mr.gtex_gwas.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
  )
  psych_prefix = '~/Desktop/tmp/ukb_idp/mr_s_bxcan_psych/MR_local.s_bxcan_psych_'
  gg_prefix = '~/Desktop/tmp/ukb_idp/mr_s_bxcan_gtexgwas/MR_local.s_bxcan_gtexgwas_'
  df_mr = rbind(
    load_mr(dd_psych, df, prefix = psych_prefix) %>% mutate(source = 'psychiatric'),
    load_mr(dd_gg, df, prefix = gg_prefix) %>% mutate(source = 'gtex_gwas')
  )
    
  df_mr_entries = df_mr %>% select(phenotype, IDP, idp_type) %>% distinct()
  df_mr_entries = left_join(df_mr_entries, idp_meta %>% select(IDP, notes), by = 'IDP')
  message('distinct pairs passing the criteria: ', nrow(df_mr_entries))
  # idx = 67
  # plot_mr(df_mr_entries$idp_type[idx], df_mr_entries$IDP[idx], df_mr_entries$phenotype[idx], prefix = gg_prefix)
  # df_mr_entries[idx, ] %>% inner_join(df_mr, by = c('idp_type', 'IDP', 'phenotype'))
}

if(isTRUE(common_pheno_s)) {
  common_pairs = list(
    Height = c('UKB_50_Standing_height', 'GIANT_HEIGHT'),
    BMI = c('UKB_21001_Body_mass_index_BMI', 'GIANT_2015_BMI_EUR'),
    Neuroticism = c('UKB_20127_Neuroticism_score', 'Neuroticism_CTG'),
    Depression = c('SSGAC_Depressive_Symptoms', 'MDD_PGC_2018'),
    `Alzheimer's disease` = c('IGAP_Alzheimer', 'AD_Jansen_2019'),
    `Parkinson's disease` = c('UKB_20002_1262_self_reported_parkinsons_disease', 'PD_Nalls_2019')
  )
  df_pairs = list()
  plist = list()
  # plist2 = list()
  for(cc in names(common_pairs)) {
    df_sub = df %>% filter(phenotype %in% common_pairs[[cc]])
    tmp = df_sub %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'zscore')
    # tmp2 = df_sub %>% mutate(pip0 = pmax(1e-3, pip)) %>% reshape2::dcast(IDP + idp_type + model ~ phenotype, value.var = 'pip0')
    ss = df_gwas %>% filter(phenotype_id %in% common_pairs[[cc]])
    ss = ss[match(colnames(tmp)[4:5], ss$phenotype_id), ]
    
    # bigger sample size is on y-axis
    if(ss$sample_size[1] > ss$sample_size[2]) {
      colnames(tmp)[4:5] = c('y', 'x')
      # colnames(tmp2)[4:5] = c('y', 'x')
    } else {
      colnames(tmp)[4:5] = c('x', 'y')
      # colnames(tmp2)[4:5] = c('x', 'y')
    }
    
    
    ss = ss[ order(ss$sample_size, decreasing = T),  ]
    df_pairs[[length(df_pairs) + 1]] = tmp %>% mutate(phenotype = cc)
    plist[[length(plist) + 1]] = tmp %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, color = 'gray') + 
      xlab(ss$short_name[1]) + 
      ylab(ss$short_name[2]) + th +
      ggtitle(cc)
    # plist2[[length(plist2) + 1]] = tmp2 %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, color = 'gray') + 
    #   xlab(ss$short_name[1]) + 
    #   ylab(ss$short_name[2]) + th +
    #   ggtitle(cc) + scale_x_log10() + scale_y_log10()
  }
  df_pairs = do.call(rbind, df_pairs)
  ggsave(
    paste0(foldern, '/', 's_bxca_compare_two_gwas.png'),
    (plist[[1]] + plist[[3]] + plist[[5]]) / (plist[[2]] + plist[[4]] + plist[[6]]),
    width = 8, height = 5.5
  )
  # df_pairs %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, color = 'gray') + facet_wrap(~phenotype, scales = 'free') + th2 
}

if(isTRUE(gen_cor_s)) {
  # load ldsc gen cor
  df_cor = list()
  tags = list(dMRI = 'dmri.original.all_covar.w_pc', T1 = 't1.scaled.all_covar.w_pc')
  mids = list(gtex_gwas = '_2nd_x_gtex-gwas_x_', psychiatric = '_2nd_x_psychiatric_x_')
  for(m in names(mids)) {
    sub = df_gwas %>% filter(folder == m)
    for(i in 1 : nrow(sub)) {
      rr = sub$phenotype_id[i]
      for(nn in names(tags)) {
        tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_2nd/', nn, mids[[m]], rr, '.ldsc_rg.log')
        tmp = load_ldsc_rg(tmp)
        tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
        df_cor[[length(df_cor) + 1]] = tmp %>% 
          mutate(idp_type = nn, phenotype = rr)
      }
    }
  }
  df_cor = do.call(rbind, df_cor)
  df_cor = df_cor %>% mutate(id = paste(phenotype, IDP, idp_type)) %>% 
    filter(id %in% (df %>% mutate(id = paste(phenotype, IDP, idp_type)) %>% pull(id)))
  df_cor_all = df_cor 
  df_cor = left_join(df %>% select(IDP, idp_type, phenotype), df_cor, by = c('IDP', 'phenotype', 'idp_type'))
  
  tmp = rbind(
    df %>% select(IDP, pval, idp_type, phenotype, model),
    df_cor_all %>% select(IDP, pval, idp_type, phenotype) %>% mutate(model = 'rg')
  )
  tmp = tmp %>% group_by(model, idp_type) %>% mutate(pval_cap = pmax(min_pval, pval)) %>% mutate(p_adj = pval * n()) %>% ungroup()
  tmp = tmp %>% group_by(model, idp_type) %>% mutate(pexp = rank(pval) / (n() + 1)) %>% arrange(pval) 
  p = tmp %>% mutate(model = factor_methods(model)) %>% 
    ggplot() + 
    geom_path(aes(x = -log10(pexp), y = -log10(pval_cap), color = model)) + 
    geom_point(data = tmp %>% filter(p_adj < 0.05), aes(x = -log10(pexp), y = -log10(pval_cap), color = model)) +
    facet_wrap(~idp_type) + 
    geom_abline(slope = 1, intercept = 0) + th2 +
    scale_color_manual(values = color_code) + 
    xlab(expression(paste(-log[10], p[expected]))) + 
    ylab(expression(paste(-log[10], p[observed]))) +
    theme(legend.position = c(0.89, 0.68), legend.title = element_blank())
  ggsave(paste0(foldern, '/s_bxcan_qqplot_w_gencor.png'), p, width = 6, height = 3)
  
}

