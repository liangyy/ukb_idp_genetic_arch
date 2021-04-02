# setwd('misc_data/supplementary_materials_3rd/')

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

library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib.R')
library(UpSetR)

foldern = 'indiv_bxcan'
dir.create(foldern)

color_code = c('ridge' = 'blue', 'elastic net' = 'orange')
color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
t1_col = 'orange'  # rgb()
dmri_col = 'blue'  # rgb()
factor_idp = function(cc) {
  factor(cc, levels = c('T1', 'dMRI'))
}


do_sbxcan_compare = T
plot_i_bxcan = T
do_mr_prep = T
plot_qq = F
check_mr_result = F#T
save_df = F#T
save_df_full = F#T

pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes', 'bmi', 'height')
models = list(ridge = 'ridge', EN = 'en')
idps = list(dMRI = 'dmri', T1 = 't1') #  
# types = c('linear', 'susie')
df = list()
for(m in names(models)) {
  for(t in names(idps)) {
    fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_round_3.linear.', idps[[t]], '_', models[[m]], '.csv')
    tmp1 = read.csv(fn1)
    fn2 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_round_3.susie.', idps[[t]], '_', models[[m]], '.csv')
    tmp2 = read.csv(fn2)
    tmp = inner_join(tmp1, tmp2, by = c('IDP', 'phenotype'))
    df[[length(df) + 1]] = tmp %>% mutate(model = m, idp_type = t)
  }
}
df = do.call(rbind, df)
df$bhat = - df$bhat
df_all = df
df_all$model[df_all$model == 'EN'] = 'elastic net'
df = df_all %>% filter(phenotype %in% pheno_interest) %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))
if(isTRUE(save_df_full)) {
  saveRDS(df %>% select(-idp_id), paste0(foldern, '/dataframe_full.indiv_bxcan.rds'))
}


idp_sig = read.table('supp_table_2.tsv', header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
df = df[df$idp_id %in% idp_sig$idp_id, ]
alpha = 0.05
min_pval = 1e-30
df = df %>% group_by(model, idp_type) %>% mutate(pval_cap = pmax(min_pval, pval)) %>% mutate(p_adj = pval * n()) %>% ungroup()

if(isTRUE(save_df)) {
  saveRDS(df %>% select(-pval_cap, -p_adj, -idp_id), paste0(foldern, '/dataframe.indiv_bxcan.rds'))
}

supp3 = 'supp_table_3.tsv'
if(!file.exists(supp3)) {
  # build the phenotype meta table
  pheno_table = list()
  
  # weekly_alcohol
  idx = 1
  ukb_field = c(1568, 1578, 1588, 1598, 1608)
  target_value = 'asis'
  aggregate_method = c('sum')
  missing_values = c(0)
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # recurrent_depressive_disorder
  idx = 2
  ukb_field = c(20002, 20124, 41202, 41204)
  target_value = c(1286, 1, 'F330 ~ F334 and F338 F339', 'F330 ~ F334 and F338 F339')
  aggregate_method = c('or')
  missing_values = c(0)
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # parent_depression
  idx = 3
  ukb_field = c(20107, 20110)
  target_value = c(12, 12)
  aggregate_method = c('or')
  missing_values = c('remove')
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # parent_AD
  idx = 4
  ukb_field = c(20107, 20110)
  target_value = c(10, 10)
  aggregate_method = c('or')
  missing_values = c('remove')
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # handedness
  idx = 5
  ukb_field = c(1707)
  target_value = c(2)
  aggregate_method = c('none')
  missing_values = c('remove')
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # daily_coffee
  idx = 6
  ukb_field = c(1498)
  target_value = c('asis')
  aggregate_method = c('none')
  missing_values = c(0)
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  # daily_cigarettes
  idx = 7
  ukb_field = c(3456)
  target_value = c('asis')
  aggregate_method = c('none')
  missing_values = c(0)
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)  
  )
  
  
  # bmi
  idx = 8
  ukb_field = c(21001)
  target_value = c('asis')
  aggregate_method = c('none')
  missing_values = c('remove')
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values) 
  )
  
  
  # height
  idx = 9
  ukb_field = c(50)
  target_value = c('asis')
  aggregate_method = c('none')
  missing_values = c('remove')
  pheno_table[[length(pheno_table) + 1]] = data.frame(
    phenotype = pheno_interest[idx], 
    ukb_field = ukb_field,
    target_value = as.character(target_value),
    aggregate_method = aggregate_method,
    missing_values = as.character(missing_values)
  )
  
  
  pheno_table = do.call(rbind, pheno_table)
  
  write.table(pheno_table, supp3, sep = '\t', quote = F, col = T, row = F)
}

if(isTRUE(plot_i_bxcan)) {
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
  
  ggsave(paste0(foldern, '/indiv_bxcan_ridge_vs_en_zscore.png'), p1, width = 5.5, height = 3)
  ggsave(paste0(foldern, '/indiv_bxcan_ridge_vs_en_pip.png'), p2, width = 5.5, height = 3)
}

if(isTRUE(plot_qq)) {
  tmp = df %>% group_by(model, idp_type) %>% mutate(pexp = rank(pval) / (n() + 1)) %>% arrange(pval) 
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
  ggsave(paste0(foldern, '/indiv_bxcan_qqplot.png'), p, width = 5.5, height = 3)
}

# compare to S-BrainXcan
if(isTRUE(do_sbxcan_compare)) {
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
  phenos = c('UKB_50_Standing_height', 'UKB_21001_Body_mass_index_BMI')
  df_sb = load_sbxcan('~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_3rd/', phenos)
  df_sb = df_sb %>% mutate(zscore = p2z(pval, bhat))
  df_sb$model[ df_sb$model == 'EN' ] = 'elastic net'
  df_sb = inner_join(df_sb, data.frame(pheno = c('height', 'bmi'), phenotype = c('UKB_50_Standing_height', 'UKB_21001_Body_mass_index_BMI')), by = c('phenotype'))
  df_in = df %>% filter(phenotype %in% c('height', 'bmi'))
  df_in = df_in %>% mutate(zscore = p2z(pval, bhat))
  df_both = inner_join(
    df_sb %>% select(IDP, bhat, zscore, pip, cs95, model, idp_type, pheno), 
    df_in %>% select(IDP, bhat, zscore, pip, cs, model, idp_type, phenotype) %>% rename(cs95 = cs), 
    by = c('IDP', 'idp_type', 'model', 'pheno' = 'phenotype'),
    suffix = c('.sb', '.in'))
  
  p = df_both %>% ggplot() + geom_point(aes(x = zscore.in, y = zscore.sb, color = idp_type), alpha = 0.5) + th2 + geom_abline(slope = 1, intercept = 0, color = 'gray') + coord_equal() + facet_wrap(~model) + 
    xlab('BrainXcan z-score') + 
    ylab('S-BrainXcan z-score') +
    scale_color_manual(values = color_code2) +
    theme(legend.position = c(0.65, 0.8), legend.title = element_blank())
  ggsave(paste0(foldern, '/ukb_height_bmi_comparison.png'), p, width = 5, height = 3)
  
  tmp2 = df_both %>% mutate(indiv_BrainXcan = pip.in > 0.5, S_BrainXcan = pip.sb > 0.5) %>%
    filter(indiv_BrainXcan | S_BrainXcan) 
  tmp2_ = tmp2 %>% select(IDP, idp_type, model, indiv_BrainXcan, S_BrainXcan) 
  kk = tmp2__ = as.matrix(tmp2_[, 4:5])
  tmp2__[kk] = 1
  tmp2__[!kk] = 0
  tmp2_[, 4:5] = tmp2__
  tmp2_ = as.data.frame(tmp2_)
  
  png(paste0(foldern, '/', 'ukb_height_bmi_pip_en.png'), width = 7, height = 5, units = 'in', res = 300)
  upset(tmp2_ %>% filter(model == 'elastic net'), main.bar.color = color_code['elastic net'], point.size = 5, line.size = 0.5, text.scale = 2.75)
  dev.off()
  
  png(paste0(foldern, '/', 'ukb_height_bmi_pip_ridge.png'), width = 7, height = 5, units = 'in', res = 300)
  upset(tmp2_ %>% filter(model == 'ridge'), main.bar.color = color_code['ridge'], point.size = 5, line.size = 0.5, text.scale = 2.75)
  dev.off()
  
}

if(isTRUE(do_mr_prep)) {
  df_sig = df %>% filter(p_adj < alpha, pip > 0.5)
  df_gwas = readRDS('../selected_pheno_to_open_gwas.rds')
  df_sig = left_join(df_sig, df_gwas %>% select(phenotype, gwas_code), by = 'phenotype')
  for(i in c('T1', 'dMRI')) {
    tmp = df_sig %>% filter(idp_type == i) %>% select(phenotype, IDP, gwas_code) %>% rename(pheno = phenotype, idp = IDP, pheno_code = gwas_code) %>% distinct()
    write.table(tmp, paste0(foldern, '/indiv_bxcan_mr.', i, '.signif.tsv'), quo = F, col = T, row = F, sep = '\t')
  }
  
}

if(isTRUE(check_mr_result)) {
  idp_meta = read.delim2('supp_table_1.tsv') %>% mutate(IDP = paste0('IDP-', ukb_field))
  mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
  dd = rbind(
    read.table(paste0(foldern, '/indiv_bxcan_mr.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
    read.table(paste0(foldern, '/indiv_bxcan_mr.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
  )
  df_mr = load_mr(dd, df)
  df_mr_entries = df_mr %>% select(phenotype, IDP, idp_type) %>% distinct()
  df_mr_entries = left_join(df_mr_entries, idp_meta %>% select(IDP, notes), by = 'IDP')
  message('distinct pairs passing the criteria: ', nrow(df_mr_entries))
  # idx = 15
  # plot_mr(df_mr_entries$idp_type[idx], df_mr_entries$IDP[idx], df_mr_entries$phenotype[idx])
  # df_mr_entries[idx, ] %>% inner_join(df_mr, by = c('idp_type', 'IDP', 'phenotype'))
}
