# setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/')

load_sbxcan = function(folder, trait_list) {
  idp_type = list(dMRI = 'dmri', T1 = 't1')
  models = list(ridge = 'ridge', EN = 'en')
  df1 = list()
  for(t in trait_list) {
    for(i in names(idp_type)) {
      for(m in names(models)) {
        df1[[length(df1) + 1]] = read.csv(paste0(folder, '_', models[[m]], '_residual', '/', idp_type[[i]], '_x_', t, '_x_simagexcan.csv'), header = T) %>% mutate(idp_type = i, model = m, phenotype = t)
      }
    }
  }
  df1 = do.call(rbind, df1)
  df1
}

load_sbxcan_permz <- function(folder, trait_list) {
  models = list(ridge = 'ridge', EN = 'en')
  df1 = list()
  for(t in trait_list) {
    for(m in names(models)) {
      df1[[length(df1) + 1]] = read.csv(paste0(folder, '_', models[[m]], '_permz', '/', t, '.sbrainxcan.csv'), header = T) %>% mutate(model = m, phenotype = t)
    }
  }
  df1 = do.call(rbind, df1)
}
# 
# load_mr = function(dd, df, prefix = '~/Desktop/tmp/ukb_idp/mr_indiv_bxcan/MR.indiv_bxcan_') {
#   df_mr = list()
#   for(i in 1 : nrow(dd)) {
#     tmp = readRDS(paste0(prefix, dd$idp_type[i], '_3rd.', dd$idp[i], '_x_', dd$pheno[i], '.rds'))
#     df_mr0 = rbind(
#       tmp$idp2pheno$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'IDP -> Phenotype'),
#       tmp$pheno2idp$mr %>% filter(method %in% mr_methods) %>% mutate(direction = 'Phenotype -> IDP') 
#     ) %>% select(direction, method, nsnp, b, pval)
#     kk = df_mr0 %>% group_by(direction) %>% summarize(nsig = sum(pval < 0.05), sign = max(sum(b > 0), sum(b <= 0))) %>% ungroup()
#     if(max(kk$nsig) < 2 | sum(kk$sign[kk$nsig >= 2] == 3) == 0) {
#       next
#     }
#     tmp2 = inner_join(
#       data.frame(model = c('ridge', 'elastic net'), method = c('BrainXcan ridge', 'BrainXcan EN')), 
#       df %>% filter(phenotype == dd$pheno[i], IDP == dd$idp[i], tolower(idp_type) == dd$idp_type[i]),
#       by = 'model'
#     )
#     df_mr0 = rbind(
#       df_mr0 %>% mutate(pip = NA), data.frame(direction = NA, method = tmp2$method, nsnp = NA, b = tmp2$bhat, pval = tmp2$pval, pip = tmp2$pip)
#     )
#     tmp = df_mr0$b[7:nrow(df_mr0)]
#     tmp = tmp[df_mr0$pval[7:nrow(df_mr0)] < 5e-2]
#     if(
#       (!check_sign(df_mr0$b[1:3], tmp)) &
#       (!check_sign(df_mr0$b[4:6], tmp))
#     ) {
#       next
#     }
#     df_mr[[length(df_mr) + 1]] = df_mr0 %>% mutate(phenotype = dd$pheno[i], IDP = dd$idp[i], idp_type = dd$idp_type[i])
#   }
#   df_mr = do.call(rbind, df_mr)
#   df_mr
# }
# 
# plot_mr = function(model, idp, pheno, prefix = '~/Desktop/tmp/ukb_idp/mr_indiv_bxcan/MR.indiv_bxcan_') {
#   mr_res = readRDS(paste0(prefix, model, '.', idp, '_x_', pheno, '.rds'))
#   p1 = mr_res$idp2pheno$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
#     geom_vline(xintercept = 0, color = 'gray') + 
#     geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
#     geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
#     geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
#     th + xlab('SNP estimated effect in IDP') + ylab('SNP estimated effect in phenotype') + ggtitle('MR: IDP -> Phenotype')
#   p2 = mr_res$pheno2idp$data %>% ggplot() + geom_hline(yintercept = 0, color = 'gray') + 
#     geom_vline(xintercept = 0, color = 'gray') + 
#     geom_point(aes(x = beta.exposure, y = beta.outcome), alpha = 0.5) + 
#     geom_errorbar(aes(x = beta.exposure, ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.5) +
#     geom_errorbarh(aes(y = beta.outcome, xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.5) +
#     th + xlab('SNP estimated effect in phenotype') + ylab('SNP estimated effect in IDP') + ggtitle('MR: Phenotype -> IDP')
#   list(p1 = p1, p2 = p2)
# }
# 
# check_sign = function(target, btarget) {
#   if(length(unique(sign(btarget))) > 1 | length(unique(sign(target))) > 1) {
#     return(F)
#   } else {
#     bsign = sign(btarget[1])
#     sig_ = sign(target[1])
#     return(sig_ == bsign)
#   }
# }

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
source('rlib.R')

foldern = 's_bxcan_twasinf'
dir.create(foldern)
mydir <- '~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th'

common_pheno_s = T#T
gen_cor_s = T#T
gen_cor_qq_s = T#T
save_df_full = T#T
gen_cor_compare = T#T

plot_s_bxcan = F
plot_s_qq = F
prep_s_mr = F
check_s_mr = F
gen_cor_mr_s = F
gen_cor_check_mr = F
plot_sig_mr = F
save_df = F
correction_factor_perm = 1.1

not_psych = c('GIANT_2015_BMI_EUR', 'GIANT_2017_BMI_Active_EUR', 'GIANT_HEIGHT', 'UKB_50_Standing_height', 'UKB_21001_Body_mass_index_BMI')
color_code = c('ridge' = 'blue', 'elastic net' = 'orange', 'rg' = 'pink')
color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
t1_col = 'orange'  
dmri_col = 'blue'  
factor_idp = function(cc) {
  factor(cc, levels = c('T1', 'dMRI'))
}
load_h2 <- function(pheno, folder) {
  mydir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/genetic_cor_4th'
  fn <- Sys.glob(paste0(mydir, '/dmri_4th_x_*_x_', pheno, '.msg.log'))
  cmd <- system(glue::glue('cat {fn} | grep "Heritability of phenotype 1$" -A 6 '), intern = TRUE)
  h2 <- stringr::str_match(cmd[3], 'Total Observed scale h2: ([0-9.]+) \\(([0-9.]+)\\)')
  as.numeric(h2[, 2])
}
df_gwas = read.delim2(paste0(mydir, '/../supplementary_materials_3rd/supp_table_4.tsv'), header = T)
folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_4th', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_4th')
df_gwas$folder = 'gtex_gwas'
df_gwas$folder[25:35] = 'psychiatric'
df = list()
for(cc in c('gtex_gwas', 'psychiatric')) {
  traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
  df[[length(df) + 1]] = load_sbxcan_permz(folders[[cc]], traits) %>% mutate(source = cc)
}
df = do.call(rbind, df)
df$model[df$model == 'EN'] = 'elastic net'
df_gwas$h2 <- sapply(df_gwas$phenotype_id, load_h2)
df_gwas$h2[df_gwas$phenotype_id == 'SCZ_PGC_2020'] <- 0.7261
df_phi <- read.table('~/Downloads/psychiatric_permz_rerun 6/idps-phi.txt', header = TRUE)
# ss <- import('scipy.stats')
# ss$chi2$logsf(2500, 1)
# pchisq(2500, 1, log.p = T, lower.tail = F)
df <- df %>% left_join(df_gwas %>% select(phenotype_id, sample_size, h2), by = c('phenotype' = 'phenotype_id')) %>%
  left_join(df_phi %>% select(IDP, phi), by = 'IDP') %>%
  mutate(phi = ifelse(phi < 0, 0, phi)) %>% 
  mutate(corrected_pvalue = pchisq(z_brainxcan ^ 2, df= 1, 
                                   ncp = (phi * sample_size * h2),
                                   lower.tail = F, log.p = TRUE),
         corrected_zscore = sqrt(qchisq(corrected_pvalue, df = 1,
                                        lower.tail = F)) * sign(z_brainxcan)) %>% 
  dplyr::rename(pval_corrected = corrected_pvalue, 
                z_corrected = corrected_zscore)
# df = df %>% mutate(idp_id = paste(IDP, model)) %>% mutate(zscore = p2z(pval, bhat)) %>%
#   mutate(z_adj_perm_null = z_adj_perm_null / correction_factor_perm) %>%
#   mutate(pval_adj_perm_null = z2p(z_adj_perm_null))
if(isTRUE(save_df_full)) {
  saveRDS(df, paste0(foldern, '/dataframe_full.sbxcan_twasinf.rds'))
  # saveRDS(df_cor, paste0(foldern, '/dataframe_full.gencor.rds'))
}
df0 <- readRDS(paste0(mydir, '/s_bxcan/dataframe_full.sbxcan.rds'))
# 
df <- left_join(df, df0 %>% select(IDP, model, phenotype, pip, cs95), by = c('IDP', 'model', 'phenotype'))
# summary(tmp$zscore.x - tmp$zscore.y)

# remove ProbTrack IDPs
df = remove_probtrack_idp(df)
df <- df %>% 
  rename(pval_raw = pval, zscore_raw = z_brainxcan, idp_type = modality) %>%
  rename(pval = pval_corrected, zscore = z_corrected) %>% 
  mutate(idp_id = paste(IDP, model))

idp_sig = read.table(paste0(mydir, '/supp_table_2.tsv'), header = T, sep = '\t')
idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP, model_name))
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
  dd_sbxcan = rbind(dd_psych, dd_gg)
  psych_prefix = '~/Desktop/tmp/ukb_idp/mr_s_bxcan_psych_3rd/MR_local.s_bxcan_psych_'
  gg_prefix = '~/Desktop/tmp/ukb_idp/mr_s_bxcan_gtexgwas_3rd/MR_local.s_bxcan_gtexgwas_'
  df_mr = rbind(
    load_mr(dd_psych, df, prefix = psych_prefix) %>% mutate(source = 'psychiatric'),
    load_mr(dd_gg, df, prefix = gg_prefix) %>% mutate(source = 'gtex_gwas')
  )
  
  df_mr_entries = df_mr %>% select(phenotype, IDP, idp_type) %>% distinct()
  df_mr_entries = left_join(df_mr_entries, idp_meta %>% select(IDP, notes), by = 'IDP')
  message('distinct pairs passing the criteria: ', nrow(df_mr_entries))
  message('distinct GWASs passing the criteria: ', length(unique(df_mr_entries$phenotype)))
  message('distinct GWASs (psychiatric) passing the criteria: ', sum(!unique(df_mr_entries$phenotype) %in% not_psych))
  # idx = 67
  # plot_mr(df_mr_entries$idp_type[idx], df_mr_entries$IDP[idx], df_mr_entries$phenotype[idx], prefix = gg_prefix)
  # df_mr_entries[idx, ] %>% inner_join(df_mr, by = c('idp_type', 'IDP', 'phenotype'))
}

if(isTRUE(common_pheno_s)) {
  common_pairs = list(
    Height = c('UKB_50_Standing_height', 'GIANT_HEIGHT'),
    BMI = c('UKB_21001_Body_mass_index_BMI', 'GIANT_2015_BMI_EUR'),
    Neuroticism = c('UKB_20127_Neuroticism_score', 'Neuroticism_CTG'),
    Intelligence = c('Intelligence_CTG', 'UKB_20016_Fluid_intelligence_score'),
    Depression = c('SSGAC_Depressive_Symptoms', 'MDD_PGC_2018'),
    `Alzheimer's disease` = c('IGAP_Alzheimer', 'AD_Jansen_2019')
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
    print(ss)
    # bigger sample size is on y-axis
    if(ss$sample_size[1] > ss$sample_size[2]) {
      namex = ss$short_name[ss$phenotype_id == colnames(tmp)[5]]
      namey = ss$short_name[ss$phenotype_id == colnames(tmp)[4]]
      colnames(tmp)[4:5] = c('y', 'x')
      # colnames(tmp2)[4:5] = c('y', 'x')
    } else {
      namex = ss$short_name[ss$phenotype_id == colnames(tmp)[4]]
      namey = ss$short_name[ss$phenotype_id == colnames(tmp)[5]]
      colnames(tmp)[4:5] = c('x', 'y')
      # colnames(tmp2)[4:5] = c('x', 'y')
    }
    
    
    ss = ss[ order(ss$sample_size, decreasing = T),  ]
    df_pairs[[length(df_pairs) + 1]] = tmp %>% mutate(phenotype = cc)
    plist[[length(plist) + 1]] = tmp %>% filter(model == 'ridge') %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + 
      # geom_abline(slope = 1, intercept = 0, color = 'gray') + 
      xlab(namex) + 
      ylab(namey) + th +
      ggtitle(cc)
    # plist2[[length(plist2) + 1]] = tmp2 %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, color = 'gray') + 
    #   xlab(ss$short_name[1]) + 
    #   ylab(ss$short_name[2]) + th +
    #   ggtitle(cc) + scale_x_log10() + scale_y_log10()
  }
  df_pairs = do.call(rbind, df_pairs)
  ggsave(
    paste0(foldern, '/', 's_bxcan_compare_two_gwas_ridge_only.png'),
    (plist[[1]] + plist[[3]] + plist[[5]]) / (plist[[2]] + plist[[4]] + plist[[6]]),
    width = 8, height = 5.5
  )
  # df_pairs %>% ggplot() + geom_point(aes(x = x, y = y), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, color = 'gray') + facet_wrap(~phenotype, scales = 'free') + th2 
}

if(isTRUE(gen_cor_s)) {
  # load ldsc gen cor
  df_cor = list()
  tags = list(dMRI = 'dmri', T1 = 't1')
  mids = list(gtex_gwas = '_4th_x_gtex-gwas_x_', psychiatric = '_4th_x_psychiatric_x_')
  for(m in names(mids)) {
    sub = df_gwas %>% filter(folder == m)
    for(i in 1 : nrow(sub)) {
      rr = sub$phenotype_id[i]
      for(nn in names(tags)) {
        tmp = paste0('~/Desktop/tmp/ukb_idp/genetic_cor_4th/', nn, mids[[m]], rr, '.ldsc_rg.log')
        tmp = load_ldsc_rg(tmp)
        tmp = tmp %>% select(p2, rg, p, z, h2_obs, h2_obs_se) %>% rename(IDP = p2, pval = p, zscore = z)
        df_cor[[length(df_cor) + 1]] = tmp %>% 
          mutate(idp_type = nn, phenotype = rr)
      }
    }
  }
  df_cor = do.call(rbind, df_cor)
  if(isTRUE(save_df_full)) {
    # saveRDS(df %>% select(-pval_cap, -p_adj, -idp_id), paste0(foldern, '/dataframe.sbxcan.rds'))
    saveRDS(df_cor, paste0(foldern, '/dataframe_full.gencor.rds'))
  }
  df_cor = remove_probtrack_idp(df_cor)
  df_cor = df_cor %>% mutate(id = paste(phenotype, IDP, idp_type)) %>% 
    filter(id %in% (df %>% mutate(id = paste(phenotype, IDP, idp_type)) %>% pull(id)))
  df_cor_all = df_cor %>% left_join(df_gwas %>% select(phenotype_id, folder), by = c('phenotype' = 'phenotype_id'))
  # df_cor = left_join(df %>% select(IDP, idp_type, phenotype), df_cor, by = c('IDP', 'phenotype', 'idp_type'))
  
  if(isTRUE(save_df)) {
    saveRDS(df %>% select(-pval_cap, -p_adj, -idp_id), paste0(foldern, '/dataframe.sbxcan.rds'))
    saveRDS(df_cor_all %>% select(-id), paste0(foldern, '/dataframe.gencor.rds'))
  }
  
  if(isTRUE(gen_cor_qq_s)) {
    tmp = rbind(
      df %>% select(IDP, pval, idp_type, phenotype, model),
      df_cor_all %>% select(IDP, pval, idp_type, phenotype) %>% mutate(model = 'rg')
    )
    tmp = tmp %>% group_by(model, idp_type) %>% mutate(pval_cap = pmax(min_pval, pval)) %>% mutate(p_adj = pval * n()) %>% ungroup()
    tmp = tmp %>% group_by(model, idp_type) %>% mutate(pexp = rank(pval, ties.method = 'random') / (n() + 1)) %>% arrange(pval) 
    p = tmp %>% ungroup() %>% mutate(model = factor_methods(model)) %>% 
      ggplot() + 
      geom_path(aes(x = -log10(pexp), y = -log10(pval_cap), color = model), size = 1.5) + 
      # geom_point(data = tmp %>% filter(p_adj < 0.05), aes(x = -log10(pexp), y = -log10(pval_cap), color = model)) +
      facet_wrap(~idp_type) + 
      geom_abline(slope = 1, intercept = 0) + th2 +
      scale_color_manual(values = color_code) + 
      xlab(expression(paste(-log[10], p[expected]))) + 
      ylab(expression(paste(-log[10], p[observed]))) +
      theme(legend.position = c(0.9, 0.72), legend.title = element_blank(), legend.background = element_blank())
    ggsave(paste0(foldern, '/s_bxcan_qqplot_w_gencor.png'), p, width = 6, height = 3.5)
  }
  
  if(isTRUE(gen_cor_mr_s)) {
    df_cor_all = df_cor_all %>% group_by(idp_type) %>% mutate(pval_cap = pmax(min_pval, pval)) %>% mutate(p_adj = pval * n()) %>% ungroup()
    df_cor_sig = df_cor_all %>% filter(p_adj < alpha)
    for(ss in unique(df_cor_sig$folder)) {
      for(i in c('T1', 'dMRI')) {
        tmp = df_cor_sig %>% filter(idp_type == i, folder == ss) %>% mutate(pp = phenotype) %>% select(phenotype, IDP, pp) %>% rename(pheno = phenotype, idp = IDP, pheno_code = pp) %>% distinct()
        write.table(tmp, paste0(foldern, '/gencor_mr.', ss, '.', i, '.signif.tsv'), quo = F, col = T, row = F, sep = '\t')
      }
    }
  }
  
  if(isTRUE(gen_cor_compare)) {
    df_m = left_join(df, df_cor, by = c('phenotype', 'IDP', 'idp_type'), suffix = c('.bxcan', '.rg'))
    p = df_m %>% ggplot() + geom_point(aes(x = zscore.rg, y = zscore.bxcan, color = idp_type), alpha = 0.3) + 
      facet_wrap(~model) + 
      geom_abline(slope = 1, intercept = 0, color = 'gray') + th2 + 
      xlab('z-score of genetic correlation') + 
      ylab('z-score of S-BrainXcan') +
      scale_color_manual(values = color_code2) +
      theme(legend.position = c(0.6, 0.8), legend.title = element_blank())
    ggsave(paste0(foldern, '/s_bxcan_vs_gencor.png'), p, width = 5.5, height = 3)
  }
  
  if(isTRUE(gen_cor_check_mr)) {
    idp_meta = read.delim2('supp_table_1.tsv') %>% mutate(IDP = paste0('IDP-', ukb_field))
    mr_methods = c('Inverse variance weighted', 'Weighted median', 'MR Egger')
    dd_psych = rbind(
      read.table(paste0(foldern, '/gencor_mr.psychiatric.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
      read.table(paste0(foldern, '/gencor_mr.psychiatric.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
    )
    dd_gg = rbind(
      read.table(paste0(foldern, '/gencor_mr.gtex_gwas.T1.signif.tsv'), header = T) %>% mutate(idp_type = 't1'),
      read.table(paste0(foldern, '/gencor_mr.gtex_gwas.dMRI.signif.tsv'), header = T) %>% mutate(idp_type = 'dmri')
    )
    dd_gencor = rbind(dd_psych, dd_gg)
    psych_prefix = '~/Desktop/tmp/ukb_idp/mr_gencor_psych_3rd/MR_local.gencor_psych_'
    gg_prefix = '~/Desktop/tmp/ukb_idp/mr_gencor_gtexgwas_3rd/MR_local.gencor_gtexgwas_'
    df_mr = rbind(
      load_mr(dd_psych, df, prefix = psych_prefix) %>% mutate(source = 'psychiatric'),
      load_mr(dd_gg, df, prefix = gg_prefix) %>% mutate(source = 'gtex_gwas')
    )
    
    df_mr_entries_gencor = df_mr %>% select(phenotype, IDP, idp_type) %>% distinct()
    df_mr_entries_gencor = left_join(df_mr_entries_gencor, idp_meta %>% select(IDP, notes), by = 'IDP')
    message('distinct pairs passing the criteria: ', nrow(df_mr_entries_gencor))
    message('distinct GWASs passing the criteria: ', length(unique(df_mr_entries_gencor$phenotype)))
    message('distinct GWASs (psychiatric) passing the criteria: ', sum(!unique(df_mr_entries_gencor$phenotype) %in% not_psych))
  }
  
  if(isTRUE(gen_cor_check_mr) & isTRUE(check_s_mr) & isTRUE(plot_sig_mr)) {
    mm = full_join(
      df_mr_entries_gencor %>% select(-notes) %>% mutate(`Genetic Correlation` = T), 
      df_mr_entries %>% select(-notes) %>% mutate(`S-BrainXcan` = T), 
      by = c('phenotype', 'IDP', 'idp_type')
    )
    mmd = full_join(
      dd_sbxcan %>% select(-pheno_code) %>% mutate(`S-BrainXcan` = T),
      dd_gencor %>% select(-pheno_code) %>% mutate(`Genetic Correlation` = T), 
      by = c('pheno', 'idp', 'idp_type')
    )
    mmd = left_join(mmd, mm, by = c('pheno' = 'phenotype', 'idp' = 'IDP', 'idp_type'), suffix = c('.signif', '.mr'))
    mmd[is.na(mmd)] = F
    mmd[, 4:7] = as.matrix(mmd[, 4:7]) * 1
    png(paste0(foldern, '/', 's_bxcan_gencor_sig_mr_not_psych.png'), width = 7, height = 5, units = 'in', res = 300)
    UpSetR::upset(mmd %>% filter(pheno %in% not_psych), point.size = 5, line.size = 0.5, text.scale = 2)
    dev.off()
    png(paste0(foldern, '/', 's_bxcan_gencor_sig_mr_psych.png'), width = 7, height = 5, units = 'in', res = 300)
    UpSetR::upset(mmd %>% filter(!pheno %in% not_psych), point.size = 5, line.size = 0.5, text.scale = 2)
    dev.off()
  }
}

