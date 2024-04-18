# misc statistics for BrainXcan results
library(ggplot2)
library(dplyr)
source('rlib.R')
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

outdir = 'misc_summary_bxcan_twasinf'
dir.create(outdir)

kk = readRDS('~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/indiv_bxcan/dataframe_full.indiv_bxcan.rds') %>% filter(model == 'ridge')
kk2 = readRDS('s_bxcan_twasinf/dataframe_full.sbxcan_twasinf.rds') %>% filter(model == 'ridge')
kk3 = readRDS('~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/s_bxcan/dataframe_full.gencor.rds')
df_all = rbind(
  kk %>% select(IDP, phenotype, bhat, pval, zscore, idp_type, model, idp_type) %>% mutate(source = 'indiv_bxcan'),
  kk2 %>% select(IDP, phenotype, bhat, pval_corrected, z_corrected, modality, model, source) %>%
    rename(pval = pval_corrected, zscore = z_corrected, idp_type = modality)
) 


# remove ProbTrack IDPs
df_all = remove_probtrack_idp(df_all)

# limit to IDPs with perf > 0.1
filter_pred_perf <- function(df) {
  df = df %>% mutate(idp_id = paste(idp_type, IDP, model))
  idp_sig = read.table('~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv', header = T, sep = '\t')
  idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
  df = df[df$idp_id %in% idp_sig$idp_id, ]
  df %>% select(-idp_id)
}
df_all = filter_pred_perf(df_all)

# limit gen cor results on the same set of IDPs and format a bit
kk3 = kk3 %>% 
  filter(IDP %in% df_all$IDP) %>% 
  select(IDP, phenotype, rg, pval, zscore, idp_type) %>% 
  mutate(model = 'genetic_cor', source = 's_bxcan') %>% 
  rename(bhat = rg)
df_all = rbind(
  df_all,
  kk3
)

if(FALSE) {
  # DEPRECATED: significance criteria: per phenotype, p < 0.05 under Bonferroni-correction
  # NEW: signficance criteria: per phenotype, qvalue < 0.05
  alpha = 0.05
  
  # set p adjusted
  df_all %>% group_by(phenotype, source, model) %>% summarise(ntest = n()) %>% ggplot() + geom_point(aes(x = phenotype, y = ntest)) + facet_grid(source ~ model)
  get_qvalue <- function(pval) {
    maxp <- max(pval)
    lambda.seq <- seq(0.05, 0.95, 0.05)
    if(maxp < 0.95) {
      lambda.seq <- lambda.seq[lambda.seq <= maxp]
    }
    return(qvalue::qvalue(pval, lambda = lambda.seq)$qvalues)
  }
  df_all = df_all %>% group_by(phenotype, source, model) %>% mutate(p_adj = get_qvalue(pval)) %>% ungroup()
  idp_signif = df_all %>% filter(model == 'ridge', p_adj < alpha) %>% pull(IDP)
  message(length(unique(idp_signif)), ' out of ', length(unique(df_all$IDP)), ' IDPs are significant in at least on phenotype.')
  
  idp_signif = df_all %>% filter(model == 'ridge', source != 'indiv_bxcan', p_adj < alpha) %>% pull(IDP)
  message(length(unique(idp_signif)), ' out of ', length(unique(df_all$IDP)), ' IDPs are significant in at least on phenotype.')
  
  # correlation between gen cor and ridge
  tmp = df_all %>% filter(source != 'indiv_bxcan')
  maxz = max(abs(tmp$zscore[!is.infinite(tmp$zscore)]))
  tmp$zscore[is.infinite(tmp$zscore)] = sign(tmp$bhat[is.infinite(tmp$zscore)]) * maxz
  df_cor = list()
  for(p in unique(tmp$phenotype)) {
    kk = tmp %>% filter(phenotype == p) %>% reshape2::dcast(IDP ~ model, value.var = 'zscore')
    df_cor[[length(df_cor) + 1]] = data.frame(cor = cor(kk$genetic_cor, kk$ridge), phenotype = p)
  }
  df_cor = do.call(rbind, df_cor)
  df_cor %>% summarize(max_cor = max(cor), median_cor = median(cor), min_cor = min(cor))
  
  # nsig in ridge and gen cor
  tmp = df_all %>% filter(source != 'indiv_bxcan')
  tmp %>% mutate(signif = p_adj < alpha) %>% group_by(model) %>% summarize(nsignif = sum(signif), ntotal = n()) %>% mutate(fold = max(nsignif) / nsignif)
  
  # SCZ
  
  df_scz = kk2 = readRDS('s_bxcan_permz/dataframe_full.sbxcan_permz.rds') %>% filter(model == 'ridge') %>% filter(phenotype == 'SCZ_PGC_2020') %>% select(-pval) %>% rename(idp_type = modality, pval = pval_adj_perm_null)
  annot = load_idp_annot()
  df_scz = filter_pred_perf(df_scz)
  df_scz_sub = remove_probtrack_idp(df_scz)
  df_scz %>% group_by(model) %>% summarize(ntest = n())
  df_scz_sub %>% group_by(model) %>% summarize(ntest = n())
  # df_scz_sub = df_scz_sub %>% left_join(annot, by = 'IDP')
  df_scz_sub %>% mutate(is_pc = substr(IDP, 1, 2) == 'PC') %>% group_by(is_pc, subtype) %>% summarize(ntest = n()) 
  df_scz_sub = df_scz_sub %>% mutate(p_adj = get_qvalue(pval)) 
  df_scz_sub %>% filter(p_adj < alpha) %>% nrow()
  
  
  # side analysis 
  common_pairs = list(
    Height = c('UKB_50_Standing_height', 'GIANT_HEIGHT'),
    BMI = c('UKB_21001_Body_mass_index_BMI', 'GIANT_2015_BMI_EUR'),
    Intelligence = c('Intelligence_CTG', 'UKB_20016_Fluid_intelligence_score')
  )  
  df_cor_examples = list()
  for(nn in names(common_pairs)) {
    mm = df_cor %>% filter(phenotype %in% common_pairs[[nn]])
    mm = mm[ match(common_pairs[[nn]], mm$phenotype), ]
    df_cor_examples[[length(df_cor_examples) + 1]] = data.frame(rel_low = mm$cor[1], rel_high = mm$cor[2], phenotype = nn)
  }
  df_cor_examples = do.call(rbind, df_cor_examples)
}



# zscore vs h2
# load h2
df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
df_h2 = rbind(df1, df2) %>% select(phenotype, h2, h2_SE) %>% rename(IDP = phenotype)
df_perf = read.delim2('~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv')

height_and_bmi = c('GIANT_HEIGHT', 'bmi', 'height', 'UKB_50_Standing_height', 'GIANT_2017_BMI_Active_EUR', 'GIANT_2015_BMI_EUR', 'UKB_21001_Body_mass_index_BMI')
psy = read.table('../preprocess_psychiatric_traits/trait_list.txt')
df_ridge = df_all %>% filter(model == 'ridge') # %>% filter(!phenotype %in% height_and_bmi) %>% filter(phenotype %in% psy$V1)
df_ridge = df_ridge %>% left_join(df_h2, by = 'IDP') %>% left_join(df_perf %>% filter(model_name == 'ridge'), by = 'IDP') %>% mutate(Spearman = as.numeric(Spearman))
breaks = as.numeric(quantile(df_ridge$Spearman, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
breaks[1] = breaks[1] - 1e-10
breaks[length(breaks)] = breaks[length(breaks)] + 1e-10
p = df_ridge %>% 
  filter(source != 'indiv_bxcan') %>%
  mutate(rank_z = rank(abs(zscore))) %>% 
  mutate(hbin = cut(Spearman, breaks = breaks)) %>%
  ggplot() +
  geom_violin(aes(x = hbin, y = abs(zscore), group = hbin)) + 
  geom_boxplot(aes(x = hbin, y = abs(zscore), group = hbin), width = 0.1) + th +
  geom_hline(yintercept = median(abs(df_ridge$zscore)), linetype = 2) + 
  coord_cartesian(ylim = c(0, 10)) + 
  xlab('Spearman correlation bin') +
  ylab('Absolute value of BrainXcan z-score') +
  ggtitle('Magnitude of BrainXcan z-score increases with \n IDP prediction performance')
ggsave(paste0(outdir, '/perf_vs_zscore.png'), p, height = 5, width = 6)

annot_pc = function(x) {
  o = substr(x, 1, 2) == 'PC'
  kk = rep('Region-Specific', length(o))
  kk[o] = 'Common Factor'
  factor(kk, levels = c('Region-Specific', 'Common Factor'))
}
p = df_ridge %>%
  filter(source != 'indiv_bxcan') %>%
  mutate(is_pc = annot_pc(IDP)) %>% ggplot() +
  geom_violin(aes(x = is_pc, y = abs(zscore))) + 
  geom_boxplot(aes(x = is_pc, y = abs(zscore)), width = 0.1) + th2 + 
  coord_cartesian(ylim = c(0, 10)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12)) +
  ylab('Absolute value of BrainXcan z-score') +
  ggtitle('PCs (common factors) \n yield more significant BrainXcan z-scores')
ggsave(paste0(outdir, '/pc_vs_zscore.png'), p, height = 5, width = 5)
