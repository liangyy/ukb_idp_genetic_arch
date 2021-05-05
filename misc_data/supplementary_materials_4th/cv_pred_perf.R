# setwd('misc_data/supplementary_materials_4th/')

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
options(stringsAsFactors = F)

foldern = 'cv_pred_perf/'
dir.create(foldern)

color_vec = c('Common Factor' = 'red', 'Region-Specific' = 'black')


# elastic net
df1 = read.table('~/Desktop/tmp/ukb_idp/idp_models_4th/fourth_round.t1_w_pc.gw_elastic_net_beta.perf.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/idp_models_4th/fourth_round.dmri_w_pc.gw_elastic_net_beta.perf.tsv.gz', header = T)
# ridge 
df1r = read.table('~/Desktop/tmp/ukb_idp/idp_models_4th/fourth_round.t1_w_pc.gw_ridge_beta.perf.tsv.gz', header = T)
df2r = read.table('~/Desktop/tmp/ukb_idp/idp_models_4th/fourth_round.dmri_w_pc.gw_ridge_beta.perf.tsv.gz', header = T)
df1 = rbind(df1 %>% mutate(model = 'EN'), df1r %>% mutate(model = 'ridge')) %>%
  mutate(is_pc = substr(phenotype, 1, 2) == 'PC')
df2 = rbind(df2 %>% mutate(model = 'EN'), df2r %>% mutate(model = 'ridge')) %>%
  mutate(is_pc = substr(phenotype, 1, 2) == 'PC')
df2$pc = rep('Region-Specific', nrow(df2))
df2$pc[df2$is_pc] = 'Common Factor'
df1$pc = rep('Region-Specific', nrow(df1))
df1$pc[df1$is_pc] = 'Common Factor'

p = rbind(df2 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), df1 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))) %>% 
  mutate(PC = as.character(pc)) %>% 
  reshape2::dcast(phenotype + idp_type + PC ~ model, value.var = 'Spearman') %>% 
  ggplot() + geom_point(aes(x = ridge, y = EN, color = PC)) +
  facet_wrap(~idp_type) + th2 + geom_abline(intercept = 0, slope = 1) + coord_equal() +
  theme(legend.position = c(0.17, 0.85), legend.text = element_text(size = 9), legend.title = element_blank()) +
  scale_color_manual(values = color_vec) + 
  xlab('Prediction performance of ridge regression') + 
  ylab('Prediction performance of elastic net')
ggsave(paste0(foldern, 'cv_pred_perf_spearman.png'), p, width = 5.5, height = 4)

p = rbind(df2 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), df1 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))) %>% 
  mutate(PC = as.character(pc)) %>% 
  reshape2::dcast(phenotype + idp_type + PC ~ model, value.var = 'Pearson') %>% 
  ggplot() + geom_point(aes(x = ridge, y = EN, color = PC)) +
  facet_wrap(~idp_type) + th2 + geom_abline(intercept = 0, slope = 1) + coord_equal() +
  theme(legend.position = c(0.17, 0.862), legend.text = element_text(size = 9), legend.title = element_blank()) +
  scale_color_manual(values = color_vec) + 
  xlab('Prediction performance of ridge regression') + 
  ylab('Prediction performance of elastic net')
ggsave(paste0(foldern, 'cv_pred_perf_pearson.png'), p, width = 5.5, height = 4)

p = rbind(df2 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), df1 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))) %>% 
  mutate(PC = as.character(pc)) %>% 
  reshape2::dcast(phenotype + idp_type + PC ~ model, value.var = 'R2') %>% 
  ggplot() + geom_point(aes(x = ridge, y = EN, color = PC)) +
  facet_wrap(~idp_type) + th2 + geom_abline(intercept = 0, slope = 1) + coord_equal() +
  theme(legend.position = c(0.17, 0.825), legend.text = element_text(size = 9), legend.title = element_blank()) +
  scale_color_manual(values = color_vec) + 
  xlab('Prediction performance of ridge regression') + 
  ylab('Prediction performance of elastic net')
ggsave(paste0(foldern, 'cv_pred_perf_r2.png'), p, width = 5.5, height = 4)


# add heritability

df2h = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
df1h = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
df1h$is_pc = substr(df1h$phenotype, 1, 2) == 'PC'
df2h$is_pc = substr(df2h$phenotype, 1, 2) == 'PC'
df2h$pc = rep('Region-Specific', nrow(df2h))
df2h$pc[df2h$is_pc] = 'Common Factor'
df1h$pc = rep('Region-Specific', nrow(df1h))
df1h$pc[df1h$is_pc] = 'Common Factor'

tmp = rbind(
  inner_join(
    df1 %>% filter(model == 'ridge') %>% select(-pc), 
    df1h, by = c('phenotype')
  ) %>% 
    mutate(idp_type = factor('T1', levels = c('T1', 'dMRI'))), 
  inner_join(
    df2 %>% filter(model == 'ridge') %>% select(-pc), 
    df2h, by = c('phenotype')
  ) %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI')))
) 
# tmp %>% mutate(pheno = factor(paste(idp_type, phenotype), levels = paste(idp_type, phenotype)[order(idp_type, h2)]), PC = as.character(is_pc)) %>% 
#   ggplot() + 
#   geom_errorbar(aes(x = pheno, ymax = h2 + 1.96 * h2_SE, ymin = h2 - 1.96 * h2_SE), color = 'gray') + 
#   geom_point(aes(x = pheno, y = h2, color = PC)) +
#   geom_point(aes(x = pheno, y = sign(Pearson) * Pearson^2, color = PC)) +
#   theme(axis.text.x = element_blank(), legend.position = c(0.4, 0.8)) +
#   scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
#   th2 + 
#   facet_grid(.~idp_type, scales = 'free_x', space = "free_x")
  
  
p = tmp %>% mutate(sign_pred_perf = sign(Spearman) * (Spearman ^ 2)) %>% ggplot() +
  geom_errorbarh(aes(y = sign_pred_perf, xmax = h2 + 1.96 * h2_SE, xmin = h2 - 1.96 * h2_SE), color = 'gray') + 
  geom_point(aes(x = h2, y = sign_pred_perf, color = pc)) +
  theme(legend.position = c(0.15, 0.8), legend.title = element_blank()) +
  scale_color_manual(values = color_vec) +
  th2 + facet_wrap(~idp_type) + ylab('Signed squared \n Spearman correlation') +
  xlab('Heritability')
ggsave(paste0(foldern, 'cv_pred_perf_vs_h2_ridge.png'), p, width = 7, height = 3.5)

tmp = rbind(
  inner_join(
    df1 %>% filter(model == 'EN') %>% select(-pc), 
    df1h, by = c('phenotype')
  ) %>% 
    mutate(idp_type = factor('T1', levels = c('T1', 'dMRI'))), 
  inner_join(
    df2 %>% filter(model == 'EN') %>% select(-pc), 
    df2h, by = c('phenotype')
  ) %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI')))
) 

p = tmp %>% mutate(sign_pred_perf = sign(Spearman) * (Spearman ^ 2)) %>% ggplot() +
  geom_errorbarh(aes(y = sign_pred_perf, xmax = h2 + 1.96 * h2_SE, xmin = h2 - 1.96 * h2_SE), color = 'gray') + 
  geom_point(aes(x = h2, y = sign_pred_perf, color = pc)) +
  theme(legend.position = c(0.15, 0.8), legend.title = element_blank()) +
  scale_color_manual(values = color_vec) +
  th2 + facet_wrap(~idp_type) + ylab('Signed squared \n Spearman correlation') +
  xlab('Heritability')
ggsave(paste0(foldern, 'cv_pred_perf_vs_h2_en.png'), p, width = 7, height = 3.5)


dd = rbind(df1 %>% mutate(idp_type = 'T1'), df2 %>% mutate(idp_type = 'dMRI'))
p1 = dd %>% filter(model == 'ridge') %>% select(-model, -is_pc, -pc) %>% reshape2::melt(id.var = c('phenotype', 'idp_type')) %>% ggplot() + geom_histogram(aes(x = value), bins = 100) + facet_wrap(~variable, scales = 'free', ncol = 1) + geom_vline(xintercept = 0) + th2 + xlab('Prediction performance') +
  theme(plot.margin = unit(c(.1,1,.1,1), "cm"))

p2 = dd %>% filter(model == 'EN') %>% select(-model, -is_pc, -pc) %>% reshape2::melt(id.var = c('phenotype', 'idp_type')) %>% ggplot() + geom_histogram(aes(x = value), bins = 100) + facet_wrap(~variable, scales = 'free', ncol = 1) + geom_vline(xintercept = 0) + th2 + xlab('Prediction performance') +
  theme(plot.margin = unit(c(.1,1,.1,1), "cm"))

ggsave(paste0(foldern, 'cv_pred_perf_hist_ridge.png'), p1, width = 6, height = 4.5)
ggsave(paste0(foldern, 'cv_pred_perf_hist_en.png'), p2, width = 6, height = 4.5)

dd %>% group_by(model) %>% summarize(mean(Spearman > 0), mean(Pearson > 0), mean(R2 > 0))
dd$model_name = rep('ridge', nrow(dd))
dd$model_name[dd$model == 'EN'] = 'elastic net'
dd %>% select(-is_pc, -pc, -model) %>% select(idp_type, phenotype, model_name, Spearman, Pearson, R2) %>% 
  rename(IDP = phenotype, IDP_type = idp_type) %>% 
  mutate(is_kept = Spearman > 0.1) %>% 
  write.table('supp_table_2.tsv', quote = F, sep = '\t', row.names = F, col.names = T)

idp_sig = read.table('supp_table_2.tsv', header = T, sep = '\t')
idp_sig %>% filter(is_kept) %>% group_by(IDP_type, model_name) %>% summarize(count = n())

library(UpSetR)
tmp = mat = idp_sig %>% reshape2::dcast(IDP + IDP_type ~ model_name, value.var = 'is_kept') %>% mutate(total = T)
mat = as.matrix(mat[ , c(-1, -2)])
mat_ = mat
mat_[mat] = 1
mat_[!mat] = 0
tmp[, 3:5] = mat_
t1_col = 'orange'  # rgb()
dmri_col = 'blue'  # rgb()

png(paste0(foldern, '/', 'kept_models_t1.png'), width = 7, height = 5, units = 'in', res = 300)
upset(tmp %>% filter(IDP_type == 'T1'), main.bar.color = t1_col, point.size = 5, line.size = 0.5, text.scale = 2)
dev.off()

png(paste0(foldern, '/', 'kept_models_dmri.png'), width = 7, height = 5, units = 'in', res = 300)
upset(tmp %>% filter(IDP_type == 'dMRI'), main.bar.color = dmri_col, point.size = 5, line.size = 0.5, text.scale = 2)
dev.off()
