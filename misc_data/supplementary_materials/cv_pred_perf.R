# setwd('misc_data/supplementary_materials/')

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

load_perf = function(wk_dir) {
  files = dir(wk_dir)
  df1 = list()
  for(f in files) {
    tmp = read.table(paste0(wk_dir, f), header = T)
    df1[[length(df1) + 1]] = tmp
  }
  df1 = do.call(rbind, df1)
  df1$phenotype = stringr::str_replace(df1$phenotype, 'x', '-')
  df1
}
# elastic net
df1 = load_perf('~/Desktop/tmp/ukb_idp/gw_elastic_net_t1_2nd/t1.scaled.all_covar.w_pc/')
df2 = load_perf('~/Desktop/tmp/ukb_idp/gw_elastic_net_dmri_2nd/dmri.original.all_covar.w_pc/')
# ridge 
df1r = read.table('~/Desktop/tmp/ukb_idp/gw_ridge_2nd/t1.scaled.all_covar.w_pc.perf.tsv.gz', header = T)
df2r = read.table('~/Desktop/tmp/ukb_idp/gw_ridge_2nd/dmri.original.all_covar.w_pc.perf.tsv.gz', header = T)
df1 = rbind(df1 %>% mutate(model = 'EN'), df1r %>% mutate(model = 'ridge')) %>%
  mutate(is_pc = substr(phenotype, 1, 2) == 'PC')
df2 = rbind(df2 %>% mutate(model = 'EN'), df2r %>% mutate(model = 'ridge')) %>%
  mutate(is_pc = substr(phenotype, 1, 2) == 'PC')

p = rbind(df2 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), df1 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))) %>% 
  mutate(PC = as.character(is_pc)) %>% 
  reshape2::dcast(phenotype + idp_type + PC ~ model, value.var = 'Spearman') %>% 
  ggplot() + geom_point(aes(x = ridge, y = EN, color = PC)) +
  facet_wrap(~idp_type) + th2 + geom_abline(intercept = 0, slope = 1) + coord_equal() +
  theme(legend.position = c(0.15, 0.8), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) 
ggsave('cv_pred_perf.png', p, width = 5.5, height = 4)


# add heritability

df2h = read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/dmri.original.all_covar.w_pc.tsv.gz', header = T)
df1h = read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/t1.scaled.all_covar.w_pc.tsv.gz', header = T)
df1h$is_pc = substr(df1h$phenotype, 1, 2) == 'PC'
df2h$is_pc = substr(df2h$phenotype, 1, 2) == 'PC'

tmp = rbind(
  inner_join(
    df1 %>% filter(model == 'ridge') %>% select(-is_pc), 
    df1h, by = c('phenotype')
  ) %>% 
    mutate(idp_type = factor('T1', levels = c('T1', 'dMRI'))), 
  inner_join(
    df2 %>% filter(model == 'ridge') %>% select(-is_pc), 
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
  geom_point(aes(x = h2, y = sign_pred_perf, color = is_pc)) +
    theme(legend.position = c(0.15, 0.75)) +
    scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
    th2 + facet_wrap(~idp_type) + ylab('Signed squared \n Spearman correlation') 
ggsave('cv_pred_perf_vs_h2.png', p, width = 7, height = 3.5)
