rm(list = ls())
set.seed(2021)
library(ggplot2)
library(patchwork)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/489d1519dd45246caf4756d7722bfa25/raw/9bbb39b80243325b7930083063566fae4af85d48/fast_linear_regression')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
df1 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.no_pc.parquet')
df2 = arrow::read_parquet('../second_round_idp_preprocessing/output/t1.scaled.all_covar.no_pc.parquet')
df12 = arrow::read_parquet('../second_round_idp_preprocessing/output/dmri.original.all_covar.w_pc.parquet')
df22 = arrow::read_parquet('../second_round_idp_preprocessing/output/t1.scaled.all_covar.w_pc.parquet')

sim_func = function(df, nrepeat = 10) {
  tmp = list()
  for(i in 1 : nrepeat) {
    index = sample(ncol(df) - 1, 2, replace = F)
    d1 = index[1] + 1
    d2 = index[2] + 1
    label = rep(F, ncol(df) - 1)
    label[d1 - 1] = T
    label[d2 - 1] = T
    y = df[[d1]] * 1 + df[[d2]] * 1 + rnorm(nrow(df), sd = 20) #  
    pp = fast_linear_regression(y, df[, -1])
    tmp[[length(tmp) + 1]] = data.frame(pval = pp$pval, label = label, idx = i)
  }
  tmp = do.call(rbind, tmp)
  
  tmp = tmp %>% group_by(idx) %>%
    mutate(pexp = rank(pval) / (n() + 1)) %>% ungroup()
  
  p = tmp %>% 
    ggplot() + 
    geom_line(aes(x = -log10(pexp), y = -log10(pval), group = idx), alpha = 0.5) +
    geom_point(data = tmp %>% filter(label), aes(x = -log10(pexp), y = -log10(pval))) +
    geom_abline(slope = 1, intercept = 0) + th +
    theme(legend.position = c(0.2, 0.8)) +
    geom_hline(yintercept = -log10(0.05 / (ncol(df) - 1)), linetype = 2) # +
  # geom_hline(data = tmp %>% filter(label), aes(yintercept = -log10(pval)))
  p
}
# nrepeat = 10
p1 = sim_func(df1, nrepeat = 5) + ggtitle('dMRI: Without adjustment')
p2 = sim_func(df2, nrepeat = 5) + ggtitle('T1: Without adjustment')
p12 = sim_func(df12, nrepeat = 5) + ggtitle('dMRI: With adjustment')
p22 = sim_func(df22, nrepeat = 5) + ggtitle('T1: With adjustment')

# (p1 + p12)

ggsave('sim_cor_dmri.png', p1, width = 5, height = 5)
ggsave('sim_cor_t1.png', p2, width = 5, height = 5)
ggsave('sim_cor_dmri_adj.png', p12, width = 5, height = 5)
ggsave('sim_cor_t1_adj.png', p22, width = 5, height = 5)
# mod = susieR::susie(as.matrix(df[, -1]), y)
# plot(mod$pip)
# summary(mod)$cs
