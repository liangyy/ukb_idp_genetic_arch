# visualization test
library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(theme_classic(base_size = 10))
source('rmd/rlib.R')
traits = read.table('misc_data/preprocess_psychiatric_traits/trait_list.txt')$V1
models = c('ridge', 'elastic_net')
tags = list(t1 = 't1.scaled.all_covar.w_pc')
df = list()
for(rr in traits) {
  for(nn in names(tags)) {
    for(mm in models) {
      tmp = paste0('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/', tags[[nn]], '.gw_', mm, '_beta_x_', rr, '_x_simagexcan.csv')
      df[[length(df) + 1]] = read.csv(tmp) %>% 
        mutate(idp_type = nn, phenotype = rr, model = mm)
    }
  }
}

df = do.call(rbind, df)
df = df %>% mutate(zscore = p2z(pval, bhat))
# df = df %>% filter(is.na(stringr::str_match(IDP, 'PC')[, 1]))

tmp = df[df$phenotype == 'SCZ_PGC_2020', ]

tmp_ridge = tmp[ tmp$model == 'ridge', ]
vis_data = readRDS('misc_data/vis_data_t1.rds')
vis_data$table$IDP = paste0('IDP-', vis_data$table$FieldID)

# sub = vis_data$table[ vis_data$table$IDP %in% tmp_ridge$IDP,  ]

source('rmd/rlib_vis.R')

types = c('ho', 'cere', 'first')
scores = c('zscore', 'pip')
for(s in scores) {
  for(i in types) {
    message('Working on ', i, ' ', s)
    q = vis_assoc(vis_data, tmp_ridge, type = i, score = s)
    ggsave(paste0('misc_data/exploring_vis/', s, '_', i, '.png'), q, height = 4, width = 10)
  }
}


r1 = vis_region(vis_data, tmp_ridge, type = 'ho')
ggsave('misc_data/exploring_vis/region.ho.png', r1, height = 6, width = 13)
r2 = vis_region(vis_data, tmp_ridge, type = 'cere')
ggsave('misc_data/exploring_vis/region.cere.png', r2, height = 5, width = 10)
r3 = vis_region(vis_data, tmp_ridge, type = 'first')
ggsave('misc_data/exploring_vis/region.first.png', r3, height = 5, width = 10)

# # + geom_tile(aes(x = Var1, y = Var2, fill = value)) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') + scale_x_continuous()
# bb %>% filter(Var3 == mid3 / 2) %>% ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') + scale_x_continuous()
# tt %>% filter(Var1 == mid1 / 2) %>% ggplot() + geom_tile(aes(x = Var2, y = Var3, fill = value)) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
# tt %>% filter(Var2 == mid2 / 2) %>% ggplot() + geom_tile(aes(x = Var1, y = Var3, fill = value)) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
# 
# tt0 %>% filter(Var3 == mid3 / 2) %>% ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = position))
# tt0 %>% filter(Var1 == mid1 / 2) %>% ggplot() + geom_tile(aes(x = Var2, y = Var3, fill = position)) 
# tt0 %>% filter(Var2 == mid2 / 2) %>% ggplot() + geom_tile(aes(x = Var1, y = Var3, fill = position)) 
