# visualization test
library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(theme_classic(base_size = 10))
source('rmd/rlib.R')
traits = read.table('misc_data/preprocess_psychiatric_traits/trait_list.txt')$V1
models = c('ridge', 'elastic_net')
tags = list(dmri = 'dmri.original.all_covar.w_pc')
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

dmri = readRDS('misc_data/download_some_matching_files/annot_dmri_idps.rds') %>% mutate(IDP = paste0('IDP-', FieldID)) 

tmp_ridge_fa = tmp[ tmp$model == 'elastic_net' & tmp$IDP %in% (dmri %>% filter(measure == 'FA') %>% pull(IDP)), ]
tmp_ridge_icvf = tmp[ tmp$model == 'elastic_net' & tmp$IDP %in% (dmri %>% filter(measure == 'ICVF') %>% pull(IDP)), ]
vis_data = readRDS('misc_data/vis_data_dmri.rds')
vis_data$table$IDP = paste0('IDP-', vis_data$table$FieldID)

# sub = vis_data$table[ vis_data$table$IDP %in% tmp_ridge_fa$IDP,  ]

source('rmd/rlib_vis.R')

types = c('tbss')
scores = c('zscore', 'pip')
for(s in scores) {
  for(i in types) {
    message('Working on ', i, ' ', s)
    q = vis_assoc(vis_data, tmp_ridge_icvf, type = i, score = s)
    ggsave(paste0('misc_data/exploring_vis/', s, '_', i, '.icvf.png'), q, height = 4, width = 10)
  }
}


r1 = vis_region(vis_data, tmp_ridge_icvf, type = 'tbss')
ggsave('misc_data/exploring_vis/region.tbss.icvf.png', r1, height = 6, width = 13)

