library(dplyr)
library(ggplot2)
library(patchwork)


second_dir = '../second_round_idp_preprocessing'
source(paste0(second_dir, '/scripts/rlib.R'))
source('rlib.R')

# fig size
ww = 12
hh = 10
# END

runall = T
if(isTRUE(runall)) {
  
  outdir = 'explore_25901'
  dir.create(outdir)
  output_tags = c('t1.scaled', 't1.original', 't1.scaled.all_covar', 'raw')
  
  # explore 25901 
  d0 = arrow::read_parquet(paste0(second_dir, '/output/t1.scaled.parquet'))
  d1 = arrow::read_parquet(paste0(second_dir, '/output/t1.original.parquet'))
  d2 = arrow::read_parquet(paste0(second_dir, '/output/t1.scaled.all_covar.parquet'))
  df_idp = data.table::fread('~/Desktop/tmp/ukb_idp/idp_phenotypes/archive/2020-05-07_idp-phenotypes_cleaned.txt', sep = '\t', data.table = F)
  
  idps = read.delim2('../supplementary_materials/supp_table_1.tsv')
  tag = 'Cerebellum'
  
  kk = list()
  for(t1_mat in list(d0, d1, d2, df_idp)) {
    idps_cort = idps %>% filter(t1_or_dmri == 'T1', t1_anatomy_group == tag)
    cort_all = as.matrix(t1_mat[, paste0('IDP-', idps_cort$ukb_field)])
    kk[[length(kk) + 1]] = do_all(cort_all)
  }
  ggsave(paste0(outdir, '/', output_tags[1], '.png'), plot_all(kk[[1]]$ps, 't1.scaled'), width = ww, height = hh)
  ggsave(paste0(outdir, '/', output_tags[2], '.png'), plot_all(kk[[2]]$ps, 't1.original'), width = ww, height = hh)
  ggsave(paste0(outdir, '/', output_tags[3], '.png'), plot_all(kk[[3]]$ps, 't1.scaled.all_covar'), width = ww, height = hh)
  ggsave(paste0(outdir, '/', output_tags[4], '.png'), plot_all(kk[[4]]$ps, 'phenotype matrix'), width = ww, height = hh)
  
}



