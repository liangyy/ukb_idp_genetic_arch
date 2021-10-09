# setwd('misc_data/supplementary_materials_4th/')

library(dplyr)

# look into scz
df = rbind(
  read.csv('~/Desktop/tmp/ukb_idp/brainxcan_pipeline/psychiatric/SCZ_PGC_2020.sbrainxcan.csv') %>% 
    mutate(phenotype = 'SCZ_PGC_2020', type = 'residual'),
  read.csv('~/Desktop/tmp/ukb_idp/brainxcan_pipeline/gtex_gwas/pgc.scz2.sbrainxcan.csv') %>% 
    mutate(phenotype = 'pgc.scz2', type = 'residual'),
  read.csv('~/Desktop/tmp/ukb_idp/brainxcan_pipeline/original/pgc.scz2.sbrainxcan.csv') %>% 
    mutate(phenotype = 'pgc.scz2', type = 'original'),
  read.csv('~/Desktop/tmp/ukb_idp/brainxcan_pipeline/original/SCZ_PGC_2020.sbrainxcan.csv') %>% 
    mutate(phenotype = 'SCZ_PGC_2020', type = 'original')
)

myfunc = function(idp) {
  df_hip = df %>% filter(IDP == idp)
  mm = readRDS('s_bxcan/dataframe_full.sbxcan.rds') %>% 
    filter(phenotype %in% df_hip$phenotype) %>% 
    filter(IDP == df_hip$IDP[1]) %>% 
    filter(model == 'elastic net') %>% 
    mutate(type = 'residual', notes = df_hip$notes[1]) 
  rbind(
    df_hip %>% select(IDP, phenotype, bhat, pval, type, notes) %>% mutate(model = 'ridge'),
    mm %>% select(IDP, phenotype, bhat, pval, type, notes, model)
  ) %>% select(-notes) %>% arrange(model, type, phenotype)
}
tmp = myfunc('IDP-25020')


df_amygdala = df %>% filter(region == 'amygdala')
tmp = myfunc('IDP-25888'); tmp


df_thalamus = df %>% filter(region == 'thalamus')
tmp = myfunc('IDP-25012'); tmp

