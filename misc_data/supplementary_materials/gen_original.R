# setwd('misc_data/supplementary_materials/')
# generate data from for data in original scale


load_sbxcan = function(folder, trait_list, idp_type) {
  models = list(ridge = 'gw_ridge', EN = 'gw_elastic_net')
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

library(dplyr)
library(patchwork)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('../../rmd/rlib.R')
source('../../rmd/rlib_calc.R')

foldern = 'gen_original'
dir.create(foldern)

save_df_full = T

not_psych = c('GIANT_2015_BMI_EUR', 'GIANT_2017_BMI_Active_EUR', 'GIANT_HEIGHT', 'UKB_50_Standing_height', 'UKB_21001_Body_mass_index_BMI')
color_code = c('ridge' = 'blue', 'elastic net' = 'orange', 'rg' = 'pink')
color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
t1_col = 'orange'  
dmri_col = 'blue'  
factor_idp = function(cc) {
  factor(cc, levels = c('T1', 'dMRI'))
}

df_gwas = read.delim2('supp_table_4.tsv', header = T)
folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_2nd')
df_gwas$folder = 'gtex_gwas'
df_gwas$folder[25:35] = 'psychiatric'
df = list()

idp_type = list(
  dMRI = 'dmri.original.all_covar.no_pc', 
  T1 = 't1.scaled.all_covar.no_pc'
)

for(cc in c('gtex_gwas')) {
  traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
  df[[length(df) + 1]] = load_sbxcan(folders[[cc]], traits, idp_type) %>% mutate(source = cc)
}
df = do.call(rbind, df)
df$model[df$model == 'EN'] = 'elastic net'
df = df %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))
if(isTRUE(save_df_full)) {
  saveRDS(df %>% select(-idp_id), paste0(foldern, '/dmri_original_t1_scaled_all_covar_no_pc.rds'))
  # saveRDS(df_cor, paste0(foldern, '/dataframe_full.gencor.rds'))
}

