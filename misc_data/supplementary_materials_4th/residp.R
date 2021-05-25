# resIDP vs IDP adjusted by PC

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('rlib.R')

outdir = 'residp'
dir.create(outdir)

# load IDP adjusted by PC (S-BrainXcan)
{
  fn = paste0(outdir, '/dataframe_full.sbxcan_w_PCadj.rds')
  if(file.exists(fn)) {
    kk2 = readRDS(fn)
  } else {
    load_sbxcan = function(folder, trait_list) {
      idp_type = list(dMRI = 'dmri', T1 = 't1')
      models = list(ridge = 'ridge', EN = 'en')
      df1 = list()
      for(t in trait_list) {
        for(i in names(idp_type)) {
          for(m in names(models)) {
            df1[[length(df1) + 1]] = read.csv(paste0(folder, '_', models[[m]], '/', idp_type[[i]], '_x_', t, '_x_simagexcan.csv'), header = T) %>% mutate(idp_type = i, model = m, phenotype = t)
          }
        }
      }
      df1 = do.call(rbind, df1)
      df1
    }
    
    source('../../rmd/rlib.R')
    source('../../rmd/rlib_calc.R')
    
    not_psych = c('GIANT_2015_BMI_EUR', 'GIANT_2017_BMI_Active_EUR', 'GIANT_HEIGHT', 'UKB_50_Standing_height', 'UKB_21001_Body_mass_index_BMI')
    color_code = c('ridge' = 'blue', 'elastic net' = 'orange', 'rg' = 'pink')
    color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')
    t1_col = 'orange'  
    dmri_col = 'blue'  
    factor_idp = function(cc) {
      factor(cc, levels = c('T1', 'dMRI'))
    }
    
    df_gwas = read.delim2('../supplementary_materials_3rd/supp_table_4.tsv', header = T)
    folders = list(gtex_gwas = '~/Desktop/tmp/ukb_idp/simagexcan/results_gtex_gwas_4th', psychiatric = '~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_4th')
    df_gwas$folder = 'gtex_gwas'
    df_gwas$folder[25:35] = 'psychiatric'
    df = list()
    for(cc in c('gtex_gwas', 'psychiatric')) {
      traits = df_gwas %>% filter(folder == cc) %>% pull(phenotype_id)
      df[[length(df) + 1]] = load_sbxcan(folders[[cc]], traits) %>% mutate(source = cc)
    }
    df = do.call(rbind, df)
    df$model[df$model == 'EN'] = 'elastic net'
    df = df %>% filter(substr(test, 1, 9) == 'adj_covar')
    df = df %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))
    saveRDS(df %>% select(-idp_id), fn)
    kk2 = df
  }
}

kk = readRDS('s_bxcan/dataframe_full.sbxcan.rds')
kk = remove_probtrack_idp(kk)
kk = filter_pred_perf(kk)

tmp1 = inner_join(
  kk %>% select(IDP, model, phenotype, zscore), 
  kk2 %>% select(IDP, model, phenotype, zscore), 
  by = c('IDP', 'model', 'phenotype'), suffix = c('.resIDP', '.PCadj')
)
# p = tmp1 %>% ggplot() +
#   geom_point(aes(x = zscore.resIDP, y = zscore.PCadj), alpha = 0.3) + 
#   geom_abline(slope = 1, intercept = 0, color = 'blue') + 
#   coord_equal() + 
#   xlab('S-BrainXcan zscore of residual IDP') + 
#   ylab('S-BrainXcan zscore of IDP \n when fitting with PC jointly') + facet_wrap(~model) + 
#   th2
# ggsave(paste0(outdir, '/residp_vs_idp_pcadj.png'), p, width = 5, height = 3)

# load IDP adjusted by PC (indiv-BrainXcan)
{
  fn = paste0(outdir, '/dataframe_full.indiv_bxcan_w_PCadj.rds')
  if(file.exists(fn)) {
    kk4 = readRDS(fn)
  } else {
    
    pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes', 'bmi', 'height')
    models = list(ridge = 'gw_ridge_beta', EN = 'gw_elastic_net_beta')
    idps = list(dMRI = 'dmri', T1 = 't1') #  
    # types = c('linear', 'susie')
    df = list()
    for(m in names(models)) {
      for(t in names(idps)) {
        fn1 = paste0('~/Desktop/tmp/ukb_idp/data/imagexcan_linear_adj.fourth_round.', idps[[t]], '.', models[[m]], '.csv')
        tmp1 = read.csv(fn1)
        df[[length(df) + 1]] = tmp1 %>% mutate(model = m, idp_type = t)
      }
    }
    df = do.call(rbind, df)
    # df$bhat = - df$bhat
    df_all = df
    df_all$model[df_all$model == 'EN'] = 'elastic net'
    df = df_all %>% filter(phenotype %in% pheno_interest) %>% mutate(idp_id = paste(idp_type, IDP, model)) %>% mutate(zscore = p2z(pval, bhat))
    saveRDS(df %>% select(-idp_id), fn)
    kk4 = df
  }
}

kk3 = readRDS('indiv_bxcan/dataframe_full.indiv_bxcan.rds')
tmp2 = inner_join(
  kk3 %>% select(IDP, model, phenotype, zscore), 
  kk4 %>% select(IDP, model, phenotype, zscore), 
  by = c('IDP', 'model', 'phenotype'), suffix = c('.resIDP', '.PCadj')
)
p = rbind(
  tmp1 %>% mutate(bxcan = 'S-BrainXcan'), 
  tmp2 %>% mutate(bxcan = 'BrainXcan')
) %>% ggplot() +
  geom_point(aes(x = zscore.resIDP, y = zscore.PCadj), alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') + 
  coord_equal() + 
  xlab('S-BrainXcan zscore of residual IDP') + 
  ylab('S-BrainXcan zscore of IDP \n when fitting with PC jointly') + facet_grid(bxcan~model) + 
  th2
ggsave(paste0(outdir, '/residp_vs_idp_pcadj.png'), p, width = 5, height = 5)