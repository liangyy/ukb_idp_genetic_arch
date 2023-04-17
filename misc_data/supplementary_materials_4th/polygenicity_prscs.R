# setwd('misc_data/supplementary_materials_4th/')

backup_dir = '~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/'

# setwd(backup_dir)
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

options(stringsAsFactors = F)

source('rlib.R')

foldern = '~/Documents/repo/GitHub/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/polygenicity_prscs'
dir.create(foldern)

color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')

# load h2
{ 
df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T)
df1$is_pc = substr(df1$phenotype, 1, 2) == 'PC'
df2$is_pc = substr(df2$phenotype, 1, 2) == 'PC'
df1$pc = rep('Region-Specific', nrow(df1))
df1$pc[df1$is_pc] = 'Common Factor'
df2$pc = rep('Region-Specific', nrow(df2))
df2$pc[df2$is_pc] = 'Common Factor'
df_h2 = rbind(
  df1 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), 
  df2 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))
)
df_h2 = df_h2 %>% rename(IDP = phenotype)
}

pcs = df_h2 %>% filter(is_pc) %>% select(IDP, idp_type) 

# load polygenicity & h2
{
  df_poly_h2 = read.csv('~/Downloads/Table_S3.xlsx - ..csv')
  df_poly_h2 = df_poly_h2 %>% mutate(stable_Me = Me > 0 & Me > 1.96 * Me_SE)
  df_poly_h2_0 <- df_poly_h2
  # remove TBSS IDPs
  df_poly_h2 = remove_probtrack_idp(df_poly_h2)
}

# load prediction performance
{
  df_perf = read.csv('~/Downloads/Table_S4.xlsx - ..csv')
  df_perf_prscs = rbind(
    read.table('~/Downloads/prscs_prediction_perf-dmri.txt', header = TRUE),
    read.table('~/Downloads/prscs_prediction_perf-t1.txt', header = TRUE)) %>%
    mutate(model_name = 'PRS-CS') %>%
    rename(Spearman = Spearman_prscs)
  df_perf = rbind(
    df_perf %>% select(IDP, Spearman, model_name),
    df_perf_prscs
  )
  model_name_rename = data.frame(
    model_name = c('elastic net', 'ridge', 'PRS-CS'),
    new_name = c('en', 'ridge', 'prscs')
  )
  df_perf = df_perf %>% left_join(model_name_rename, by = 'model_name')
  df_perf = df_perf %>% reshape2::dcast(IDP ~ new_name, value.var = 'Spearman')
  df_perf = df_perf %>% 
    mutate(
      Spearman.ridge_vs_prscs = ridge - prscs,
      Spearman.en_vs_prscs = en - prscs,
      Spearman.ridge_vs_en = ridge - en
    )
}

plot_h2_vs_me = T
if(isTRUE(plot_h2_vs_me)) {
  tmp = inner_join(df_poly_h2, df_perf, by = 'IDP')
  is_pc = substr(tmp$IDP, 1, 2) == 'PC'
  tmp$pc = rep('Region-Specific', nrow(tmp))
  tmp$pc[is_pc] = 'Common Factor'
  tmp %>% filter(stable_Me) %>% 
    ggplot() + geom_point(aes(x = h2, y = Me))
  pp = tmp %>% filter(stable_Me) %>% 
    ggplot() + geom_point(aes(x = Me, y = Spearman.ridge_vs_prscs, color = pc)) + 
    scale_x_log10() +
    xlab('Estimated polygenicity (Me)') +
    ylab('Difference in Spearman correlation \n between ridge and PRS-CS') +
    th +
    theme(legend.position = c(0.8, 0.18), legend.title = element_blank()) +
    scale_color_manual(values = c('Common Factor' = 'red', 'Region-Specific' = 'black'))
  ggsave(paste0(foldern, '/pred_perf_vs_Me.png'), pp, width = 5.5, height = 5)
}

# # count number
# {
#   df_poly = df_poly %>% mutate(signif = Ma_est - 1.96 * Ma_err > 0 | Ma_est + 1.96 * Ma_err < 0, pos_signif = Ma_est - 1.96 * Ma_err > 0)
#   df_poly_signif = df_poly %>% filter(pos_signif == T)
#   message('Significant under alpha = 0.05: ', sum(df_poly$signif), ' out of ', nrow(df_poly))
#   message('Significant under alpha = 0.05 (positive only): ', sum(df_poly$pos_signif), ' out of ', nrow(df_poly))
#   df_poly_signif %>% summarize(median(Ma_est), min(Ma_est), max(Ma_est))
# }
