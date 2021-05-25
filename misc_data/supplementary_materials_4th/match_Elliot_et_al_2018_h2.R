# setwd('misc_data/supplementary_materials_4th/')
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
options(stringsAsFactors = F)

source('rlib.R')

outdir = 'match_Elliot_et_al_2018_h2'
dir.create(outdir)

map_file = 'Elliot_2018_to_UKB_field.tsv'
if(!file.exists(map_file)) {
  source('../supplementary_materials_3rd/rlib_match_Elliot_et_al_2018_h2.R')
  # h2 from Elliot et al 2018
  df_old = read.csv('../SuppTable_2_Elliot_et_al_2018.csv')
  # # our
  # df_new = rbind(
  #   read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/dmri.original.all_covar.w_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC'),
  #   read.table('~/Desktop/tmp/ukb_idp/heritability_2nd_round/t1.scaled.all_covar.w_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC')
  # )
  
  # map ids
  df_map = readRDS('../download_some_matching_files/cleanup_annot_our_idps.rds')
  # map another
  df_map2 = read.delim2('../supplementary_materials_3rd/supp_table_1.tsv', header = T)
  df_map = left_join(df_map2, df_map %>% select(FieldID, n2018_id), by = c('ukb_field' = 'FieldID'))
  
  # matched 
  tmp = df_map[ df_map$n2018_id %in% df_old$IDP, ]
  # not matched
  tmp_not = df_map[ !df_map$n2018_id %in% df_old$IDP, ]
  
  # FAST ROIs are not matched (by manually check)
  df_fast = df_old[ !is.na(stringr::str_match(df_old$IDP, '^IDP_T1_FAST_ROIs')), ]
  
  dd = mymatch(df_fast$IDP, tmp_not[, c('anatomy', 'left_or_right')])
  kk = df_fast[ is.na(dd$anatomy), ]
  
  dd[ is.na(dd$anatomy), ] =  mymatch(kk$IDP, tmp_not[, c('anatomy', 'left_or_right')], find_func = hard_find)
  df_fast = cbind(df_fast, dd)
  tmp_not$left_or_right = as.character(tmp_not$left_or_right)
  df_fast$lr = as.character(df_fast$lr)
  tmp_not = left_join(tmp_not, df_fast, c('left_or_right' = 'lr', 'anatomy'))
  
  
  # save Elliot et al IDP ID <-> UKB Field ID map
  df = rbind(
    tmp %>% select(ukb_field, n2018_id), 
    tmp_not %>% select(ukb_field, IDP) %>% rename(n2018_id = IDP)
  )
  write.table(df, map_file, sep = '\t', row = F, col = T, quo = F)
} else {
  df_map = read.table(map_file, header = T)
}

# h2 from Elliot et al 2018
df_old2 = read.csv('../SuppTable_2_Elliot_et_al_2018.csv')
# h2 from IDP GWAS 2021
df_old = read.csv('../idp_gwas_2021_supp_tab_1.csv') %>% mutate(IDP = paste0('IDP-', UKB.ID))
df_old = df_old %>% select(IDP, Heritability, Heritability.SE.)
colnames(df_old) = c('IDP', 'h2.2021', 'h2_SE.2021')
# our
df_new = rbind(
  read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_w_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC') %>% mutate(idp_type = 'dMRI'),
  read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_w_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC') %>% mutate(idp_type = 'T1')
)
# our no PC
df_new2 = rbind(
  read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.dmri_no_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC'),
  read.table('~/Desktop/tmp/ukb_idp/heritability_4th_round/fourth_round.t1_no_pc.tsv.gz', header = T) %>% filter(substr(phenotype, 1, 2) != 'PC')
)

df_old2 = inner_join(df_old2, df_map, by = c('IDP' = 'n2018_id'))
df_merge = inner_join(
  df_new %>% select(phenotype, idp_type, h2, h2_SE) %>% rename(IDP = phenotype),
  df_old2 %>% mutate(IDP = paste0('IDP-', ukb_field)) %>% select(IDP, h2, se..h2.) %>% rename(h2_SE = se..h2.),
  by = 'IDP',
  suffix = c('.Ours', '.2018')
)
df_merge = inner_join(
  df_merge,
  df_old,
  by = 'IDP'
)
df_merge = inner_join(
  df_new2 %>% select(phenotype, h2, h2_SE) %>% rename(IDP = phenotype, h2.Ours_no_PC = h2, h2_SE.Ours_no_PC = h2_SE),
  df_merge,
  by = 'IDP',
)

df_merge = remove_probtrack_idp(df_merge)

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1, col = 'red', ...)
}
png(paste0(outdir, '/h2_pairs.png'))
pairs(df_merge[, c('h2.Ours_no_PC', 'h2.Ours', 'h2.2018', 'h2.2021')], panel = my_line)
dev.off()
p1 = df_merge %>% ggplot() + 
  geom_errorbar(aes(x = h2.2021, ymax = h2.Ours_no_PC + 1.96 * h2_SE.Ours_no_PC, ymin = h2.Ours_no_PC - 1.96 * h2_SE.Ours_no_PC), alpha = 0.1) + 
  geom_errorbarh(aes(y = h2.Ours_no_PC, xmax = h2.2021 + 1.96 * h2_SE.2021, xmin = h2.2021 - 1.96 * h2_SE.2021), alpha = 0.1) +
  geom_point(aes(x = h2.2021, y = h2.Ours_no_PC)) + 
  geom_abline(slope = 1, intercept = 0) + th + 
  ylab('Heritability \n BrainXcan pipeline w/o PC adjustment') +
  xlab('Heritability \n Smith et al 2021')

p2 = df_merge %>% ggplot() + 
  geom_errorbar(aes(x = h2.2021, ymax = h2.Ours + 1.96 * h2_SE.Ours, ymin = h2.Ours - 1.96 * h2_SE.Ours), alpha = 0.1) + 
  geom_errorbarh(aes(y = h2.Ours, xmax = h2.2021 + 1.96 * h2_SE.2021, xmin = h2.2021 - 1.96 * h2_SE.2021), alpha = 0.1) +
  geom_point(aes(x = h2.2021, y = h2.Ours)) + 
  geom_abline(slope = 1, intercept = 0) + th + 
  ylab('Heritability \n BrainXcan pipeline with PC adjustment') +
  xlab('Heritability \n Smith et al 2021') 

ggsave(paste0(outdir, '/h2_smith_vs_bxcan_w_pc.png'), p2, height = 4, width = 4)
ggsave(paste0(outdir, '/h2_smith_vs_bxcan_no_pc.png'), p1, height = 4, width = 4)

p3 = df_merge %>% ggplot() + 
  geom_errorbar(aes(x = h2.Ours_no_PC, ymax = h2.Ours + 1.96 * h2_SE.Ours, ymin = h2.Ours - 1.96 * h2_SE.Ours), alpha = 0.1) + 
  geom_errorbarh(aes(y = h2.Ours, xmax = h2.Ours_no_PC + 1.96 * h2_SE.Ours_no_PC, xmin = h2.Ours_no_PC - 1.96 * h2_SE.Ours_no_PC), alpha = 0.1) +
  geom_point(aes(x = h2.Ours_no_PC, y = h2.Ours)) + 
  geom_abline(slope = 1, intercept = 0) + th + 
  ylab('Heritability \n BrainXcan pipeline with PC adjustment') +
  xlab('Heritability \n BrainXcan pipeline w/o PC adjustment') + coord_equal()

ggsave(paste0(outdir, '/h2_bxcan_w_pc_vs_no_pc.png'), p3, height = 4, width = 4)
