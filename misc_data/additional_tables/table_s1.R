# table S1
library(dplyr)
outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- read.delim2('/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_1_renamed.tsv')
dmris <- c('FA', 'ICVF', 'ISOVF', 'OD')
grays <- c('Cerebellum', 'Cortical')
sbc <- 'Subcortical'
map <- data.frame(
  idp = c(
    dmris,
    paste0('w-', dmris),
    paste0('Gray-', grays),
    paste0(c('Gray-', ''), sbc)
  ),
  pc1_name = c(
    paste0('PC-', dmris, '-TBSS-1'),
    paste0('PC-', dmris, '-ProbTrack-1'),
    paste0('PC-', grays, '-1'),
    paste0('PC-', sbc, c('_GMvol', '_vol'), '-1')
  )
)
df <- df %>% 
  # filter(substr(subtype, 1, 1) != 'w') %>%
  mutate(IDP = paste0('IDP-', ukb_field)) %>%
  left_join(map, by = c('subtype' = 'idp')) %>%
  select(IDP, modality, subtype, pc1_name, region, side, measurement_type,
         dmri_measure, t1_anatomy_group, notes, ukb_link) 
  
write.csv(df, paste0(outdir, '/table_s1.csv'), row.names = FALSE)
