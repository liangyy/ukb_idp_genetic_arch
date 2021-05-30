library(dplyr)
source('rlib.R')
df = read.delim2('supp_table_1.tsv')
idp = load_idp_annot()
df = left_join(df %>% mutate(IDP = paste0('IDP-', ukb_field)), idp %>% select(IDP, subtype), by = 'IDP') %>% 
  select(-IDP)
sum(is.na(df$subtype))
write.table(
  df %>% rename(
    modality = t1_or_dmri,
    region = anatomy,
    side = left_or_right
  ), 
  'supp_table_1_renamed.tsv', quote = F, row.names = F, sep = '\t'
)

