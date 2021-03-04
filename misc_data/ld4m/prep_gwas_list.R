setwd('misc_data/ld4m/')
library(dplyr)

df_tab = read.csv('ld4m_table1.csv') %>% mutate(trait = tolower(Trait))
df_data = read.csv('ld4m_gwas_list.csv') %>% mutate(trait = tolower(TRAIT))

df_map = list(
  `Number of children` = 'Number Children (Pooled)',
  `College` = 'College Education',
  `FVC` = 'Forced Vital Capacity (FVC)',
  `CVD including HT` = 'Cardiovascular Diseases',
  `BP systolic` = 'Systolic Blood Pressure',
  `FEV1/FVC` = 'FEV1-FVC Ratio',
  `WBC count` = 'White Blood Cell Count',
  `AID` = 'Auto Immune Traits',
  `BMD heel` = 'Heel T Score',
  `Type II diabetes` = 'Type 2 Diabetes',
  `RBC distribution width` = 'Red Blood Cell Distribution Width',
  `RBC count` = 'Red Blood Cell Count',
  `Sunburn` = 'Sunburn Occasion'
)
df_map = data.frame(trait = tolower(names(df_map)), new_name = tolower(unlist(df_map)))
df_map$new_name = as.character(df_map$new_name)

df_tab = left_join(df_tab, df_map, by = 'trait')
df_tab$new_name[is.na(df_tab$new_name)] = df_tab$trait[is.na(df_tab$new_name)]

df_sub = df_tab %>% inner_join(df_data, by = c('new_name' = 'trait')) 
write.table(df_sub$FILE, 'ld4m_external_traits.txt', col = F, row = F, quo = F)

df_sub2 = df_tab %>% left_join(df_data, by = c('new_name' = 'trait')) 
df_map2 = list(
  `Age at first birth F` = 'PASS_AgeFirstBirth',
  `Schizophrenia` = 'PASS_Schizophrenia',
  `RA` = 'PASS_Rheumatoid_Arthritis'
)
df_map2 = data.frame(trait = tolower(names(df_map2)), FILE = unlist(df_map2))
kk = df_sub2$trait[df_sub2$trait %in% df_map2$trait]
df_map2 = df_map2[ match(kk, df_map2$trait), ]
df_sub2$FILE[df_sub2$trait %in% kk] = df_map2$FILE
saveRDS(df_sub2, 'ld4m_map.rds')