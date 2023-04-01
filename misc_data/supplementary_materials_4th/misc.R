library(dplyr)
library(ggplot2)
# https://academic.oup.com/cercor/article/30/4/2307/5669893
# compare heritability estimates
mm <- readxl::read_xlsx('~/Documents/repo/overleaf/603ea4a4adb746d6e98dec93/tables/Table_S2.xlsx')
df0 <- readxl::read_xlsx('~/Documents/repo/overleaf/603ea4a4adb746d6e98dec93/tables/Table_S3.xlsx')
df0 <- df0 %>% left_join(mm, by = 'IDP')
df00 <- df0 %>% filter(subtype == 'Subcortical')
df1 <- readxl::read_xlsx('~/Downloads/tables1_hsq_alldatasets_maf001_bhz241.xlsx', sheet = 'hsq_allDatasets_maf001')
df1 <- df1 %>% filter(dataset == 'UKB-22110')
kk <- inner_join(df00, df1, by = c('region' = 'pheno')) %>%
  mutate(h2.ref = as.numeric(`V(G)/Vp`))
kk %>% ggplot() +
  geom_point(aes(x = h2, y = h2.ref))
