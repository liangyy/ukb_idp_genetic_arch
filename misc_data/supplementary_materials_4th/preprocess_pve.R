# setwd('misc_data/supplementary_materials_4th/')
library(dplyr)
library(ggplot2)
source('rlib.R')
theme_set(theme_classic(base_size = 12))

# df_pve = data.frame(
#   modality = c('FA', 'ICVF', 'ISOVF', 'OD', 'Gray-Cortical', 'Gray-Subcortical', 'Subcortical', 'Gray-Cerebellum', 'w-FA', 'w-ICVF', 'w-ISOVF', 'w-OD')
# )

dmris = c('FA', 'ICVF', 'ISOVF', 'OD')
df = list()
for(dd in dmris) {
  temp = readRDS(paste0('../fourth_round_idp_preprocessing/output/', dd, '.pca_results.rds'))
  df[[length(df) + 1]] = data.frame(
    modality = c(dd, paste0('w-', dd)), 
    pve = c(temp$TBSS$pve, temp$ProbTrack$pve)
  )
}
t1s = c('Cortical', 'Subcortical_vol', 'Subcortical_GMvol', 'Cerebellum')
for(tt in t1s) {
  temp = readRDS(paste0('../fourth_round_idp_preprocessing/output/', tt, '.pca_results.rds'))
  df[[length(df) + 1]] = data.frame(
    modality = pc_map(tt), 
    pve = temp$pve
  )
}
df = do.call(rbind, df)
p = df %>% ggplot() + geom_bar(aes(x = modality, y = pve), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave('preprocess_pve.png', p, height = 3, width = 6)


# count the number of IDPs per subtype
dfc = list()
dmri_subtype = yaml::read_yaml('../fourth_round_idp_preprocessing/output/dmri_covar.yaml')
for(n in names(dmri_subtype)) {
  dfc[[length(dfc) + 1]] = data.frame(subtype = n, nidp = length(dmri_subtype[[n]]$x))
}
t1_subtype = yaml::read_yaml('../fourth_round_idp_preprocessing/output/t1_covar.yaml')
for(n in names(t1_subtype)) {
  dfc[[length(dfc) + 1]] = data.frame(subtype = n, nidp = length(t1_subtype[[n]]$x))
}
dfc = do.call(rbind, dfc)
summary(dfc)
