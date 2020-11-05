# setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/process_dmri/')

library(dplyr)
library(ggplot2)

# For the dMRI
# We perform PCA on the dMRI matrix.

# For each PCA analysis, we keep both the aggregated PC and the residual.
# The number of PCs to regress out depends on the overall PVE. 
# We aim at PVE >= 50%

pve_cutoff = 0.5

standardize = function(x) {
  apply(x, 2, function(y) { (y - mean(y)) / sd(y) })
}
split_by_pca = function(x, pve_cutoff = 0.5) {
  x = standardize(x)
  res = svd(x)
  pve = cumsum(res$d^2 / sum(res$d^2))
  npc = sum(pve <= pve_cutoff) + 1
  pc_mat = res$u[, 1 : npc, drop = F]
  res = x - pc_mat %*% (t(pc_mat) %*% x)
  list(residual = res, pc = pc_mat)
}


myplot_magic = function(df) {
  mytmp = df # 
  mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% mutate(order = order(Var1)) %>% group_by(measure.1) %>% summarize(pos = mean(order))
  mytmp %>% ggplot() +
    theme_bw() + 
    geom_raster(aes(x = as.character(Var1), y = as.character(Var2), fill = value)) + scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    annotate(geom = "text", x = mytmp_g$pos, y = -20, label = mytmp_g$measure.1, size = 4) + coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
    theme(
      plot.margin = unit(c(0.1, 0.1, 2, 0.1), "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}

# to work on washington
# we need to use reticulate to load parquet using pandas from python
system('export RETICULATE_PYTHON=/vol/bmd/yanyul/miniconda3/envs/ukb_idp/bin/python') 
library(reticulate)
pd = import('pandas')
df = pd$read_parquet('/vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet')

# df = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_final-phenotypes.parquet')
df_annot = readRDS('../download_some_matching_files/annot_dmri_idps.rds')
df = df[, c('individual', paste0('IDP-', as.character(df_annot$FieldID)))]
mat = as.matrix(df[, -1])

# PCA starts here
res = split_by_pca(mat, pve_cutoff = pve_cutoff)
mat_res = res$residual
pc_all = res$pc

corr = cor(mat_res)
df_corr = corr %>% reshape2::melt() %>% 
  inner_join(df_annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, type, position, measure, lr), by = c('Var1' = 'idp')) %>%
  inner_join(df_annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, type, position, measure, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
df_corr = df_corr %>% mutate(
  label1 = paste0(Var1, '\n', measure.1),
  label2 = paste0(Var2, '\n', measure.2)
)
p1 = myplot_magic(df_corr %>% filter(type.1 == 'mean', type.2 == 'mean') )
p2 = myplot_magic(df_corr %>% filter(type.1 == 'weighted_mean_in_tract', type.2 == 'weighted_mean_in_tract') )
p3 = myplot_magic(df_corr %>% filter(type.1 == 'weighted_mean_in_tract', type.2 == 'mean') )
ggsave('dmri_residual_mean_corr.png', p1)
ggsave('dmri_residual_weighted_mean_corr.png', p2)
ggsave('dmri_residual_cross_corr.png', p3)


mat_res = as.data.frame(mat_res)
colnames(mat_res) = colnames(df)[-1]
pc_all = as.data.frame(pc_all)
colnames(pc_all) = paste0('PC-', 1 : ncol(pc_all))
df_out = data.frame(individual = df$individual)
df_out = cbind(df_out, mat_res, pc_all)
tmp = r_to_py(df_out)
tmp$to_parquet('/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_dMRI.parquet', index = F)
