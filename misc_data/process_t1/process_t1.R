# setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/process_t1/')

library(dplyr)
library(ggplot2)

# For the T1
# We perform PCA on the T1 matrix.

# For each PCA analysis, we keep both the aggregated PC and the residual.
# The number of PCs to regress out depends on the overall PVE. 
# We aim at PVE >= 30%

pve_cutoff = 0.3

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
  # mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% mutate(order = order(Var1)) %>% group_by(matter_type.1) %>% summarize(pos = mean(order))
  mytmp %>% ggplot() +
    theme_bw() + 
    geom_raster(aes(x = as.character(Var1), y = as.character(Var2), fill = value)) + scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    # annotate(geom = "text", x = mytmp_g$pos, y = -20, label = mytmp_g$matter_type.1, size = 4) + 
    coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
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
# the RDS in the below is generated in ../download_some_matching_files/explore_dmri.R
df_annot = readRDS('../process_t1/t1_meta.rds')
df = df[, c('individual', paste0('IDP-', as.character(df_annot$FieldID)))]
mat = as.matrix(df[, -1])

# PCA starts here
res = split_by_pca(mat, pve_cutoff = pve_cutoff)
mat_res = res$residual
pc_all = res$pc

corr = cor(mat_res)
df_corr = corr %>% reshape2::melt() %>% 
  inner_join(df_annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% 
  select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var1' = 'idp')) %>%
  inner_join(df_annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% 
  select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
p0 = myplot_magic(df_corr)
ggsave('t1_residual_corr.png', p1)


mat_res = as.data.frame(mat_res)
colnames(mat_res) = colnames(df)[-1]
pc_all = as.data.frame(pc_all)
colnames(pc_all) = paste0('PC-', 1 : ncol(pc_all))
df_out = data.frame(individual = df$individual)
df_out = cbind(df_out, mat_res, pc_all)
tmp = r_to_py(df_out)
tmp$to_parquet('/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_T1.parquet', index = F)
