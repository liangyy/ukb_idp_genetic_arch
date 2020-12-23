library(optparse)

option_list <- list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Mode of the PCA adjustment: w_pc/no_pc",
              metavar="character"),
  make_option(c("-i", "--input_matrix"), type="character", default=NULL,
              help="Input IDP matrix.",
              metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL,
              help="Prefix of the output files.",
              metavar="character"),
  make_option(c("-p", "--pve_cutoff"), type="numeric", default=NULL,
              help="PVE cutoff for the adjustment.",
              metavar="character"),
  make_option(c("-a", "--annot_type"), type="character", default=NULL,
              help="Type of annotation: t1/dmri",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# working directory: misc_data/second_round_idp_preprocessing
source('scripts/rlib.R')
library(ggplot2)
library(dplyr)

# opt = list(input_matrix = 'tmp_original.parquet', pve_cutoff = 0.2, annot_type = 't1')

df_data = arrow::read_parquet(opt$input_matrix)
# x = as.matrix(df_data[, -1])

message(paste0('type = ', opt$annot_type))
if(opt$annot_type == 't1') {
  annot = readRDS('../process_t1/t1_meta.rds')
  normalized = !is.na(stringr::str_match(as.character(annot$Field), '\\(normalised for head size\\)$'))
  annot = annot[ !normalized, ]
} else {
  annot = readRDS('../download_some_matching_files/annot_dmri_idps.rds')
}

message(paste0('mode = ', opt$mode))
if(opt$mode == 'w_pc') {
  res = split_by_pca(as.matrix(df_data[, -1]), pve_cutoff = opt$pve_cutoff, skip = F)
  mat_res = res$residual
  pc_all = res$pc
  mat_res = as.data.frame(mat_res)
  colnames(mat_res) = colnames(df_data)[-1]
  pc_all = as.data.frame(pc_all)
  colnames(pc_all) = paste0('PC-', 1 : ncol(pc_all))
  df_out = data.frame(individual = df_data$individual)
  df_out = cbind(df_out, mat_res, pc_all)
  # save the pca results
  saveRDS(list(pc_loadings = res$pc_loadings, pve = res$pve), paste0(opt$output_prefix, '.pca_results.rds'))
} else if(opt$mode == 'no_pc') {
  res = split_by_pca(as.matrix(df_data[, -1]), pve_cutoff = opt$pve_cutoff, skip = T)
  mat_res = res$residual
  mat_res = as.data.frame(mat_res)
  colnames(mat_res) = colnames(df_data)[-1]
  df_out = data.frame(individual = df_data$individual)
  df_out = cbind(df_out, mat_res)
} else {
  message(paste0('Wrong mode = ', opt$mode))
  quit()
}

corr = cor(mat_res)
df_corr = corr %>% reshape2::melt() %>% 
  inner_join(annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% 
               select(idp), by = c('Var1' = 'idp')) %>%
  inner_join(annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% 
               select(idp), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
p0 = myplot_magic(df_corr)
# p0
ggsave(paste0(opt$output_prefix, '_residual_corr.png'), p0)
# plot the first two IDP after the PC adjustment as sanity check
p1 = data.frame(idp1 = mat_res[, 1], idp2 = mat_res[, 2]) %>% ggplot() + geom_point(aes(x = idp1, y = idp2))
ggsave(paste0(opt$output_prefix, '_first_two_idps.png'), p1)

# do inverse normalization before saving
message('Performing inverse normalization on residuals.')
df_out[, -1] = apply(df_out[, -1], 2, inv_norm)
arrow::write_parquet(df_out, paste0(opt$output_prefix, '.parquet'))


