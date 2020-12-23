library(optparse)

option_list <- list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Mode of the scaling: scaled/regress/original.",
              metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="Type of the data: t1/dmri.",
              metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file name.",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# working directory: misc_data/second_round_idp_preprocessing
source('scripts/rlib.R')

# opt = list(type = 't1', mode = 'original')
df_idp = data.table::fread('~/Desktop/tmp/ukb_idp/idp_phenotypes/archive/2020-05-07_idp-phenotypes_cleaned.txt', sep = '\t', data.table = F)

message(paste0('type = ', opt$type))
if(opt$type == 't1') {
  annot = readRDS('../process_t1/t1_meta.rds')
  normalized = !is.na(stringr::str_match(as.character(annot$Field), '\\(normalised for head size\\)$'))
  annot = annot[ !normalized, ]
} else {
  annot = readRDS('../download_some_matching_files/annot_dmri_idps.rds')
}

df_selected = df_idp[, c('individual', paste0('IDP-', annot$FieldID))]

size_factor = df_idp[, 'IDP-25000']

message(paste0('mode = ', opt$mode))
if(opt$mode == 'scaled') {
  df_selected[, -1] = df_selected[, -1] * size_factor
  # df_selected[, -1] = apply(df_selected[, -1], 2, inv_norm)
} else if(opt$mode == 'regress') {
  # df_selected[, -1] = apply(df_selected[, -1], 2, inv_norm)
  # size_factor_inrt = inv_norm(size_factor)
  df_selected[, -1] = apply(df_selected[, -1], 2, function(y){ regress_out_univariate(y = y, x = size_factor) })
  # df_selected[, -1] = apply(df_selected[, -1], 2, inv_norm)
} else if(opt$mode == 'original') {
  # df_selected[, -1] = apply(df_selected[, -1], 2, inv_norm)
} else {
  message(paste0('Wrong mode = ', opt$mode))
  quit()
}
# outliers = rep(F, nrow(df_selected))
# for(i in 2 : ncol(df_selected)) {
#   sd_tmp = sd(df_selected[, i])
#   mean_tmp = mean(df_selected[, i])
#   out_of_range = (df_selected[, i] > mean_tmp + sd_tmp * 4) | (df_selected[, i] < mean_tmp - sd_tmp * 4)
#   outliers[out_of_range] = T 
# }
# df_selected = df_selected[!outliers, ]
df_selected$individual = as.character(df_selected$individual)
arrow::write_parquet(df_selected, opt$output)

