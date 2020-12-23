library(optparse)

option_list <- list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Mode of the regressing out: non_idp_covar/all_covar.",
              metavar="character"),
  make_option(c("-i", "--input_matrix"), type="character", default=NULL,
              help="Input IDP matrix.",
              metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file name.",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# working directory: misc_data/second_round_idp_preprocessing
source('scripts/rlib.R')

# opt = list(input_matrix = 'tmp.parquet', mode = 'all_covar')

df_data = arrow::read_parquet(opt$input_matrix, as_data_frame = T)
df_covar = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_covariates.parquet', as_data_frame = T)
df_covar$individual = as.character(df_covar$individual)
df_idp = data.table::fread('~/Desktop/tmp/ukb_idp/idp_phenotypes/archive/2020-05-07_idp-phenotypes_cleaned.txt', sep = '\t', data.table = F)
df_idp$individual = as.character(df_idp$individual)
non_idp_covar = c(paste0('pc', 1 : 10), 'age_recruitment', 'sex', 'age_recruitment_2', 'age_sex', 'age_2_sex')
idp_covar = paste0('IDP-', 25756 : 25759)
df_non_idp = df_covar[, c('individual', non_idp_covar)]
df_idp = df_idp[, c('individual', idp_covar)]

df_non_idp = df_non_idp[rowSums(is.na(df_non_idp[, -1])) == 0, ]
df_idp = df_idp[rowSums(is.na(df_idp[, -1])) == 0, ]
df_data = df_data[rowSums(is.na(df_data[, -1])) == 0, ]

indiv_common = intersect(intersect(df_non_idp$individual, df_idp$individual), df_data$individual)
df_data = df_data[match(indiv_common, df_data$individual), ]
df_idp = df_idp[match(indiv_common, df_idp$individual), ]
df_non_idp = df_non_idp[match(indiv_common, df_non_idp$individual), ]

message(paste0('mode = ', opt$mode))
if(opt$mode == 'all_covar') {
  covar_mat = as.matrix(cbind(df_idp[, -1], df_non_idp[, -1]))
} else if(opt$mode == 'non_idp_covar') {
  covar_mat = as.matrix(df_non_idp[, -1])
} else {
  message(paste0('Wrong mode = ', opt$mode))
  quit()
}

df_data[, -1] = regress_out_matrix(as.matrix(df_data[, -1]), as.matrix(covar_mat))

# message('Performing inverse normalization on residuals.')
# df_data[, -1] = apply(df_data[, -1], 2, inv_norm)
arrow::write_parquet(df_data, opt$output)
