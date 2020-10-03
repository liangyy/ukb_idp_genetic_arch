# This is a wrapper script for running snpnet on UKB IDPs.
# Input:
# 1. Genotype in plink2 PGEN format.
# 2. A phenotype table
# Output:
# Predicted phenotype via K-fold validation.
# Procedure:
# 1. Split the data into K-fold
# 2. For each fold, take the rest K-1 folds and do a inner round of validation.
#   2.1 In the inner round of validation, split the data into two parts, 
#       training and validation sets. 
#   2.2 Train a sequence of model along the regularization path using training
#       data and compute validation performance. To reduce the runtime, apply 
#       early stopping as implemented in snpnet.
#   2.3 Train a full model using the whole K-1 folds with 
#       the sequence of lambdas from the max to the best lambda values.
#   2.4 Make prediction on the held out fold.
# 3. Repeat step 2 for each of the fold to obtain Ypred for each fold.
# 4. Calculate R2, Pearson's correlation, Spearman's correlation 
#    between Ypred and Yobs.

library(optparse)

option_list <- list(
    make_option(c("-g", "--genotype"), type="character", default=NULL,
                help="Genotype prefix (PGEN format).",
                metavar="character"),
    make_option(c("-p", "--phenotype_table"), type="character", default=NULL,
                help="The phenotype table TSV with header.",
                metavar="character"),
    make_option(c("-n", "--nfold"), type="integer", default=5,
                help="The number of fold to use.",
                metavar="character"),
    make_option(c("-n2", "--inner_nfold"), type="integer", default=5,
                help="The number of fold to use for the inner loop of training.",
                metavar="character"),
    make_option(c("-i", "--indiv_col"), type="character", default=NULL,
                help="Name of individual ID column in phenotype table.",
                metavar="character"),
    make_option(c("-e", "--pheno_list"), type="character", default=NULL,
                help="The list of phenotype to work with.",
                metavar="character"),
    make_option(c("-s", "--snpnet_config"), type="character", default=NULL,
                help="A YAML file specifying the configuration of snpnet.",
                metavar="character"),
    make_option(c("-r", "--random_seed"), type="integer", default=1,
                help="Random seed to make the partition.",
                metavar="character"),
    make_option(c("-o", "--output_prefix"), type="character", default=NULL,
                help="Output prefix of: 1. a table of performance; 2. Yobs (for book-keeping)",
                metavar="character"),
    make_option(c("-l", "--log_level"), type="integer", default=9,
                help="Logging level",
                metavar="character"),
    make_option(c("-a", "--alpha"), type="numeric", default=1,
                help="The alpha value in glmnet",
                metavar="character"),
    make_option(c("-m", "--mem"), type="numeric", default=10000,
                help="Memory usage for plink2 in MB",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# helper functions
source('gw_lasso_helper.R')

# logging config
logging::basicConfig(level = opt$log_level)

# set random seed
set.seed(opt$random_seed)

# load snpnet configuration
snpnet_config = yaml::read_yaml(opt$snpnet_config)

logging::loginfo('Loading phenotype table.')
df_phenotype = load_phenotype(
  opt$phenotype_table, 
  indiv_col = opt$indiv_col,
  pheno_list = opt$pheno_list
)
# df_phenotype: first 2 columns are FID and IID and they are followed by phenotypes.
nsample = nrow(df_phenotype)

logging::loginfo('Splitting into folds.')
partitions = get_partitions(nsample, opt$nfold)

logging::loginfo('Looping over phenotypes.')
pred_list = list()
for(pheno in colnames(df_phenotype)[c(-1, -2)]) {
  df_out = data.frame(yobs = df_phenotype[[pheno]], ypred = NA)
  for(k in 1 : opt$nfold) {
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Initialization.'))
    
    # build inner training split col
    train_idx = which(partitions != k)
    test_idx = which(partitions == k)
    ntrain = length(train_idx)
    inner_partitions = get_partitions(ntrain, opt$inner_nfold)
    valid_idx = train_idx[inner_partitions == 1]
    split_col = rep(NA, nsample)
    split_col[train_idx] = 'train'
    split_col[valid_idx] = 'val'
    split_col[test_idx] = 'test'
    if(sum(is.na(split_col))) {
      stop('The split_col construction is failed.')
    }
    df_phenotype$split = split_col
    
    # build refit col
    split_col[split_col == 'val'] = 'train'
    split_col[split_col == 'test'] = 'val'
    df_phenotype$split_refit = split_col
    
    # write the temporary phenotype file
    tmp_pheno_file = paste0(opt$output_prefix, '_', pheno, '_f', k, '.phe')
    write.table(df_phenotype, tmp_pheno_file, col = T, row = F, quo = F, sep = '\t')
    
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Inner fit (early stopping applied).'))
    snpnet_config$early.stopping = TRUE
    inner_fit = snpnet::snpnet(
      genotype.pfile = opt$genotype, 
      phenotype.file = tmp_pheno_file, 
      phenotype = pheno, 
      split.col = "split", 
      configs = snpnet_config,
      alpha = opt$alpha,
      mem = opt$mem
    )
    max_idx <- sum(!is.na(inner_fit$metric.val))
    
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Full fit.'))
    snpnet_config$early.stopping = FALSE
    full_fit = snpnet::snpnet(
      genotype.pfile = opt$genotype, 
      phenotype.file = tmp_pheno_file, 
      phenotype = pheno, 
      split.col = "split_refit", 
      configs = snpnet_config,
      lambda = inner_fit$full.lams[1:max_idx],
      alpha = opt$alpha,
      mem = opt$mem
    )
    
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Predict.'))
    full_pred = snpnet::predict_snpnet(
      fit = full_fit,			       
      new_genotype_file = opt$genotype, 
      new_phenotype_file = tmp_pheno_file, 
      phenotype = pheno, 
      split_col = "split_refit", 
      split_name = 'val',
      configs = snpnet_config
    )
    test_pred = full_pred$prediction$val
    ypred_test = test_pred[, ncol(test_pred)]  # this is from the best lambda
    df_out$ypred[test_idx] = ypred_test
    
    # clean up intermediate file 
    system(paste0('rm ', tmp_pheno_file))
  }
  pred_list[[pheno]] = df_out
}

logging::loginfo('Calculate performance.')
out = list()
for(pheno in names(pred_list)) {
  df_perf = eval_perf(pred_list[[pheno]]$ypred, pred_list[[pheno]]$yobs)
  df_perf$phenotype = pheno
  out[[length(out) + 1]] = df_perf
}
out = do.call(rbind, out)

logging::loginfo('Writing results to disk.')
saveRDS(pred_list, paste0(opt$output_prefix, '.yval.rds'))
write.table(
  out, paste0(opt$output_prefix, '.performance.tsv'), 
  col = T, row = F, quo = F, sep = '\t'
)

logging::loginfo('Done.')

