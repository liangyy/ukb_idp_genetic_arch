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
    make_option(c("-f", "--family_col"), type="character", default=NULL,
                help="Name of family ID column in phenotype table.",
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
                metavar="character"),
    make_option(c("-d", "--mode"), type="character", default='cv_performance',
                help="mode = cv_performance or model_training. Default = cv_performance.",
                metavar="character"),
    make_option(c("-x", "--remove_missing_in_phenotype"), action='store_true',
                help="If used, remove missing values for each of the phenotype")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

# helper functions
source('gw_lasso_helper.R')
library(dplyr)

# logging config
logging::basicConfig(level = opt$log_level)

# set random seed
set.seed(opt$random_seed)

# sanity check on mode
# we only allow: cv_performance and model_training for now!
if(opt$mode == 'cv_performance') {
  mode = 'cv_performance'
}  else if(opt$mode == 'model_training') {
  mode = 'model_training'
} else {
  message('Wrong mode = ', opt$mode)
  quit()
}
if(opt$remove_missing_in_phenotype) {
  message('--remove_missing_in_phenotype is effective')
}

# load snpnet configuration
snpnet_config = yaml::read_yaml(opt$snpnet_config)

logging::loginfo('Loading phenotype table.')
df_phenotype = load_phenotype(
  opt$phenotype_table, 
  indiv_col = opt$indiv_col,
  family_col = opt$family_col,
  pheno_list = opt$pheno_list
)
# df_phenotype: first 2 columns are FID and IID and they are followed by phenotypes.
nsample = nrow(df_phenotype)

if(mode == 'cv_performance') {
  logging::loginfo('Splitting into folds.')
  partitions = get_partitions(nsample, opt$nfold)
} else {
  # add trivial partition
  partitions = rep(2, nsample)
  # set nfold to 1 
  opt$nfold = 1
  # load snp meta info
  df_snp = data.table::fread(
    paste0('zstdcat ', opt$genotype, '.pvar.zst'), 
    header = T, 
    sep = '\t'
  )
  colnames(df_snp)[which(colnames(df_snp) == '#CHROM')] = 'CHR'
}


logging::loginfo('Looping over phenotypes.')
pred_list = list()
beta_list = list()
fid_iid = paste0(df_phenotype$FID, '_', df_phenotype$IID)
for(pheno in colnames(df_phenotype)[c(-1, -2)]) {
  myy = df_phenotype[[pheno]]
  myx = fid_iid
  mynsample = length(myy)
  if(opt$remove_missing_in_phenotype) {
    is.missing = is.na(myy)
    if(all(is.missing)) {
      message('All values are missing for phenotype ', pheno, '. Skip!')
      next
    }
    myy = myy[!is.missing]
    myx = fid_iid[!is.missing]
    if(mode == 'cv_performance') {
      partitions = get_partitions(mynsample, opt$nfold)
    } else {
      partitions = rep(2, mynsample)
    }
  }  
  df_out = data.frame(indiv = myx, yobs = myy, ypred = NA)
  for(k in 1 : opt$nfold) {
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Initialization.'))
    
    # build inner training split col
    train_idx = which(partitions != k)
    test_idx = which(partitions == k)
    ntrain = length(train_idx)
    inner_partitions = get_partitions(ntrain, opt$inner_nfold)
    valid_idx = train_idx[inner_partitions == 1]
    split_col = rep(NA, mynsample)
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
    cache_file = paste0(opt$output_prefix, '_', pheno, '_f', k, '.rds')
    write.table(df_phenotype, tmp_pheno_file, col = T, row = F, quo = F, sep = '\t')
    
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Inner fit (early stopping applied).'))
    snpnet_config[['early.stopping']] = TRUE
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
   
    # need some different handling when running on full data.
    to_predict = 'val'
    to_predict_idx = test_idx
    to_split_col = 'split_refit'
    if(mode == 'model_training') {
      to_predict = 'train'
      to_predict_idx = c(train_idx, valid_idx)
      to_split_col = NULL
    }    
 
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Full fit.'))
    snpnet_config[['early.stopping']] = FALSE
    full_fit = snpnet::snpnet(
      genotype.pfile = opt$genotype, 
      phenotype.file = tmp_pheno_file, 
      phenotype = pheno, 
      split.col = to_split_col, 
      configs = snpnet_config,
      lambda = inner_fit$full.lams[1:max_idx],
      alpha = opt$alpha,
      mem = opt$mem
    )
    
    logging::loginfo(paste0('Working on ', pheno, ': ', k, ' / ', opt$nfold, ' fold. Predict.'))
    full_pred = tryCatch(
      {
        snpnet::predict_snpnet(
          fit = full_fit,			       
          new_genotype_file = opt$genotype, 
          new_phenotype_file = tmp_pheno_file, 
          phenotype = pheno, 
          split_col = to_split_col, 
          split_name = to_predict,
          configs = snpnet_config
        )
      }, error = function(e) {
        list(
          prediction = list(
            val = matrix(mean(df_out$yobs[to_predict_idx]), ncol = length(inner_fit$full.lams[1:max_idx]), nrow = length(to_predict_idx))
          )
        )
      }
    )
    # save some intermediate results
    saveRDS(list(inner_fit = inner_fit, full_fit = full_fit, full_pred = full_pred), cache_file)
     
    # record results
    # prediction on held-out data (cv_performance) or full data (model_training)
    test_pred = full_pred$prediction[[to_predict]]
    opt_idx = min(which.max(inner_fit$metric.val), ncol(test_pred)) 
    ypred_test = test_pred[, opt_idx]  # this is from the best lambda
    if(!is.null(names(ypred_test))) {
      df_out$ypred[match(names(ypred_test), df_out$indiv)] = as.numeric(ypred_test)
    } else {
      # degenerate model
      df_out$ypred[ to_predict_idx ] = as.numeric(ypred_test)[1]
    }
    # save model weights
    if(mode == 'model_training') {
      mod = full_fit$beta[[opt_idx]]
      snpid_infos = parse_snp(names(mod))
      beta_out = data.frame(snpid = snpid_infos$snpid, alt = snpid_infos$allele, weight = as.numeric(mod))
      print(head(beta_out))
      beta_out = inner_join(beta_out, df_snp %>% select(ID, REF, ALT, CHR), by = c('snpid' = 'ID'))
      # the sign of the weights should be relative to df_snp as always
      beta_out$weight[beta_out$alt == beta_out$REF] = - beta_out$weight[beta_out$alt == beta_out$REF]
      beta_out = beta_out %>% select(-alt)
      beta_out = beta_out[ beta_out$weight != 0, ]
      beta_list[[length(beta_list) + 1]] = beta_out %>% mutate(phenotype = pheno)
    }
    
 
    # clean up intermediate file 
    system(paste0('rm ', tmp_pheno_file))
  }
  pred_list[[pheno]] = df_out
}

if(mode == 'cv_performance') {
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
} else {
  logging::loginfo('Saving model weights.')
  beta_list = do.call(rbind, beta_list)
  gz1 <- gzfile(paste0(opt$output_prefix, '.weights.tsv.gz'), "w")
  write.table(beta_list, gz1, quo = F, row = F, col = T, sep = '\t')
  close(gz1)
  logging::loginfo('Saving predicted values (in sample prediction).')
  saveRDS(pred_list, paste0(opt$output_prefix, '.y_insample.rds'))
}


logging::loginfo('Done.')

