load_phenotype = function(phenotype_file, indiv_col, pheno_list) {
  df_pheno = data.table::fread(phenotype_file, sep = '\t', header = T, data.table = F)
  # replace '-' with 'x' in the string
  colnames(df_pheno) = fix_str(colnames(df_pheno))
  phenotypes = load_list(pheno_list)
  df_indiv = data.frame(FID = df_pheno[[indiv_col]], IID = df_pheno[[indiv_col]])
  df_pheno = df_pheno[, phenotypes, drop = F]
  cbind(df_indiv, df_pheno)
}

fix_str = function(str_) {
  stringr::str_replace(str_, '-', 'x')
}

load_list = function(fn) {
  fix_str(read.table(fn, header = F, stringsAsFactors = F)$V1)
}

get_partitions = function(ntotal, nfold) {
  size = floor(ntotal / nfold)
  rest = ntotal - nfold * size
  out = c()
  for(i in 1 : nfold) {
    tmp = rep(i, size)
    if(i <= rest) {
      tmp = c(tmp, i)
    }
    out = c(out, tmp)
  }
  sample(out, length(out), replace = FALSE)
}

eval_perf = function(ypred, yobs) {
  # return R2, Pearson Cor, Spearman Cor
  total_error = mean( (yobs - mean(yobs)) ^ 2 )
  total_residual_error = mean( (yobs - ypred) ^ 2 )
  R2 = 1 - total_residual_error / total_error
  pearson = cor(ypred, yobs, method = 'pearson')
  spearman = cor(ypred, yobs, method = 'spearman')
  data.frame(R2 = R2, Pearson = pearson, Spearman = spearman)
}
