load_phenotype = function(phenotype_file, indiv_col, pheno_list, family_col = NULL) {
  df_pheno = data.table::fread(phenotype_file, sep = '\t', header = T, data.table = F)
  # replace '-' with 'x' in the string
  colnames(df_pheno) = fix_str(colnames(df_pheno))
  phenotypes = load_list(pheno_list)
  if(is.null(family_col)) {
    family_col = indiv_col
  }
  df_indiv = data.frame(FID = df_pheno[[family_col]], IID = df_pheno[[indiv_col]])
  df_pheno = df_pheno[, phenotypes, drop = F]
  # remove phenotypes with constant value
  pheno_sd = apply(df_pheno, 2, sd)
  df_pheno = df_pheno[, pheno_sd > 1e-10, drop = F]
  cbind(df_indiv, df_pheno)
}

fix_str = function(str_) {
  stringr::str_replace_all(str_, '-', 'x')
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

parse_snp = function(str) {
  str = as.character(str)
  tmp = strsplit(str, '_')
  snpid = unlist(lapply(tmp, function(x) { paste0(x[-length(x)], collapse = '_') }))
  allele = unlist(lapply(tmp, function(x) { x[length(x)] }))
  list(snpid = snpid, allele = allele)
}
