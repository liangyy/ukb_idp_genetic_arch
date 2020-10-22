# this script is motivated by `rmd/export_mat_fac.R`
# we regress out some IDP PCs from IDP residuals (Owen prepared) so that the 
# IDPs are less correlated.

library(dplyr)

df = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_final-phenotypes.parquet')

# IDP residual matrix
mat = as.matrix(df[, -1])

# center and standardize IDP matrix
mat_centered = apply(mat, 2, function(x){(x - mean(x)) / sd(x)})
# calculate X' X / N for PCA
cor_mat2 = t(mat_centered) %*% mat_centered / nrow(mat_centered)
# EVD
eig_res = eigen(cor_mat2)

# loop over a list of # PCs
for(top_n in c(0, 1, 2, 5, 10, 20, 50)) {
  if(top_n == 0) {
    pve = 0
    mat_res2 = mat_centered
  } else {
    pve = sum(eig_res$values[1:top_n]) / sum(eig_res$values)
    message('Working on nPC = ', top_n, ' and PVE = ', signif(pve, digits = 3))
    
    top_n_pc = mat_centered %*% eig_res$vectors[, 1 : top_n]
    tmp = qr(top_n_pc)
    Q_ = qr.Q(tmp)
    
    # mat_res = mat - Q_ %*% (t(Q_) %*% mat)
    mat_res2 = mat_centered - Q_ %*% (t(Q_) %*% mat_centered)
  }
  
  # standardize
  mat_res2 = apply(mat_res2, 2, function(x){(x - mean(x)) / sd(x)})
  
  mat_out = as.data.frame(mat_res2)
  mat_out$individual = df$individual
  mat_out = mat_out[, c(ncol(mat_out), 1 : (ncol(mat_out) - 1))]
  arrow::write_parquet(mat_out, paste0('/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.regress_out_', top_n, 'PCs.parquet'))
  
  # paste0('~/Desktop/tmp/ukb_idp/regress_out_idp_pcs/2020-05-18_final-phenotypes.regress_out_', top_n, 'PCs.parquet')
}
