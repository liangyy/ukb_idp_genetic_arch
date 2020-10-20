library(dplyr)
library(ggplot2)
library(patchwork)
library("gridExtra")
theme_set(theme_bw(base_size = 15))

source('https://gist.githubusercontent.com/liangyy/489d1519dd45246caf4756d7722bfa25/raw/90c572c0a287f2cada53811c7cd51ec14fade488/fast_linear_regression')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('https://gist.githubusercontent.com/liangyy/4c647634fe00b3f042ebd1599dda65c7/raw/9977562b65d0fb63a693fa7fa60035a37641ad2f/qqplot_by_group')


set.seed(1)
df = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_final-phenotypes.parquet')
mat = as.matrix(df[, -1])

plist = c()
clist = c()
for(n in 1:10) {
  select_idx = sample(ncol(mat), 1); y = rnorm(nrow(mat)) * 20 + rowSums(mat[, select_idx, drop = F])
  message(select_idx)
  res = fast_linear_regression(as.numeric(y), mat, covariate = matrix(1, nrow = nrow(mat), ncol = 1))
  plist = c(plist, res$pval)
  clist = c(clist, rep(paste0('repeat', n), ncol(mat)))
}

p = qqplot_by_group(plist, clist) + th + theme(legend.position = 'none'); p


plist = c()
clist = c()
for(n in 1:10) {
  y = rnorm(nrow(mat)) * 20 
  res = fast_linear_regression(as.numeric(y), mat, covariate = matrix(1, nrow = nrow(mat), ncol = 1))
  plist = c(plist, res$pval)
  clist = c(clist, rep(paste0('repeat', n), ncol(mat)))
}

p2 = qqplot_by_group(plist, clist) + th + theme(legend.position = 'none'); p2
p + p2

# do SVD
mat_centered = apply(mat, 2, function(x){(x - mean(x)) / sd(x)})
# cor_mat = cor(mat)
cor_mat2 = t(mat_centered) %*% mat_centered / nrow(mat_centered)
eig_res = eigen(cor_mat2)

p = list()
for(top_n in c(0, 1, 2, 5, 10, 20, 40, 50)) {
  pve = sum(eig_res$values[1:top_n]) / sum(eig_res$values)
  top_n_pc = mat_centered %*% eig_res$vectors[, 1 : top_n]
  tmp = qr(top_n_pc)
  Q_ = qr.Q(tmp)
  
  # mat_res = mat - Q_ %*% (t(Q_) %*% mat)
  mat_res2 = mat_centered - Q_ %*% (t(Q_) %*% mat_centered)
  
  plist = c()
  clist = c()
  for(n in 1:10) {
    select_idx = sample(ncol(mat_res2), 1); y = rnorm(nrow(mat_res2)) * 20 + rowSums(mat_res2[, select_idx, drop = F])
    # message(select_idx)
    res = fast_linear_regression(as.numeric(y), mat_res2, covariate = matrix(1, nrow = nrow(mat_res2), ncol = 1))
    plist = c(plist, res$pval)
    clist = c(clist, rep(paste0('repeat', n), ncol(mat_res2)))
  }
  
  p[[length(p) + 1]] = qqplot_by_group(plist, clist) + th + theme(legend.position = 'none') + ggtitle(paste0('nPC = ', top_n, ' PVE = ', signif(pve, digits = 3)))
  
  
}

out = do.call(grid.arrange, c(p, ncol = 4))
ggsave('export_mat_fac.png', out, width = 15, height = 8)
