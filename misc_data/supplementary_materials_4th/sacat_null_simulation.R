# setwd('misc_data/supplementary_materials_4th/')
library(dplyr)
library(ggplot2)
theme_set(theme_classic(base_size = 15))
set.seed(2020)
devtools::source_gist("38431b74c6c0bf90c12f")

outdir = 'sacat'
dir.create(outdir)

acat_simple_signed = function(p_vec, sign_stat)
{
  TT = sum( sign(sign_stat) * tan( (0.5 - p_vec) *pi ) )
  ( .5 - atan(TT / length(p_vec)) / pi )
}
acat_two_dir2 = function(p_vec, sign_stat) {
  p1 = acat_simple_signed(p_vec, sign_stat)
  p2 = acat_simple_signed(p_vec, - sign_stat)
  # 1 - (1 - min(p1, p2)) ^ 2
  2 * min(p1, p2)
}
n = 10000
k = 20
corr_between_test = c(0, 0.1, 0.5)
png(paste0(outdir, '/sacat_null.png'), width = 6, height = 9, units = "in", res = 300)
par(mfrow = c(3, 2))
for(r in corr_between_test) {
  R = matrix(r, ncol = k, nrow = k) 
  diag(R) = 1
  R[1, 2:k] = -R[1, 2:k]
  R[2:k, 1] = -R[2:k, 1]
  Z = mvtnorm::rmvnorm(n, sigma = R)
  P = pnorm(abs(Z), lower.tail = F) * 2
  acat_s5 = apply(cbind(P, Z), 1, function(x) { 
    ss = length(x) / 2
    acat_two_dir2(x[1 : ss], x[(ss + 1) : (ss * 2)])
  })
  hist(acat_s5, xlab = 'p-value combined by SACAT', main = paste0('r = ', r))
  qqunif(acat_s5, main = paste0('r = ', r))
}
dev.off()
