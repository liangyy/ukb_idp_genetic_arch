inv_norm = function(x, offset = 1, ...) {
  r = rank(x, ...)
  g = r / (length(r) + offset)
  o = qnorm(g)
  return(o)
}

regress_out_univariate = function(y, x) {
  y = y - mean(y)
  x = x - mean(x)
  b = sum(y * x) / sum(x * x)
  y - x * b
}

regress_out_matrix = function(y_mat, covar_mat) {
  x = apply(covar_mat, 2, function(k) {k - mean(k)})
  y = apply(y_mat, 2, function(k) {k - mean(k)})
  tmp = qr(x)
  Q = qr.Q(tmp)
  y - Q %*% (t(Q) %*% y)
}

standardize = function(x) {
  apply(x, 2, function(y) { (y - mean(y)) / sd(y) })
}

split_by_pca = function(x, pve_cutoff = 0.5, skip = F) {
  x = standardize(x)
  # x = apply(x, 2, inv_norm)
  if(skip == T) {
    resi = x
    loading = NULL
    pc_mat = NULL
    pve = NULL
  } else {
    res = svd(x)
    pve = cumsum(res$d^2 / sum(res$d^2))
    npc = sum(pve <= pve_cutoff) + 1
    pc_mat = res$u[, 1 : npc, drop = F]
    pc_mat = apply(pc_mat, 2, function(k) {k - mean(k)})
    loading = res$v[, 1 : npc, drop = F]
    resi = x - pc_mat %*% (t(pc_mat) %*% x)
    pve = pve[1 : npc]
  }
  list(residual = resi, pc = pc_mat, pc_loadings = loading, pve = pve)
}

myplot_magic = function(df) {
  mytmp = df 
  mytmp %>% ggplot() +
    theme_bw() + 
    geom_raster(aes(x = as.character(Var1), y = as.character(Var2), fill = value)) + scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    # annotate(geom = "text", x = mytmp_g$pos, y = -20, label = mytmp_g$matter_type.1, size = 4) + 
    coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
    theme(
      plot.margin = unit(c(0.1, 0.1, 2, 0.1), "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}
