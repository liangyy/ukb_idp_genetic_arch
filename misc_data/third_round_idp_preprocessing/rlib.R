get_upper_tri = function(mat) {
  return(mat[upper.tri(mat)])
}
plot_image = function(mat) {
  p = mat %>% reshape2::melt() %>% ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradientn(limits = c(-1, 1), colors = c("blue", "white", "red")) + 
    theme(axis.text.x = element_blank())
  p
}
do_all = function(mat, skip_pca = F) {
  tmp1 = cor(mat)
  p1 = plot_image(tmp1)
  if(skip_pca == T) {
    message('No PCA')
    p2 = ggplot() + geom_text(x = 1, y = 1, label = 'placeholder')
    p3 = ggplot() + geom_text(x = 1, y = 1, label = 'placeholder')
    p4 = ggplot() + geom_text(x = 1, y = 1, label = 'placeholder')
    res = list(residual = standardize(mat), pc_loadings = NA, pve = NA, pc = NULL)
  } else {
    res = split_by_pca(mat, pve_cutoff = 0.1, skip = F)
    tmp2 = cor(res$residual)
    p2 = plot_image(tmp2)
    p4 = data.frame(pc_loading = res$pc_loadings, idp = colnames(mat)) %>% 
      ggplot() + geom_point(aes(x = pc_loading, y = idp))
    p3 = data.frame(before = get_upper_tri(tmp1), after = get_upper_tri(tmp2)) %>% 
      ggplot() + geom_histogram(aes(x = before, fill = 'before'), alpha = 0.5, binwidth = 0.1) +
      geom_histogram(aes(x = after, fill = 'after'), alpha = 0.5, binwidth = 0.1)
  }
  list(ps = list(p1, p2, p3, p4), res = res)
}
plot_all = function(ps, title) {
  p1 = ps[[1]] + ggtitle('Before PC adj')
  p2 = ps[[2]] + ggtitle('After PC adj')
  p3 = ps[[3]] + ggtitle('Histogram of corr')
  p4 = ps[[4]] + ggtitle('PC loading')
  ((p1 + p2) / (p3 + p4)) + plot_annotation(title = title)
}