get_upper_tri = function(mat) {
  return(mat[upper.tri(mat)])
}
get_corr = function(mat) {
  tmp1 = cor(mat)
  df = data.frame(corr = get_upper_tri(tmp1)) 
  df
}