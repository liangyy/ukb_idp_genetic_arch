load_gcta_hsq = function(fn) {
  cmd = paste0('cat ', fn, ' | grep "V(G)/Vp"')
  # system(cmd)
  tmp = tryCatch({
    data.table::fread(cmd = cmd, data.table = F, header = F)
  }, warning = function(error) {
    c(NA, NA, NA)
  }
  )
  # if(is.null(tmp)) {
  #   tmp = c(NA, NA, NA)
  # }
  data.frame(h2 = as.numeric(tmp[2]), h2_SE = as.numeric(tmp[3]))
}
