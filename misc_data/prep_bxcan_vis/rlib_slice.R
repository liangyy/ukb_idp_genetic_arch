slice_3d = function(tensor, dim_, idx) {
  if(dim_ == 1) {
    return(tensor[idx, , ])
  } else if(dim_ == 2) {
    return(tensor[, idx, ])
  } else if(dim_ == 3) {
    return(tensor[, , idx])
  } 
}
get_range = function(img3d) {
  dd = list()
  for(d in 1 : 3) {
    for(i in 1 : dim(img3d)[d]) {
      dd[[length(dd) + 1]] = data.frame(dim = d, i = i, code = unique(as.numeric(slice_3d(img3d, dim_ = d, idx = i))))
    }
  }
  dd = do.call(rbind, dd)
  list(range = dd %>% group_by(dim, code) %>% summarize(min_i = min(i), max_i = max(i)) %>% ungroup(), dd = dd)
}
check_slice = function(img3d, idxs, dims = c(1, 2, 3)) {
  dd = list()
  code_all = unique(as.numeric(img3d))
  total = length(code_all)
  for(i in 1 : length(dims)) {
    d = dims[i]
    idx = idxs[i]
    dd[[length(dd) + 1]] = data.frame(dim = d, i = idx, code = unique(as.numeric(slice_3d(img3d, dim_ = d, idx = idx))))
  }
  dd = do.call(rbind, dd)
  tmp = dd %>% group_by(dim, i) %>% summarize(n = length(unique(code))) %>% ungroup()
  code_appear = unique(dd$code)
  tmp = rbind(
    tmp, data.frame(dim = 'all', i = 'all', n = total),
    data.frame(dim = 'slice', i = 'slice', n = length(code_appear))
  )
  # print(code_all)
  # print(code_appear)
  list(report = tmp, missing = code_all[! code_all %in% code_appear])
}
