
parse_id = function(x) {
  lr_map = list(L = 'left', R = 'right', V = 'vermis')
  kk = stringr::str_remove(x, '^IDP_T1_FAST_ROIs_')
  lr_ = substr(kk, 1, 1)
  if(lr_ %in% names(lr_map)) {
    lr = lr_map[[lr_]]
    kk = stringr::str_remove(kk, paste0('^', lr_, '_'))
  } else {
    lr = NA
  }
  list(lr = lr, name = kk)
}

str_in = function(str, pat) {
  !is.na(stringr::str_match(str, paste0('^', pat)))
}

match_two_strlist = function(s1, s2) {
  match = T
  nl = min(length(s1), length(s2))
  for(j in 1 : nl) {
    if(!str_in(s1[j], s2[j])) {
      match = F
      break
    }
  }
  match
}

ana_tostr = function(xx) {
  # tolower(unlist(strsplit(pool[ll], ' ')))
  res = strsplit(xx, ', ')[[1]]
  if(length(res) > 1) {
    ll = res[2]
    ll = stringr::str_remove(ll, ' ')
  } else {
    ll = NULL
  }
  res = res[1]
  res = tolower(unlist(strsplit(res, ' ')))
  c(res, ll)
}

put_cere_in_last = function(xx) {
  o = c()
  last = NULL
  for(x in xx) {
    if(x != 'cerebellum') {
      o = c(o, x)
    } else {
      last = x
    }
  }
  c(o, last)
}

find_x = function(pool, name) {
  res = NA
  name_strs = put_cere_in_last(tolower(unlist(strsplit(name, '_'))))
  for(ll in 1 : length(pool)) {
    ana_strs = put_cere_in_last(ana_tostr(pool[ll]))
    match = match_two_strlist(ana_strs, name_strs)
    if(match == T) {
      res = pool[ll]
      break
    } 
  }
  res
}

mymatch = function(x, y, find_func = find_x) {
  # find matched y for each x
  # y[, 1]: anatomy; y[, 2] = left_or_right
  o = list()
  for(xi in x) {
    res = parse_id(xi)
    lr = res$lr
    name = res$name
    if(is.na(lr)) {
      pool_y = y[is.na(y[[2]]), ]
    } else {
      pool_y = y[!is.na(y[[2]]) & (y[[2]] == lr), ]
    }
    o[[length(o) + 1]] = data.frame(anatomy = find_func(pool_y[[1]], name), lr = lr)
  }
  do.call(rbind, o)
}

hard_find = function(pool, name) {
  res = NA
  hard_map = data.frame(
    id1 = c('inf_temp_gyrus_tempocc', 'mid_temp_gyrus_tempocc', 'latocc_cortex_sup', 'latocc_cortex_inf'), 
    id2 = c('inferior temporal gyrus, temporooccipital part', 'middle temporal gyrus, temporooccipital part', 'lateral occipital cortex, superior division', 'lateral occipital cortex, inferior division')
  )
  if(name %in% hard_map$id1) {
    res = pool[ pool == hard_map$id2[ hard_map$id1 == name ] ]
  } else {
    res = NA
  }
  res
}
