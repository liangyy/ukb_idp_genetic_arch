# check slicing def
# the goal is to have all code being shown under the slicing

datadir = 'bxcan_vis'

source('rlib_slice.R')
meta_list = readRDS(paste0(datadir, '/meta_plot.rds'))
tags = names(meta_list)
tbss = F
for(tag in tags) {
  logging::loginfo(paste0('Tag = ', tag))
  title = meta_list[[tag]]$full_name
  tag = 'TBSS-ICVF'
  save_name = tag
  if(substr(tag, 1, 4) == 'TBSS') {
    if(isTRUE(tbss)) {
      next
    } else {
      tbss = T
      save_name = 'TBSS'
    }
  }
  vis = readRDS(paste0(datadir, '/', save_name, '.rds'))
  curr_dims_ = meta_list[[tag]]$slide_position
  curr_dims_ = floor(dim(vis$img) / curr_dims_)
  kk = get_range(vis$img)
  df = kk$range
  print(paste('---------', tag, '---------'))
  print(paste0('dim = ', paste0(curr_dims_, collapse = ', ')))
  print(df %>% group_by(dim) %>% summarize(max(min_i), min(max_i)))
  re = check_slice(vis$img, curr_dims_)
  tmp = df %>% filter(code %in% re$missing)
  print(df %>% filter(code %in% re$missing))
}

fix_dims = list(
  `Subcortical-vol` = c(100, 115, 60),
  `Subcortical-GMvol` = c(100, 115, 60),
  `Cortical` = c(c(50, 50, 45), c(130, 120, 85)),
  `Cerebellum` = c(91, 44, 30),
  `TBSS` = c(c(60, 109, 91), c(120, 130, 40))
  # `TBSS` = c(120, 130, 40)
)
tbss = F
for(tag in tags) {
  logging::loginfo(paste0('Tag = ', tag))
  title = meta_list[[tag]]$full_name
  save_name = tag
  if(substr(tag, 1, 4) == 'TBSS') {
    if(isTRUE(tbss)) {
      next
    } else {
      tbss = T
      save_name = 'TBSS'
    }
  }
  vis = readRDS(paste0(datadir, '/', save_name, '.rds'))
  curr_dims_ = fix_dims[[save_name]]
  dims = 1 : 3
  if(length(curr_dims_) > 3) {
    dims = c(dims, dims)
  }
  kk = get_range(vis$img)
  df = kk$range
  print(paste('---------', tag, '---------'))
  print(paste0('dim = ', paste0(curr_dims_, collapse = ', ')))
  print(df %>% group_by(dim) %>% summarize(max(min_i), min(max_i)))
  re = check_slice(vis$img, curr_dims_, dims)
  tmp = df %>% filter(code %in% re$missing)
  print(df %>% filter(code %in% re$missing))
}



