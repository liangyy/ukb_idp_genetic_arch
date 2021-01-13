vis_assoc = function(vis_data, brainxcan_table, type, score = 'zscore') {
  # type = 'ho'
  sub = vis_data$table[ vis_data$table$IDP %in% brainxcan_table$IDP,  ]
  sub1 = sub[ !is.na(sub$db_name) & sub$db_name == type, ]
  img = vis_data$db[[type]]
  img2 = img
  img2[] = NA
  for(i in 1 : nrow(sub1)) {
    img2[img == sub1$color_code[i]] = brainxcan_table[[score]][ brainxcan_table$IDP == sub1$IDP[i] ] 
  }
  # tt0 = reshape2::melt(img) %>% filter(value != 0) %>% left_join(vis_data$table[ vis_data$table$db_name == type, c('color_code', 'position', 'lr')], by = c('value' = 'color_code'))
  tt = reshape2::melt(img2) %>% filter(value != 0)
  bb = reshape2::melt(vis_data$bg) %>% filter(value != 0) %>% rename(bg = value)
  bb = left_join(bb, tt, by = c('Var1', 'Var2', 'Var3'))
  mid1 = vis_data$bg@dim_[2]
  mid2 = vis_data$bg@dim_[3]
  mid3 = vis_data$bg@dim_[4]
  if(type == 'ho') {
    d1 = d2 = d3 = 2
  } else if(type == 'cere') {
    d1 = 2
    d2 = 4
    d3 = 3
  } else if(type == 'first') {
    d1 = d2 = d3 = 2
  } else if(type == 'tbss') {
    d1 = d2 = d3 = 2
  }
  
  tmp = rbind(
    bb %>% filter(Var1 == floor(mid1 / d1)) %>% mutate(direction = 'a1') %>% rename(x = Var2, y = Var3, ref = Var1),
    bb %>% filter(Var2 == floor(mid2 / d2)) %>% mutate(direction = 'a2') %>% rename(x = Var1, y = Var3, ref = Var2),
    bb %>% filter(Var3 == floor(mid3 / d3)) %>% mutate(direction = 'a3') %>% rename(x = Var1, y = Var2, ref = Var3)
  )
  
  p = tmp %>% ggplot() + 
    geom_raster(aes(x, y, fill = value)) +
    scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0, na.value = 'transparent') +
    geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
    scale_alpha(range = c(0.2, 0.8)) +
    coord_equal() + facet_wrap(~direction, labeller = label_both) +
    guides(color = guide_legend("-bg"), alpha = FALSE)
  
  # p1 = bb %>% filter(Var1 == floor(mid1 / d1)) %>% ggplot() + 
  #   geom_raster(aes(Var2, Var3, fill = value)) +
  #   scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #   geom_tile(aes(Var2, Var3, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # p2 = bb %>% filter(Var2 == floor(mid2 / d2)) %>% ggplot() + 
  #   geom_raster(aes(Var1, Var3, fill = value)) +
  #   scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #   geom_tile(aes(Var1, Var3, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # p3 = bb %>% filter(Var3 == floor(mid3 / d3)) %>% ggplot() + 
  #   geom_raster(aes(Var1, Var2, fill = value)) +
  #   scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #   geom_tile(aes(Var1, Var2, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # list(a1 = p1, a2 = p2, a3 = p3)
  p
}


vis_region = function(vis_data, brainxcan_table, type) {
  # type = 'ho'
  sub = vis_data$table[ vis_data$table$IDP %in% brainxcan_table$IDP,  ]
  sub1 = sub[ !is.na(sub$db_name) & sub$db_name == type, ]
  img = vis_data$db[[type]]
  tt0 = reshape2::melt(img) %>% filter(value != 0) %>% left_join(sub1[, c('color_code', 'position', 'lr')], by = c('value' = 'color_code'))
  bb = reshape2::melt(vis_data$bg) %>% filter(value != 0) %>% rename(bg = value)
  bb = left_join(bb, tt0, by = c('Var1', 'Var2', 'Var3'))
  mid1 = vis_data$bg@dim_[2]
  mid2 = vis_data$bg@dim_[3]
  mid3 = vis_data$bg@dim_[4]
  if(type == 'ho') {
    d1 = d2 = d3 = 2
  } else if(type == 'cere') {
    d1 = 2
    d2 = 4
    d3 = 3
  } else if(type == 'first') {
    d1 = d2 = d3 = 2
  } else if(type == 'tbss') {
    d1 = d2 = d3 = 2
  }
  # bb$pos_lr = paste(bb$position, bb$lr)
  # bb$pos_lr[is.na(bb$position)] = NA
  # mycolors = c(RColorBrewer::brewer.pal(name="Paired", n = 12), RColorBrewer::brewer.pal(name="Dark2", n = 8))
  
  tmp = rbind(
    bb %>% filter(Var1 == floor(mid1 / d1)) %>% mutate(direction = 'a1') %>% rename(x = Var2, y = Var3, ref = Var1),
    bb %>% filter(Var2 == floor(mid2 / d2)) %>% mutate(direction = 'a2') %>% rename(x = Var1, y = Var3, ref = Var2),
    bb %>% filter(Var3 == floor(mid3 / d3)) %>% mutate(direction = 'a3') %>% rename(x = Var1, y = Var2, ref = Var3)
  )
  
  p = tmp %>% ggplot() + 
    geom_raster(aes(x, y, fill = position)) +
    # scale_fill_manual(values = mycolors) +
    geom_tile(aes(x, y, alpha = -bg), fill = "grey20") +
    scale_fill_discrete(na.value = "transparent") +
    scale_alpha(range = c(0.2, 0.8)) +
    coord_equal() + facet_wrap(~direction, labeller = label_both) + theme(legend.position = 'bottom', legend.title = element_blank()) +
    guides(color = guide_legend("-bg"), alpha = FALSE)
  # 
  # p1 = bb %>% filter(Var1 == floor(mid1 / d1)) %>% ggplot() + 
  #   geom_raster(aes(Var2, Var3, fill = pos_lr)) +
  #   scale_fill_manual(values = mycolors) +
  #   geom_tile(aes(Var2, Var3, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # p2 = bb %>% filter(Var2 == floor(mid2 / d2)) %>% ggplot() + 
  #   geom_raster(aes(Var1, Var3, fill = value)) +
  #   scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #   geom_tile(aes(Var1, Var3, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # p3 = bb %>% filter(Var3 == floor(mid3 / d3)) %>% ggplot() + 
  #   geom_raster(aes(Var1, Var2, fill = value)) +
  #   scale_fill_gradient2(name = score, low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #   geom_tile(aes(Var1, Var2, alpha = -bg), fill = "grey20") +
  #   scale_alpha(range = c(0.2, 0.8)) +
  #   coord_equal() 
  # 
  # list(a1 = p1, a2 = p2, a3 = p3)
}
