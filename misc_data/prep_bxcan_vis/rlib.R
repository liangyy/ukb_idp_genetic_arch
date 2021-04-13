add_idp = function(df_label, df_idp, idp_col = 't1_anatomy_group', map = NULL) {
  if(!is.null(map)) {
    tmp = left_join(df_label %>% select(position), map, by = 'position')
    df_label$position[!is.na(tmp$new_name)] = tmp$new_name[ !is.na(tmp$new_name) ]
  }
  kk = inner_join(df_label, df_idp, by = c('measurement_type', idp_col, 'position' = 'anatomy', 'lr' = 'left_or_right'))
  kk = kk %>% mutate(IDP = paste0('IDP-', ukb_field))
  kk = kk %>% distinct()
  if(nrow(df_label) > nrow(kk)) {
    message('The IDP matching has some issues. Please fix before proceeding.')
    return(NULL)
  }
  return(kk %>% select(color_code, IDP))
}

to_str_na = function(ss) {
  ss[is.na(ss)] = 'NA'
  ss
}

is_same = function(l1, l2) {
  n1 = length(intersect(l1, l2))
  n2 = length(union(l1, l2))
  if(n1 == n2 & n1 == length(l1) & n1 == length(l2)) {
    return(T)
  } else {
    return(F)
  }
}

prep_data = function(load_func, idps, measurement_type_, t1_anatomy_group_ = NULL, map = NULL) {
  message(measurement_type_, ' ', t1_anatomy_group_)
  dd = load_func()
  dd$label = dd$label %>% 
    mutate(
      measurement_type = measurement_type_, 
      t1_anatomy_group = t1_anatomy_group_,
      lr = to_str_na(lr)
    )
  message('Number of IDPs in vis = ', nrow(dd$label))
  IDPs = idps %>% filter(measurement_type == measurement_type_) 
  if(!is.null(t1_anatomy_group_)) {
    IDPs = IDPs %>% filter(t1_anatomy_group == t1_anatomy_group_)
  }
  IDPs = IDPs %>% pull(IDP)
  message('Number of IDPs in meta = ', length(IDPs))
  dd$label = add_idp(dd$label, idps, idp_col = if(! is.null(t1_anatomy_group_)) {'t1_anatomy_group' } else {NULL}, map = map)
  if(!is_same(dd$label$IDP, IDPs)) {
    message('Something is wrong')
  }
  list(df = dd, IDPs = IDPs)
}

# message('Subcortical total volumes (based on FIRST)')
load_first = function() {
  # df_first = read.table('~/Desktop/tmp/ukb_image/subcortical_labels.txt', comment.char = '#', sep = ',')
  df_first = read.table('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/subcortical_labels.txt', comment.char = '#', sep = ',')
  stem = 'Brain-Stem /4th Ventricle'
  # tmp = df_first[ df_first$V2 == stem, ]
  df_first = df_first[ df_first$V2 != stem, ]
  lr = unlist(lapply(strsplit(as.character(df_first$V2), '-'), function(x) {tolower(x[1])}))
  position = unlist(lapply(strsplit(as.character(df_first$V2), '-'), function(x) {paste0(tolower(x[-1][1]), collapse = ' ')}))
  df_first = cbind(df_first, data.frame(position = position, lr = lr))
  # df_first = rbind(
  #   df_first, 
  #   cbind(tmp, data.frame(position = 'brain stem + 4th ventricle', lr = NA))
  # )
  colnames(df_first)[1:2] = c('color_code', 'name_in_db')
  rownames(df_first) = NULL
  # kk3 = readNIfTI('~/Desktop/tmp/ukb_image/output_name_all_fast_firstseg.nii.gz')
  kk3 = readNIfTI('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/output_name_all_fast_firstseg.nii.gz')
  list(label = df_first, img = kk3)
}

# message('Cortical (based on FAST (Harvard + Oxford))')
load_ho = function() {
  # json_file <- "/Users/yanyul/Documents/repo/github/brain-coloring/acr2full.json"
  json_file <- "/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/acr2full.json"
  json_data <- jsonlite::fromJSON(paste(readLines(json_file), collapse=""))
  df_tmp = data.frame(name = names(json_data), full_name = unlist(json_data))
  # json_file <- "/Users/yanyul/Documents/repo/github/brain-coloring/labelmapper"
  json_file <- "/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/labelmapper"
  json_data2 <- jsonlite::fromJSON(paste(readLines(json_file), collapse=""))
  df_c = data.frame(idx = 1 : length(json_data2), name = json_data2)
  df_c = left_join(df_c, df_tmp, by = 'name')
  df_c$idx = df_c$idx - 1
  colnames(df_c)[1] = 'color_code'
  colnames(df_c)[2] = 'id_in_db'
  colnames(df_c)[3] = 'name_in_db'
  df_c = df_c[2:49, ]
  df_c$position = tolower(df_c$name_in_db)
  lr = c(rep('right', nrow(df_c)), rep('left', nrow(df_c)))
  code = df_c$color_code
  old_code_right = df_c$color_code
  old_code_left = df_c$color_code
  df_c$color_code = df_c$color_code * 2
  df_c2 = df_c
  df_c2$color_code = df_c2$color_code - 1
  df_c = rbind(df_c, df_c2)
  df_c$lr = lr
  rownames(df_c) = NULL
  # source: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;7ae2372b.1304
  # kk1 = readNIfTI('~/Downloads/HarvardOxford-Cortical-Lateralized-20130419/HarvardOxford/HarvardOxford-cortl-maxprob-thr25-1mm.nii.gz')
  kk1 = readNIfTI('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/HarvardOxford-cortl-maxprob-thr25-1mm.nii.gz')
  # tagsXML = XML::xmlParse('~/Downloads/HarvardOxford-Cortical-Lateralized-20130419/HarvardOxford-Cortical-Lateralized.xml')
  tagsXML = XML::xmlParse('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/HarvardOxford-Cortical-Lateralized.xml')
  tmp = XML::xmlToList(tagsXML)
  df = list()
  for(i in 1 : length(tmp$data)) {
    label = tmp$data[[i]]$text
    color_code = as.numeric(tmp$data[[i]]$.attrs['index']) + 1
    df[[length(df) + 1]] = data.frame(label = label, color_code)
  }
  df = do.call(rbind, df)
  df$lr = unlist(lapply(strsplit(df$label, ' '), function(x) { tolower(x[1]) }))
  df$position = unlist(lapply(strsplit(df$label, ' '), function(x) { paste0(tolower(x[-1]), collapse = ' ') }))
  df = df[, c('color_code', 'label', 'position', 'lr')]
  colnames(df)[2] = 'name_in_db'
  list(label = df, img = kk1)
}

# message('Cerebellum (based on the Diedrichsen cerebellar atlas)')
load_cere = function() {
  # annot = read.table('~/Downloads/Cerebellum-MNIflirt-MRICroN/Cerebellum-MNIflirt.nii.txt')
  annot = read.table('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/Cerebellum-MNIflirt.nii.txt')
  lr = unlist(lapply(strsplit(annot$V2, '_'), function(x) {tolower(x[1])}))
  # lr[ lr == 'vermis' ] = NA
  position = unlist(lapply(strsplit(annot$V2, '_'), function(x) {paste0(tolower(x[-1]), collapse = ' ')}))
  position = stringr::str_replace(position, 'crus', 'crus ')
  position = stringr::str_replace(position, 'i iv', 'i-iv')
  position = paste0(position, ' cerebellum')
  df = data.frame(color_code = annot$V1, name_in_db = annot$V2, position = position, lr = lr)
  # kk3 = readNIfTI('/usr/local/fsl/data/atlases/Cerebellum/Cerebellum-MNIflirt-maxprob-thr0-1mm.nii.gz')
  kk3 = readNIfTI('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/Cerebellum-MNIflirt-maxprob-thr0-1mm.nii.gz')
  list(label = df, img = kk3)
}

# message('TBSS based')
load_tbss = function() {
  kk = readNIfTI('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz')
  # orthographic(kk)
  # table(kk)
  labels = XML::xmlParse('/usr/local/fsl/data/atlases/JHU-labels.xml')
  tmp = XML::xmlToList(labels)
  df = list()
  for(i in 1 : length(tmp$data)) {
    label = tmp$data[[i]]$text
    color_code = as.numeric(tmp$data[[i]]$.attrs['index']) 
    df[[length(df) + 1]] = data.frame(label = label, color_code)
  }
  df = do.call(rbind, df)
  df = df[-1, ]
  lr = unlist(lapply(strsplit(df$label, ' '), function(x) {
    if(x[length(x)] == '') {
      x = x[-length(x)]
    }
    tmp = x[length(x)]
    if(tmp == 'L') {
      return('left')
    } else if(tmp == 'R') {
      return('right')
    } else {
      return(NA)
    }
  }))
  to_remove = c(' \\(include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus\\)', ' \\(include optic radiation\\)', ' \\(column and body of fornix\\)', ' \\(could be a part of anterior internal capsule\\)', ' \\(can not be resolved with current resolution\\)', ' \\(a part of MCP\\)')
  to_change = data.frame(target = c('\\) / ', '\\(', '\\)'), to = c('+', '', ''))
  position = unlist(lapply(strsplit(df$label, ' '), function(x) {
    tmp = x
    if(x[length(x)] == '') {
      x = x[-length(x)]
    }
    if(x[length(x)] %in% c('L', 'R')) {
      tmp = x[-length(x)]
    }
    tmp = paste0(tmp, collapse = ' ')
    for(pat in to_remove) {
      tmp = stringr::str_replace(tmp, pat, '')
    }
    for(i in 1 : nrow(to_change)) {
      tmp = stringr::str_replace(tmp, to_change$target[i], to_change$to[i])
    }
    tolower(tmp)
  }))
  df = cbind(df, data.frame(position = position, lr = lr))
  df = df[, c('color_code', 'label', 'position', 'lr')]
  colnames(df)[2] = 'name_in_db'
  list(label = df, img = kk)
}
