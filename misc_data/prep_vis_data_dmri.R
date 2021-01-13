library(dplyr)
library(oro.nifti)
options(stringsAsFactors = F)
# 
dmri = readRDS('misc_data/download_some_matching_files/annot_dmri_idps.rds')
# pos1 = unique(dmri$position[dmri$type == 'weighted_mean_in_tract'])
# pos2 = unique(dmri$position[dmri$type == 'mean'])
# intersect(pos1, pos2)

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

df_tbss = load_tbss()

dmri_tbss = dmri %>% filter(type == 'mean') %>% select(position, lr) %>% distinct()


tmp = left_join(dmri_tbss, df_tbss$label, by = c('position', 'lr'))
tmp[(is.na(tmp$color_code)), ]

idp = readRDS('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/download_some_matching_files/annot_dmri_idps.rds')

tmp = rbind(
  df_tbss$label %>% mutate(db_name = 'tbss', type = 'mean')
)
idp$lr[is.na(idp$lr)] = 'NA'
tmp$lr[is.na(tmp$lr)] = 'NA'
idp$type[is.na(idp$type)] = 'NA'
tmp = left_join(idp, tmp, by = c('position', 'lr', 'type'))

kk = readNIfTI('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')

saveRDS(list(table = tmp, db = list(tbss = df_tbss$img), bg = kk), 'misc_data/vis_data_dmri.rds')
