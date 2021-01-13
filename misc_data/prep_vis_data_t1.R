# form the table
library(oro.nifti)
library(dplyr)
options(stringsAsFactors = F)
## FIRST
load_first = function() {
  # df_first = read.table('~/Desktop/tmp/ukb_image/subcortical_labels.txt', comment.char = '#', sep = ',')
  df_first = read.table('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/subcortical_labels.txt', comment.char = '#', sep = ',')
  stem = 'Brain-Stem /4th Ventricle'
  tmp = df_first[ df_first$V2 == stem, ]
  df_first = df_first[ df_first$V2 != stem, ]
  lr = unlist(lapply(strsplit(as.character(df_first$V2), '-'), function(x) {tolower(x[1])}))
  position = unlist(lapply(strsplit(as.character(df_first$V2), '-'), function(x) {paste0(tolower(x[-1][1]), collapse = ' ')}))
  df_first = cbind(df_first, data.frame(position = position, lr = lr))
  df_first = rbind(
    df_first, 
    cbind(tmp, data.frame(position = 'brain stem + 4th ventricle', lr = NA))
  )
  colnames(df_first)[1:2] = c('color_code', 'name_in_db')
  rownames(df_first) = NULL
  # kk3 = readNIfTI('~/Desktop/tmp/ukb_image/output_name_all_fast_firstseg.nii.gz')
  kk3 = readNIfTI('/Users/yanyul/Desktop/tmp/ukb_image/data_for_t1/output_name_all_fast_firstseg.nii.gz')
  list(label = df_first, img = kk3)
}
df_first = load_first()

## FAST (Harvard + Oxford)
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
df_ho = load_ho()

## the Diedrichsen cerebellar atlas
load_cerebellar = function() {
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
df_cere = load_cerebellar()

tmp = rbind(
  df_first$label %>% mutate(db_name = 'first', matter_type = 'NA'), 
  df_ho$label %>% mutate(db_name = 'ho', matter_type = 'grey matter'), 
  df_cere$label %>% mutate(db_name = 'cere', matter_type = 'grey matter')
)
idp = readRDS('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/process_t1/t1_meta.rds')
idp$lr[is.na(idp$lr)] = 'NA'
tmp$lr[is.na(tmp$lr)] = 'NA'
idp$matter_type[is.na(idp$matter_type)] = 'NA'
tmp = left_join(idp, tmp, by = c('position', 'lr', 'matter_type'))

kk = readNIfTI('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')

saveRDS(list(table = tmp, db = list(first = df_first$img, ho = df_ho$img, cere = df_cere$img), bg = kk), 'vis_data_t1.rds')



