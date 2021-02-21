# Here I download some matching files provided by groups working intensively on UKB IDP (Nature 2018) 
# And I clean it up so that it matches with what we have for UKB ImageXcan project
# setwd('misc_data/download_some_matching_files/')
options(stringsAsFactors = F)
library(dplyr)
df_list = data.table::fread('cat list.txt | awk \'{print $1,$2,$3,$4}\'', data.table = F)
df_list$V1[19] = 19
df_info = readLines('list_last_col.txt')
df_info = stringr::str_remove(df_info, ' $')
df_idcol = readLines('list_id_col.txt')
df_list$info = df_info
df_list$n2018_id = df_idcol
df_list = df_list[ df_list$info != '', ]
df_ukb = read.delim2('~/Downloads/Data_Dictionary_Showcase.tsv')


remove_special_char = function(x, char = c('\\-', '_', '\\+', '\\(', '\\)', ',', 'noise', '\\.', 'the')) {
  for(cc in char) {
    # print(cc)
    x = stringr::str_remove(x, cc)
  }
  x
}

bag_of_char = function(s) {
  paste0(sort(strsplit(s, '')[[1]]), collapse = '')
}

str_washer = function(ss) {
  ss = lapply(strsplit(ss, ' '), function(x) { remove_special_char(x)})
  ss = unlist(lapply(ss, function(x) {paste0(tolower(x), collapse = '')}))
  ss = sapply(ss, function(x){remove_special_char(x, char = c('intaskfmridata$', 'fromdmridata$', 'fromt1brainimage$', '^total', 'calculatedwithbianca', 'brainimage$', 'images$'))})
  sapply(ss, function(x){bag_of_char(x)})
}

df_list$ss = str_washer(df_list$info)
length(unique(df_list$ss))
df_ukb$footprint = str_washer(df_ukb$Notes)
df_ukb$footprint2 = str_washer(df_ukb$Field)


df_ukb_sub = df_ukb[ df_ukb$footprint %in% df_list$ss, ]
df_ukb_sub$good_match = df_ukb_sub$footprint
df_ukb_sub2 = df_ukb[ df_ukb$footprint2 %in% df_list$ss, ]
df_ukb_sub2$good_match = df_ukb_sub2$footprint2
df_ukb_sub = rbind(df_ukb_sub, df_ukb_sub2); 
df_ukb_sub = df_ukb_sub[!duplicated(df_ukb_sub$good_match), ]
df_list[!df_list$ss %in% df_ukb_sub$good_match, ]
df_list = left_join(df_list, df_ukb_sub, by = c('ss' = 'good_match'))

modalities = c('T2_FLAIR', 'FS', 'SWI', 'T1', 'dMRI', 'rfMRI_CONN', 'rfMRI_AMP', 'tfMRI')
df_mod = list()
for(mm in modalities) {
  tmp = read.table(paste0(mm, '.txt'), header = F)
  df_mod[[length(df_mod) + 1]] = tmp %>% mutate(modality = mm)
}
df_mod = do.call(rbind, df_mod)
df_list = left_join(df_list, df_mod, by = 'V1')
table(df_list$modality)
colnames(df_list)[1:3] = c('code_in_their_paper', 'id_in_their_paper', 'id2_in_their_paper')
df_list = df_list %>% select(-footprint, -footprint2) %>% rename(footprint = ss)
df_list$modality = unlist(lapply(strsplit(df_list$id_in_their_paper, '_'), function(x){x[[1]]}))
df_list$modality[1:16] = 'covariate'
saveRDS(df_list, 'cleanup_annot_our_idps.rds')
