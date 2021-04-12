# form the table
library(oro.nifti)
library(dplyr)
options(stringsAsFactors = F)
source('rlib.R')

outdir = 'bxcan_vis'
dir.create(outdir)

idps = read.delim2('../supplementary_materials_3rd/supp_table_1.tsv')
idps$left_or_right[is.na(idps$left_or_right)] = 'NA'
idps = idps %>% mutate(IDP = paste0('IDP-', ukb_field))

t1_tags = c('Subcortical-vol', 'Cortical', 'Cerebellum')
dmri_tag = 'TBSS'
dmri_measures = c('ICVF', 'ISOVF', 'FA', 'OD')

message('--------------------- T1 IDPs --------------------')
results = list()
categories = list()
ii = 1

measurement_type = 'Subcortical volumes (FIRST)'
t1_anatomy_group = 'Subcortical'
tmp = prep_data(load_first, idps, measurement_type, t1_anatomy_group)
df_subcortical_vol = tmp$df
categories[[t1_tags[ii]]] = list(
  IDP = tmp$IDPs,
  vis_rds = paste0(t1_tags[ii], '.rds'),
  slide_position = c(2, 2, 2)
)
results[[t1_tags[ii]]] = df_subcortical_vol
ii = ii + 1


measurement_type = 'Regional grey matter volumes (FAST)'
t1_anatomy_group = 'Cortical'
tmp = prep_data(load_ho, idps, measurement_type, t1_anatomy_group)
df_cortical = tmp$df
categories[[t1_tags[ii]]] = list(
  IDP = tmp$IDPs,
  vis_rds = paste0(t1_tags[ii], '.rds'),
  slide_position = c(2, 2, 2)
)
results[[t1_tags[ii]]] = df_cortical
ii = ii + 1


measurement_type = 'Regional grey matter volumes (FAST)'
t1_anatomy_group = 'Cerebellum'
tmp = prep_data(load_cere, idps, measurement_type, t1_anatomy_group)
df_cerebellum = tmp$df
categories[[t1_tags[ii]]] = list(
  IDP = tmp$IDPs,
  vis_rds = paste0(t1_tags[ii], '.rds'),
  slide_position = c(2, 4, 3)
)
results[[t1_tags[ii]]] = df_cerebellum


for(tag in names(results)) {
  outfile = paste0(outdir, '/', tag, '.rds')
  # df_color = add_idp(results[[tag]]$label, idps, idp_col = 't1_anatomy_group')
  df_color = results[[tag]]$label
  message('Number of labels', ' (tag = ', tag, ') saved = ', nrow(df_color))
  saveRDS(list(df_color = df_color, img = results[[tag]]$img), outfile)
}

message('--------------------- dMRI IDPs --------------------')
measurement_type_ = 'dMRI skeleton (TBSS-style measurement)'
t1_anatomy_group = NULL
tmp = prep_data(load_tbss, idps, measurement_type_, t1_anatomy_group)
df_tbss = tmp$df
for(dm in dmri_measures) {
  idp_sub = idps %>% 
    filter(measurement_type == measurement_type_, dmri_measure == dm) %>% 
    pull(IDP)
  idp_ext = tmp$IDPs[tmp$IDPs %in% idp_sub]
  if(length(idp_ext) != length(idp_sub)) {
    message('dMRI wrong matching ', length(idp_ext), ' != ', length(idp_sub))
  }
  categories[[paste0(dmri_tag, '-', dm)]] = list(
    IDP = idp_ext,
    vis_rds = paste0(dmri_tag, '.rds'),
    slide_position = c(2, 2, 2)
  )
}

outfile = paste0(outdir, '/', dmri_tag, '.rds')
# df_color = add_idp(
#   df_tbss$label, 
#   idps, 
#   idp_col = NULL
# )
df_color = df_tbss$label
message('Number of labels saved = ', nrow(df_color))
saveRDS(list(df_color = df_color, img = df_tbss$img), outfile)



message('----------------- background image -----------------')
metadata_bg = readNIfTI('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
outfile = paste0(outdir, '/bg_img.rds')
saveRDS(metadata_bg, outfile)


message('---------------- saving categories -----------------')
outfile = paste0(outdir, '/meta_plot.rds')
saveRDS(categories, outfile)



