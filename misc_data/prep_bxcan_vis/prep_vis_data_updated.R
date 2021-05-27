# this is an updated version of prep_vis_data.R 
# prep_vis_data2.R and prep_vis_data3.R are not affected by this update
# this update is due to the bad slicing which will miss some IDPs
# with this, we update the slide_position in meta data
# see details in check_slice.R for how we come up with the new slicing
# new slide_position = data.frame(dim = c(...), idx = c(...))
fix_dims = list(
  `Subcortical-vol` = data.frame(dim = 1:3, idx = c(100, 115, 60)),
  `Subcortical-GMvol` = data.frame(dim = 1:3, idx = c(100, 115, 60)),
  `Cortical` = data.frame(dim = c(1:3, 1:3), idx = c(c(50, 50, 45), c(130, 120, 85))),
  `Cerebellum` = data.frame(dim = 1:3, idx = c(91, 44, 30)),
  `TBSS` = data.frame(dim = c(1:3, 1:3), idx = c(c(60, 109, 91), c(120, 130, 40)))
  # `TBSS` = c(120, 130, 40)
)

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

t1_tags = c('Subcortical-vol', 'Subcortical-GMvol', 'Cortical', 'Cerebellum')
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
  slide_position = fix_dims[[t1_tags[ii]]],
  full_name = 'Subcortical Total Volume',
  type = 'T1'
)
results[[t1_tags[ii]]] = df_subcortical_vol
ii = ii + 1


measurement_type = 'Regional grey matter volumes (FAST)'
t1_anatomy_group = 'Subcortical'
map = data.frame(position = 'accumbens', new_name = 'ventral striatum')
tmp = prep_data(load_first, idps, measurement_type, t1_anatomy_group, map)
df_subcortical_gm = tmp$df
categories[[t1_tags[ii]]] = list(
  IDP = tmp$IDPs,
  vis_rds = paste0(t1_tags[ii], '.rds'),
  slide_position = fix_dims[[t1_tags[ii]]],
  full_name = 'Subcortical Gray Matter Volume',
  type = 'T1'
)
results[[t1_tags[ii]]] = df_subcortical_gm
ii = ii + 1


measurement_type = 'Regional grey matter volumes (FAST)'
t1_anatomy_group = 'Cortical'
tmp = prep_data(load_ho, idps, measurement_type, t1_anatomy_group)
df_cortical = tmp$df
categories[[t1_tags[ii]]] = list(
  IDP = tmp$IDPs,
  vis_rds = paste0(t1_tags[ii], '.rds'),
  slide_position = fix_dims[[t1_tags[ii]]],
  full_name = 'Cortical Gray Matter Volume',
  type = 'T1'
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
  slide_position = fix_dims[[t1_tags[ii]]],
  full_name = 'Cerebellum Gray Matter Volume',
  type = 'T1'
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
dmri_full_name_map = list(
  FA = 'Fractional Anisotropy',
  ICVF = 'Intra-Cellular Volume Fraction',
  ISOVF = 'Isotropic or free water Volume Fraction',
  OD = 'Orientation Dispersion Index'
)
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
    slide_position = fix_dims[[dmri_tag]],
    full_name = dmri_full_name_map[[dm]],
    type = 'dMRI'
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



