# setwd('misc_data/supplementary_materials/')
library(dplyr)

annot_t1 = readRDS('../process_t1/t1_meta.rds')
normalized = !is.na(stringr::str_match(as.character(annot_t1$Field), '\\(normalised for head size\\)$'))
annot_t1 = annot_t1[ !normalized, ]

annot_dmri = readRDS('../download_some_matching_files/annot_dmri_idps.rds')

# we want the following columns in the final table
# ukb_field, t1_or_dmri, anatomy, left_or_right, notes, measurement_type, ukb_link

# working on T1
sub = annot_t1 %>% select(FieldID, modality, position, lr, Notes, matter_type, Link)
colnames(sub) = c('ukb_field', 't1_or_dmri', 'anatomy', 'left_or_right', 'notes', 'measurement_type', 'ukb_link')
# these are not for a specific anatomical region
sub$anatomy[ sub$ukb_field < 25011 ] = NA
# measurement_type == NA -> FIRST subcortical region
# measurement_type == grey matter -> FAST segmentation
# measurement_type for ukb_field < 25011 -> NA (neither FIRST nor FAST)
sub$measurement_type[ is.na(sub$measurement_type) ] = 'Subcortical volumes (FIRST)'
sub$measurement_type[ sub$measurement_type == 'grey matter' ] = 'Regional grey matter volumes (FAST)'
sub$measurement_type[ sub$ukb_field < 25011 ] = NA
# we need to exclude 25025 from "Regional grey matter volumes (FAST)" since its classification is FAST in brain_mri.pdf but not so in the ukb website
sub$measurement_type[ sub$ukb_field == 25025 ] = NA
sub$measurement_type[ is.na(sub$measurement_type) ] = 'Other T1 measurements'

# t1 visualization category
df_vis_t1 = readRDS('../vis_data_t1.rds')$table
# add a visualization classification
cate = sub$measurement_type
df_vis_t1 = inner_join(sub %>% select(ukb_field), df_vis_t1 %>% mutate(FieldID = as.numeric(FieldID)), by = c('ukb_field' = 'FieldID'))
cate[ df_vis_t1$db_name == 'cere' ] = 'Cerebellum'
cate[ cate == 'Regional grey matter volumes (FAST)' ] = 'Cortical'
cate[ cate == 'Subcortical volumes (FIRST)' ] = 'Subcortical'
cate[ cate == 'Other T1 measurements' ] = 'Total'
sub$t1_anatomy_group = cate


# working on dMRI
sub2 = annot_dmri %>% select(FieldID, modality, position, lr, Notes, type, Link, measure)
# all columns as above plus dmri_measure
colnames(sub2) = c('ukb_field', 't1_or_dmri', 'anatomy', 'left_or_right', 'notes', 'measurement_type', 'ukb_link', 'dmri_measure')
# measurement_type == 'mean' -> dMRI skeleton (TBSS-style measurement)
# measurement_type == 'weighted_mean_in_tract' -> dMRI weighted means (probabilistic-tractography-based measurement)
sub2$measurement_type[ sub2$measurement_type == 'mean'] = 'dMRI skeleton (TBSS-style measurement)'
sub2$measurement_type[ sub2$measurement_type == 'weighted_mean_in_tract'] = 'dMRI weighted means (probabilistic-tractography-based measurement)'

# combine the two
# add extra column dmri_measure for sub (T1) 
sub$dmri_measure = NA
sub2$t1_anatomy_group = NA
df = rbind(sub, sub2[, colnames(sub)])

final_order = c('ukb_field', 't1_or_dmri', 'anatomy', 'left_or_right', 'measurement_type', 'dmri_measure', 't1_anatomy_group', 'notes', 'ukb_link')
df = df[, final_order]

# ad hoc fix
# some cortical annotations are wrong
# see https://docs.google.com/spreadsheets/d/1Jhxw-DkN8kusdnG7eZAXecvJdLJEddefFKMBIDouC4k/edit#gid=0
df_fix = read.csv('supp_table_1_t1_fix.csv')
df_fix = df_fix %>% filter(X.1 != '')
df_fix$X.1 = as.character(df_fix$X.1)
map = data.frame(x0 = c('brainstem', 'subcortical'), x1 = c('Brainstem', 'Subcortical'))
df_fix = df_fix %>% inner_join(map, by = c('X.1' = 'x0'))
df_fix$x1 = as.character(df_fix$x1)
df = df %>% left_join(df_fix %>% select(ukb_field, x1), by = 'ukb_field')
df$t1_anatomy_group[ !is.na(df$x1) ] = df$x1[ !is.na(df$x1) ]
df = df %>% select(-x1)

df = df %>% filter(is.na(dmri_measure) | dmri_measure %in% c('FA', 'ICVF', 'ISOVF', 'OD'))

write.table(df, 'supp_table_1.tsv', quote = F, row.names = F, sep = '\t')

