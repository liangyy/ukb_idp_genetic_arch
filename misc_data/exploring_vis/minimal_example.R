# Some dependent libraries
library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(theme_classic(base_size = 10))

# Load helper functions
# Misc helper functions
rlib_file = '/Users/yanyul/Documents/repo/github/ukb_idp_genetic_arch/rmd/rlib.R'
# Main functions for visualization
rlib_file_for_vis = '/Users/yanyul/Documents/repo/github/ukb_idp_genetic_arch/rmd/rlib_vis.R'
source(rlib_file)
source(rlib_file_for_vis)

# The S-BrainXcan results to visualize
sbrainxcan_file = '/Users/yanyul/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/dmri.original.all_covar.w_pc.gw_elastic_net_beta_x_SCZ_PGC_2020_x_simagexcan.csv'

# The file containing some meta information of IDPs: 
#   IDP ID, information of the measurement, etc
idp_meta_file = '/Users/yanyul/Documents/repo/github/ukb_idp_genetic_arch/misc_data/download_some_matching_files/annot_dmri_idps.rds'

# To visualize dMRI, we have the following types of measurements: 
# FA, MD, MO, OD, L1, L2, L3, ICVF, ISOVF
# For visualization, we need to select one type to visualize.
measure_to_show = 'ICVF'

# The dataset required for visualizing the brain
vis_data_file = '/Users/yanyul/Documents/repo/github/ukb_idp_genetic_arch/misc_data/vis_data_dmri.rds'

# Set output directory
out_prefix = '/Users/yanyul/Documents/repo/github/ukb_idp_genetic_arch/misc_data/exploring_vis/minimal.'

# Load S-BrainXcan data
df = read.csv(sbrainxcan_file) %>% mutate(zscore = p2z(pval, bhat))

# Load IDP meta information
idp_meta = readRDS(idp_meta_file) %>% mutate(IDP = paste0('IDP-', FieldID))

# Extract the set of IDPs to visualize 
df_sub = df[ df$IDP %in% (idp_meta %>% filter(measure == measure_to_show) %>% pull(IDP)), ]

# Load visualization dataset
vis_data = readRDS(vis_data_file)
vis_data$table$IDP = paste0('IDP-', vis_data$table$FieldID)

# Generate figures using vis_assoc function
# for dMRI: set type = 'tbss'
# for T1: set type = 'cere' (for cerebellum regions), 'first' (for subcortical regions), 'ho' (for the cortical regions) 
types = c('tbss')
scores = c('zscore', 'pip')
for(s in scores) {
  for(i in types) {
    message('Working on ', i, ' ', s)
    q = vis_assoc(vis_data, df_sub, type = i, score = s)
    ggsave(paste0(out_prefix, s, '_', i, '.', measure_to_show, '.png'), q, height = 4, width = 10)
  }
}

# Visualize the regions using vis_region
r1 = vis_region(vis_data, df_sub, type = types[1])
ggsave(paste0(out_prefix, 'region.tbss.', measure_to_show, '.png'), r1, height = 6, width = 13)

