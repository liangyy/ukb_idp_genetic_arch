# here, we prepare the meta data for generating the brainxcan report
library(dplyr)

outdir = 'bxcan_vis'
dir.create(outdir)

options(stringsAsFactors = F)
gannot2 = read.delim2('../supplementary_materials_4th/supp_table_1.tsv')
gannot2 = gannot2 %>% mutate(IDP = paste0('IDP-', ukb_field)) %>% rename(region = anatomy)
# gannot2$region = as.character(gannot2$region)
gen_substype = function(gannot2) {
  xx = rep(NA, nrow(gannot2))
  xx[gannot2$t1_or_dmri == 'T1'] = gannot2$t1_anatomy_group[gannot2$t1_or_dmri == 'T1']
  xx[gannot2$t1_or_dmri == 'dMRI'] = gannot2$dmri_measure[gannot2$t1_or_dmri == 'dMRI']
  xx[gannot2$t1_or_dmri == 'dMRI'] = gannot2$dmri_measure[gannot2$t1_or_dmri == 'dMRI']
  is_g = gannot2$measurement_type == 'Regional grey matter volumes (FAST)'
  xx[is_g] = paste0('Gray-', xx[is_g])
  is_w = gannot2$measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)'
  xx[is_w] = paste0('w-', xx[is_w])
  xx
}
gannot2$subtype = gen_substype(gannot2)

write.table(
  gannot2 %>% select(IDP, t1_or_dmri, subtype, left_or_right, region, notes, ukb_link), 
  paste0(outdir, '/idp_meta_data_tmp.csv'), 
  quo = T, row = F, sep = ',', col = T
)

color_red = "#D55E00"
color_yellow = "#F0E442"
color_green = "#009E73"
color_orange = "#E69F00"
color_lightblue = "#56B4E9"
color_blue = "#0072B2"
mycolor_map = list(
  "FA" = color_red, 
  "ISOVF" = color_blue,
  "OD" = color_green,
  "ICVF" = color_orange,
  # "w-OD",
  # "w-ICVF",
  # "w-ISOVF",
  "Subcortical" = color_orange,
  "Gray-Cerebellum" = color_blue,
  "Gray-Cortical" = color_blue,
  "Gray-Subcortical" = color_blue,
  "Gray-Brainstem" = color_blue,
  "PC-Cerebellum-1" = "black",
  "PC-Cortical-1" = "black",
  # "PC-FA-ProbTrack-1",
  "PC-FA-TBSS-1" = "black",
  # "PC-ICVF-ProbTrack-1",
  "PC-ICVF-TBSS-1" = "black",
  # "PC-ISOVF-ProbTrack-1",
  "PC-ISOVF-TBSS-1" = "black",
  # "PC-OD-ProbTrack-1",
  "PC-OD-TBSS-1" = "black",
  "PC-Subcortical_GMvol-1" = "black",
  "PC-Subcortical_vol-1" = "black",
  "Brainstem" = "gray",
  "Total" = "gray"
)

yaml::write_yaml(mycolor_map, paste0(outdir, '/report_color_code.yaml'))
