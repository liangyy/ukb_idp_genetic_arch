library(dplyr)

outdir = 'bxcan_vis'
dir.create(outdir)

load_idp_annot = function() {
  idp_annot = read.csv(paste0(outdir, '/idp_meta_data_tmp.csv')) # %>% select(-ukb_link)
  pc_annot = rbind(
    pc2df(yaml::read_yaml('../../misc_data/fourth_round_idp_preprocessing/output/dmri_covar_for_pc.yaml')) %>% 
      mutate(t1_or_dmri = 'dMRI', ukb_link = NA),
    pc2df(yaml::read_yaml('../../misc_data/fourth_round_idp_preprocessing/output/t1_covar_for_pc.yaml')) %>% 
      mutate(t1_or_dmri = 'T1', ukb_link = NA)
  )
  rbind(
    idp_annot,
    pc_annot[, colnames(idp_annot)]
  )
}

pc_map = function(x) {
  maplist = list(
    `Subcortical_vol` = "Subcortical",
    Cortical = "Gray-Cortical",
    Subcortical_GMvol = "Gray-Subcortical",
    Cerebellum = "Gray-Cerebellum"
  )
  if(x %in% names(maplist)) {
    return(maplist[[x]])
  } else {
    kk = strsplit(x, '-')[[1]]
    tag = kk[1]
    type = kk[2]
    if(type == 'ProbTrack') {
      return(paste0('w-', tag))
    } else {
      return(tag)
    }
  }
}

pc2df = function(mylist) {
  dd = list()
  for(n in names(mylist)) {
    dd[[length(dd) + 1]] = data.frame(IDP = mylist[[n]]$x, subtype = pc_map(n), left_or_right = NA, region = 'PC', notes = mylist[[n]]$x)
  }
  do.call(rbind, dd)
}

kk = load_idp_annot()
write.table(kk, paste0(outdir, '/idp_meta_data.csv'), 
quo = T, row = F, sep = ',', col = T)