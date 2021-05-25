load_idp_annot = function() {
  idp_annot = read.csv('../prep_bxcan_vis/bxcan_vis/idp_meta_data.csv') %>% select(-ukb_link)
  # pc_annot = rbind(
  #   pc2df(yaml::read_yaml('../../misc_data/fourth_round_idp_preprocessing/output/dmri_covar_for_pc.yaml')) %>% 
  #     mutate(t1_or_dmri = 'dMRI'),
  #   pc2df(yaml::read_yaml('../../misc_data/fourth_round_idp_preprocessing/output/t1_covar_for_pc.yaml')) %>% 
  #     mutate(t1_or_dmri = 'T1')
  # )
  # rbind(
  #   idp_annot,
  #   pc_annot[, colnames(idp_annot)]
  # )
  idp_annot
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

remove_probtrack_idp = function(df) {
  idp_annot = load_idp_annot()
  tbss_idps = idp_annot$IDP[ substr(idp_annot$subtype, 1, 2) == 'w-' ]
  df = df[ !df$IDP %in% tbss_idps, ]
  df
}

filter_pred_perf = function(df) {
  df = df %>% mutate(idp_id = paste(idp_type, IDP, model))
  idp_sig = read.table('supp_table_2.tsv', header = T, sep = '\t')
  idp_sig = idp_sig %>% filter(is_kept) %>% mutate(idp_id = paste(IDP_type, IDP, model_name))
  df = df[df$idp_id %in% idp_sig$idp_id, ]
  df %>% select(-idp_id)
}
