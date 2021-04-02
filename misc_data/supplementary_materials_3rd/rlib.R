get_pcs_for_3rd = function() {
  t1 = get_pcs_from_file('../../submission/genetic_cor/third_round_t1.txt')
  dmri = get_pcs_from_file('../../submission/genetic_cor/third_round_dmri.txt')
  data.frame(IDP = c(t1, dmri), idp_type = c(rep('T1', length(t1)), rep('dMRI', length(dmri))))
}

get_pcs_from_file = function(fn) {
  ff = read.table(fn)$V1
  t1 = stringr::str_remove(basename(ff), '\\.parquet')
  t1 = t1[substr(t1, 1, 2) == 'PC']
}