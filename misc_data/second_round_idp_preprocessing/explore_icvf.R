# setwd('misc_data/second_round_idp_preprocessing/')

library(dplyr)
library(ggplot2)
library(patchwork)
source('scripts/rlib.R')

force_run = T
outdir = 'explore_icvf'
dir.create(outdir)

# fig size
ww = 12
hh = 10
# END

tag = 'ICVF'

outputs = c(paste0(outdir, '/', tag, '_probt_tbss.png'), paste0(outdir, '/', tag, '_probt.png'), paste0(outdir, '/', tag, '_tbss.png'))
if(sum(file.exists(outputs)) != 3) {
  doit = T
}

get_upper_tri = function(mat) {
  return(mat[upper.tri(mat)])
}
plot_image = function(mat) {
  p = mat %>% reshape2::melt() %>% ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradientn(limits = c(-1, 1), colors = c("blue", "white", "red")) + 
    theme(axis.text.x = element_blank())
  p
}
do_all = function(mat) {
  tmp1 = cor(mat)
  p1 = plot_image(tmp1)
  res = split_by_pca(mat, pve_cutoff = 0.1, skip = F)
  tmp2 = cor(res$residual)
  p2 = plot_image(tmp2)
  p4 = data.frame(pc_loading = res$pc_loadings, idp = colnames(mat)) %>% 
    ggplot() + geom_point(aes(x = pc_loading, y = idp))
  p3 = data.frame(before = get_upper_tri(tmp1), after = get_upper_tri(tmp2)) %>% 
    ggplot() + geom_histogram(aes(x = before, fill = 'before'), alpha = 0.5, binwidth = 0.1) +
    geom_histogram(aes(x = after, fill = 'after'), alpha = 0.5, binwidth = 0.1)
  list(ps = list(p1, p2, p3, p4), res = res)
}
plot_all = function(ps, title) {
  p1 = ps[[1]] + ggtitle('Before PC adj')
  p2 = ps[[2]] + ggtitle('After PC adj')
  p3 = ps[[3]] + ggtitle('Histogram of corr')
  p4 = ps[[4]] + ggtitle('PC loading')
  ((p1 + p2) / (p3 + p4)) + plot_annotation(title = title)
}

if(doit | force_run) {
  # ICVF as example!
  
  # extract ICVF
  idps = read.delim2('../supplementary_materials/supp_table_1.tsv')
  idps_icvf = idps %>% filter(t1_or_dmri == 'dMRI', dmri_measure == tag)
  
  # load dmri
  dmri_mat = arrow::read_parquet('output/dmri.original.all_covar.parquet')
  # head(dmri_mat)
  
  # ICVF all
  icvf_all = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf$ukb_field)])
  kk = do_all(icvf_all)
  
  
  # ICVF ProbTrack
  idps_icvf_prob = idps_icvf %>% filter(measurement_type == 'dMRI weighted means (probabilistic-tractography-based measurement)')
  icvf_prob = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_prob$ukb_field)])
  kk2 = do_all(icvf_prob)
  
  
  # ICVF TBSS
  idps_icvf_tbss = idps_icvf %>% filter(measurement_type == 'dMRI skeleton (TBSS-style measurement)')
  icvf_tbss = as.matrix(dmri_mat[, paste0('IDP-', idps_icvf_tbss$ukb_field)])
  kk3 = do_all(icvf_tbss)
  
  ggsave(outputs[1], plot_all(kk$ps, 'ICVF: TBSS + ProbTrack'), width = ww, height = hh)
  ggsave(outputs[2], plot_all(kk2$ps, 'ICVF: ProbTrack only'), width = ww, height = hh)
  ggsave(outputs[3], plot_all(kk3$ps, 'ICVF: TBSS only'), width = ww, height = hh)
  
}

