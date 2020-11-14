# setwd('misc_data/process_t1/')

library(ggplot2)
library(dplyr)
library(patchwork)
options(stringsAsFactors = F)

theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/489d1519dd45246caf4756d7722bfa25/raw/90c572c0a287f2cada53811c7cd51ec14fade488/fast_linear_regression')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('https://gist.githubusercontent.com/liangyy/4c647634fe00b3f042ebd1599dda65c7/raw/9977562b65d0fb63a693fa7fa60035a37641ad2f/qqplot_by_group')


idp_meta = readRDS('../download_some_matching_files/cleanup_annot_our_idps.rds')
t1_meta = idp_meta %>% filter(modality == 'T1')
df_idp = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_final-phenotypes.parquet')
df_t1 = df_idp[ , c('individual', paste0('IDP-', t1_meta$FieldID))]

t1_meta_rds = 't1_meta.rds'
if(!file.exists(t1_meta_rds)) {
  # head_size = idp_meta %>% filter(FieldID == 25000)
  # df_covar = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_covariates.parquet')
  # df_covar = inner_join(df_idp %>% select(individual), df_covar %>% mutate(individual = as.character(individual)), by = 'individual')
  # df_others = data.table::fread('../../../ptrs-ukb/output/query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
  # df_others = left_join(df_idp %>% select(individual), df_others %>% select(eid, height) %>% mutate(individual = as.character(eid)), by = 'individual')
  # df_covar = df_covar[, paste0('IDP-', head_size$FieldID, '_out')]
  # cor(df_t1$`IDP-25002`, df_covar$`IDP-25000_out`)
  # plot(df_others$height, 1 / df_covar$`IDP-25000_out`)
  # plot(df_others$height, df_t1$`IDP-25010`)
  # head(t1_meta %>% select(Field, FieldID), 12)
  detect_str = function(ss, str) {
    !is.na(stringr::str_match(ss, str)[, 1])
  }
  remove_str = function(ss, str) {
    stringr::str_remove(ss, str)
  }
  extract_str = function(ss, str) {
    stringr::str_match(ss, str)[, 2]
  }
  t1_str_parser = function(ss) {
    # about norm head size
    norm_head_str = ' \\(normalised for head size\\)$'
    is_norm_head_size = detect_str(ss, norm_head_str)
    ss = remove_str(ss, norm_head_str)
    # gray matter
    matter_str = '(grey matter|white matter|grey\\+white matter)'
    matter_type = extract_str(ss, matter_str)
    ss = remove_str(ss, matter_str)
    # is left|right|vermis
    lr_str = '\\((left|vermis|right)\\)'
    lr = extract_str(ss, lr_str)
    ss = remove_str(ss, lr_str)
    # remove "Volumne of" and "in"
    ss = remove_str(ss, 'Volume of')
    ss = remove_str(ss, ' in ')
    # remove extra space in the prefix and suffix
    ss = remove_str(ss, '^ +')
    pos = remove_str(ss, ' +$')
    pos[ pos == '' ] = 'overall'
    pos[ pos == 'brain,' ] = 'overall'
    data.frame(position = tolower(pos), lr = lr, matter_type = matter_type, normalized_by_head_size = is_norm_head_size)
  }
  ss = t1_str_parser(t1_meta$Field)
  t1_meta = cbind(t1_meta %>% select(Field, FieldID, id_in_their_paper, Path, Participants, Notes, Link, modality), ss)
  
  saveRDS(t1_meta, t1_meta_rds)
} else {
  t1_meta = readRDS(t1_meta_rds)
}

corr = cor(df_t1[, -1])

# what we want to do:
# 1. Perform PCA
# 2. Look at the how to interpret PC
# 3. Look at the residual correlation
# 4. Look at the QQ-plot with simulated one signal
source('../../rmd/rlib_calc.R')

# SVD
mat_std = standardize(df_t1[, -1])
res = svd(mat_std)
pve = res$d^2 / sum(res$d^2)
plot(cumsum(pve))
plot(res$v[, 1], res$v[, 2])

vmat = res$v
vmat = as.data.frame(vmat) 
vmat$idp = colnames(df_t1)[-1]
df_v = vmat # %>% reshape2::melt(id.var = 'idp')
df_v = df_v %>% left_join(t1_meta %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, normalized_by_head_size, position, matter_type, lr), by = 'idp')
df_v %>% ggplot() + geom_point(aes(x = V1, y = V2, color = matter_type))

# correlate PC with annotation
out = list()
library(lme4)
lmcont = lmerControl(
  check.nobs.vs.nlev = "ignore",
  check.nobs.vs.rankZ = "ignore",
  check.nlev.gtreq.5 = "ignore",
  check.nobs.vs.nRE="ignore",
  check.rankX =    c("ignore"),
  check.scaleX = "ignore",
  check.formula.LHS="ignore",
  check.conv.grad = .makeCC("warning", tol = 1e-3, relTol = NULL),
  check.conv.singular = 'ignore'
)
for(i in 1 : (ncol(vmat) - 1)) { 
  if(i %% 10 == 0) {
    message(i)
  }
  for(v in c('normalized_by_head_size', 'position', 'matter_type', 'lr')) {
    tmp = t1_meta[, v]
    tmp[is.na(tmp)] = 'NA'
    pval = corr_cat_cont(tmp, df_v[, i], control = lmcont)
    out[[length(out) + 1]] = data.frame(pc = i, feature = v, pval = pval)
  }
}
out = do.call(rbind, out)
out$pval = pmax(out$pval, 1e-20)
p1 = out %>% ggplot() + geom_point(aes(x = pc, y = -log10(pval))) + facet_wrap(~feature) 

plist = list()
# top factor per feature
for(v in c('normalized_by_head_size', 'position', 'matter_type', 'lr')) {
  tmp = out %>% filter(feature == v)
  best = tmp[ which.min(tmp$pval), ]$pc
  tmp = df_v[, c(paste0('V', best), v)]
  colnames(tmp) = c('pc', 'feature')
  tmp$feature[is.na(tmp$feature)] = 'NA'
  plist[[length(plist) + 1]] = tmp %>% ggplot() + geom_boxplot(aes(x = feature, y = pc)) + ggtitle(paste0(v, ': ', 'PC', best))
}
p2 = do.call(gridExtra::grid.arrange, plist)
p1 + p2
ggsave('explain_pc.png', p1 + p2, width = 16, height = 8)

# residual correlation
t1_meta$matter_type[is.na(t1_meta$matter_type)] = 'NA'
df_corr = corr %>% reshape2::melt() %>% 
  inner_join(t1_meta %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var1' = 'idp')) %>%
  inner_join(t1_meta %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))

myplot_magic = function(df) {
  mytmp = df # 
  # mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% mutate(order = order(Var1)) %>% group_by(matter_type.1) %>% summarize(pos = mean(order))
  mytmp %>% ggplot() +
    theme_bw() + 
    geom_raster(aes(x = as.character(Var1), y = as.character(Var2), fill = value)) + scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    # annotate(geom = "text", x = mytmp_g$pos, y = -20, label = mytmp_g$matter_type.1, size = 4) + 
    coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
    theme(
      plot.margin = unit(c(0.1, 0.1, 2, 0.1), "cm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) 
}
p0 = myplot_magic(df_corr)

pc_mat = res$u
plist2 = list()
ll = c(1, 2, 5, 10, 20, 50)
for(npc in ll) {
  message('Working on nPC = ', npc)
  tt = pc_mat[, 1 : npc]
  mat_regress = mat_std - tt %*% (t(tt) %*% mat_std)
  mat_regress_corr = cor(mat_regress)
  df_corr_tmp = mat_regress_corr %>% reshape2::melt() %>% 
    inner_join(t1_meta %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var1' = 'idp')) %>%
    inner_join(t1_meta %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, normalized_by_head_size, position, matter_type, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
  
  pp = myplot_magic(df_corr_tmp)
  pp = pp + ggtitle(paste0('nPC = ', npc))
  plist2[[length(plist2) + 1]] = pp
}

# do.call(gridExtra::grid.arrange, c(plist2, ncol = 2))
for(npc in ll) {
  ggsave(paste0('regress_out_corr_npc_', npc, '.png'), plist2[[which(ll == npc)]])
}
ggsave('regress_out_corr.png', p0)

# QQ-plot on the simulated y
mat = df_t1[, -1]
## first, on the observed data
plist = c()
clist = c()
for(n in 1:10) {
  select_idx = sample(ncol(mat), 1); y = rnorm(nrow(mat)) * 20 + rowSums(mat[, select_idx, drop = F])
  message(select_idx)
  res = fast_linear_regression(as.numeric(y), mat, covariate = matrix(1, nrow = nrow(mat), ncol = 1))
  plist = c(plist, res$pval)
  clist = c(clist, rep(paste0('repeat', n), ncol(mat)))
}

p = qqplot_by_group(plist, clist) + th + theme(legend.position = 'none'); p

## second, on the residual 
cor_mat2 = t(mat_std) %*% mat_std / nrow(mat_std)
eig_res = eigen(cor_mat2)
p = list()
for(top_n in c(0, 1, 2, 5, 10, 20, 40, 50)) {
  if(top_n == 0) {
    pve = 0
    mat_res2 = mat_std
  } else{
    pve = sum(eig_res$values[1:top_n]) / sum(eig_res$values)
    top_n_pc = mat_std %*% eig_res$vectors[, 1 : top_n]
    tmp = qr(top_n_pc)
    Q_ = qr.Q(tmp)
    
    # mat_res = mat - Q_ %*% (t(Q_) %*% mat)
    mat_res2 = mat_std - Q_ %*% (t(Q_) %*% mat_std)
  }
  
  
  plist = c()
  clist = c()
  for(n in 1:10) {
    select_idx = sample(ncol(mat_res2), 1); y = rnorm(nrow(mat_res2)) * 20 + rowSums(mat_res2[, select_idx, drop = F])
    # message(select_idx)
    res = fast_linear_regression(as.numeric(y), mat_res2, covariate = matrix(1, nrow = nrow(mat_res2), ncol = 1))
    plist = c(plist, res$pval)
    clist = c(clist, rep(paste0('repeat', n), ncol(mat_res2)))
  }
  
  p[[length(p) + 1]] = qqplot_by_group(plist, clist) + th + theme(legend.position = 'none') + ggtitle(paste0('nPC = ', top_n, ' PVE = ', signif(pve, digits = 3)))
  
  
}
out = do.call(grid.arrange, c(p, ncol = 4))
ggsave('export_mat_fac.png', out, width = 15, height = 8)

out = do.call(gridExtra::grid.arrange, c(p, ncol = 4))
ggsave('export_mat_fac.png', out, width = 15, height = 8)

