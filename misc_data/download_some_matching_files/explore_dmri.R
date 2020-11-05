# This is script to explore dMRI IDPs
# The basic goal is to have a more informative annotation (anatomy, measurement)
setwd('~/Documents/repo/github/ukb_idp_genetic_arch/misc_data/download_some_matching_files/')
options(stringsAsFactors = T)
library(ggplot2)
library(patchwork)
library(dplyr)

dmri_meta = 'annot_dmri_idps.rds'
if(file.exists(dmri_meta)) {
  df_dmri = readRDS(dmri_meta)
} else {
  df_idp = readRDS('cleanup_annot_our_idps.rds')
  df_dmri = df_idp[ df_idp$modality == 'dMRI', ]
  
  measurement_strs = c('ISOVF', 'OD', 'ICVF', 'L3', 'L2', 'L1', 'MO', 'MD', 'FA')
  modes = c('Weighted-mean TYPE in tract', 'Mean TYPE in')
  
  string_a_in_b = function(a, b) {
    !is.na(stringr::str_extract(b, a))
  }
  
  parse_weighted_mean = function(str) {
    parse_(str, '^Weighted-mean ', ' in tract')
  }
  
  parse_mean = function(str) {
    parse_(str, '^Mean ', ' in')
  }
  
  parse_ = function(str, pre, suf) {
    pat = paste0(pre, '([A-Za-z0-9]+)', suf)
    mea = stringr::str_match(str, pat)[1, 2]
    pat2 = paste0(pre, mea, suf, ' ')
    remain = stringr::str_remove(str, pat2)
    direction = NA
    pos = remain
    for(dd in c('right', 'left')) {
      dd_str = paste0('\\(', dd, '\\)')
      if(string_a_in_b(dd_str, remain)) {
        pos = stringr::str_remove(remain, paste0(' ', dd_str))
        direction = dd
        break
      }
    }
    list(lr = direction, pos = pos, measure = mea)
  }
  
  parse_dmri_field = function(x) {
    if(string_a_in_b(a = '^Weighted-mean', b = x)) {
      res = parse_weighted_mean(x)
      type = 'weighted_mean_in_tract'
    } else if(string_a_in_b(a = '^Mean', b = x)) {
      res = parse_mean(x)
      type = 'mean'
    } else {
      message('Cannot work on ', x)
    }
    if(!is.null(res$lr)) {
      lr = res$lr
    } else {
      lr = NA
    } 
    pos = res$pos
    mea = res$measure
    data.frame(type = type, position = position_cleaner(pos), measure = mea, lr = lr)
  }
  
  position_cleaner = function(ss) {
    stringr::str_remove(ss, ' on FA skeleton')
  }
  ss = sapply(df_dmri$Field, parse_dmri_field, simplify = F)
  ss = do.call(rbind, ss)
  df_dmri = cbind(df_dmri, ss)
  
  saveRDS(df_dmri, dmri_meta)
}

# load idp
df_idp = arrow::read_parquet('~/Desktop/tmp/ukb_idp/idp_phenotypes/2020-05-18_final-phenotypes.parquet')
colnames(df_idp)[-1] = unlist(lapply(strsplit(as.character(colnames(df_idp)[-1]), '-'), function(x){x[[2]]}))
dmri = df_idp[, c('individual', df_dmri$FieldID)] %>% reshape2::melt(id.var = 'individual') %>% rename(idp = variable)
dmri = dmri %>% left_join(df_dmri %>% mutate(idp = as.character(FieldID)) %>% select(idp, type, position, measure, lr), by = 'idp')
tmp = dmri %>% filter(!is.na(lr)) %>% reshape2::dcast(individual + type + position + measure ~ lr, value.var = 'value')
tmp = tmp %>% filter(!is.na(left) & !is.na(right))

p1 = tmp %>% sample_n(10000) %>% ggplot() + geom_point(aes(x = right, y = left), alpha = 0.2) + facet_wrap( ~ measure)
p2 = tmp %>% filter(measure == 'L1') %>% sample_n(10000) %>% ggplot() + geom_point(aes(x = right, y = left), alpha = 0.2) + facet_wrap( ~ position)
p1 + p2

mat = df_idp[, c('individual', df_dmri$FieldID)]
mat = as.matrix(mat[, -1])
mat_clean = apply(mat, 2, function(x) { (x - mean(x)) / sd(x)})
res = svd(mat_clean)
pve = res$d^2 / sum(res$d^2)
plot(cumsum(pve)[1:100])
plot(res$v[, 1], res$v[, 2])
vmat = res$v
vmat = as.data.frame(vmat) 
vmat$idp = colnames(mat)
df_v = vmat # %>% reshape2::melt(id.var = 'idp')
df_v = df_v %>% left_join(df_dmri %>% mutate(idp = as.character(FieldID)) %>% select(idp, type, position, measure, lr), by = 'idp')
df_v %>% ggplot() + geom_point(aes(x = V1, y = V2, color = measure))

corr_cat_cont = function(categ, cont, ...) {
  mod = lmer(cont ~ (1 | categ), ...)
  if(!isSingular(mod, tol = 1e-05)) {
    test = lmerTest::ranova(mod)
    return(test$`Pr(>Chisq)`[2])
  } else {
    return(1)
  }
}

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
  for(v in c('type', 'position', 'measure', 'lr')) {
    tmp = df_v[, v]
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
for(v in c('type', 'position', 'measure', 'lr')) {
  tmp = out %>% filter(feature == v)
  best = tmp[ which.min(tmp$pval), ]$pc
  tmp = df_v[, c(paste0('V', best), v)]
  colnames(tmp) = c('pc', 'feature')
  tmp$feature[is.na(tmp$feature)] = 'NA'
  plist[[length(plist) + 1]] = tmp %>% ggplot() + geom_boxplot(aes(x = feature, y = pc)) + ggtitle(paste0(v, ': ', 'PC', best))
}
p2 = do.call(gridExtra::grid.arrange, plist)
p1 + p2

# run flash
library(flashr)
flash_cache = 'explore_dmri_flash.rds'
if(file.exists(flash_cache)) {
  f = readRDS(flash_cache)
} else {
  data = flash_set_data(mat_clean[1:2000, ])
  f = flash_add_greedy(
    data,
    Kmax = 100,
    var_type = 'by_column',
    init_fn = 'udv_si',
    ebnm_fn = 'ebnm_pn'
  )
  f = flash_backfit(data, f)
  saveRDS(f, flash_cache)
}

fres = flash_get_ldf(f)
f_vmat = fres$f %>% as.data.frame
f_vmat$idp = colnames(mat)
f_vmat = f_vmat %>% left_join(df_dmri %>% mutate(idp = as.character(FieldID)) %>% select(idp, type, position, measure, lr), by = 'idp')
fpve = f$pve
plot(cumsum(fpve))
plot(log(pve[1:100]), log(fpve))

f_vmat %>% ggplot() + geom_point(aes(x = V1, y = V2, color = measure))
out2 = list()
# loop over all factors
for(i in 1 : 100) { 
  if(i %% 10 == 0) {
    message(i)
  }
  for(v in c('type', 'position', 'measure', 'lr')) {
    tmp = f_vmat[, v]
    tmp[is.na(tmp)] = 'NA'
    pval = corr_cat_cont(tmp, f_vmat[, i], control = lmcont)
    out2[[length(out2) + 1]] = data.frame(pc = i, feature = v, pval = pval)
  }
}
out2 = do.call(rbind, out2)
out2$pval = pmax(out2$pval, 1e-20)

p1 = out2 %>% ggplot() + geom_point(aes(x = pc, y = -log10(pval))) + facet_wrap(~feature) 

plist = list()
# top factor per feature
for(v in c('type', 'position', 'measure', 'lr')) {
  tmp = out2 %>% filter(feature == v)
  best = tmp[ which.min(tmp$pval), ]$pc
  tmp = df_v[, c(paste0('V', best), v)]
  colnames(tmp) = c('pc', 'feature')
  tmp$feature[is.na(tmp$feature)] = 'NA'
  plist[[length(plist) + 1]] = tmp %>% ggplot() + geom_boxplot(aes(x = feature, y = pc)) + ggtitle(paste0(v, ': ', 'PC', best))
}
p2 = do.call(gridExtra::grid.arrange, plist)
p1 + p2


corr = cor(mat_clean)
df_corr = corr %>% reshape2::melt() %>% 
  inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var1' = 'idp')) %>%
  inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
df_corr = df_corr %>% mutate(
  label1 = paste0(Var1, '\n', measure.1),
  label2 = paste0(Var2, '\n', measure.2)
)
myplot_magic = function(df) {
  mytmp = df # 
  mytmp_g = mytmp[!duplicated(mytmp$Var1), ] %>% mutate(order = order(Var1)) %>% group_by(measure.1) %>% summarize(pos = mean(order))
  mytmp %>% ggplot() +
    theme_bw() + 
    geom_raster(aes(x = as.character(Var1), y = as.character(Var2), fill = value)) + scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
    annotate(geom = "text", x = mytmp_g$pos, y = -20, label = mytmp_g$measure.1, size = 4) + coord_equal(expand = FALSE, clip = "off", ylim = c(0, sum(!duplicated(mytmp$Var1)))) +
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
p0 = myplot_magic(df_corr %>% filter(type.1 == 'mean', type.2 == 'mean') )

pc_mat = res$u
plist2 = list()
ll = c(1, 2, 5, 10, 20, 50)
for(npc in ll) {
  message('Working on nPC = ', npc)
  tt = pc_mat[, 1 : npc]
  mat_regress = mat_clean - tt %*% (t(tt) %*% mat_clean)
  mat_regress_corr = cor(mat_regress)
  df_corr_tmp = mat_regress_corr %>% reshape2::melt() %>% 
    inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var1' = 'idp')) %>%
    inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
  pp = myplot_magic(df_corr_tmp %>% filter(type.1 == 'mean', type.2 == 'mean'))
  pp = pp + ggtitle(paste0('nPC = ', npc))
  plist2[[length(plist2) + 1]] = pp
}

collect_corr = list()
for(tt in unique(df_corr_tmp$position.1)) {
  kk = df_corr_tmp %>% filter(position.1 == tt, type.1 == 'mean', position.2 == tt, type.2 == 'mean', (is.na(lr.1) & is.na(lr.2)) | (lr.1 == lr.2), Var1 != Var2)
  collect_corr[[length(collect_corr) + 1]] = kk %>% mutate(type_label = 'same region and LR, diff measurement')
  kk = df_corr_tmp %>% filter(position.1 == tt, type.1 == 'mean', position.2 == tt, type.2 == 'mean', (!is.na(lr.1) & !is.na(lr.2)) & (lr.1 != lr.2), Var1 != Var2)
  collect_corr[[length(collect_corr) + 1]] = kk %>% mutate(type_label = 'same region diff LR and measurement')
}
kk = df_corr_tmp %>% filter(position.1 != position.2, type.1 == 'mean', type.2 == 'mean')
collect_corr[[length(collect_corr) + 1]] = kk %>% mutate(type_label = 'anything else')
collect_corr = do.call(rbind, collect_corr)
collect_corr %>% group_by(type_label) %>% summarize(median_abs_corr = median(abs(value)))

# do.call(gridExtra::grid.arrange, c(plist2, ncol = 2))
for(npc in ll) {
  ggsave(paste0('regress_out_corr_npc_', npc, '.png'), plist2[[which(ll == npc)]])
}
ggsave('regress_out_corr.png', p0)


flash_mat = fres$l
plist3 = list()
ll = c(1, 2, 5, 10, 20, 50)
for(nf in ll) {
  message('Working on nFactor = ', nf)
  tt = flash_mat[, 1 : nf]
  tmp = qr(tt)
  qq = qr.Q(tmp)
  mat_regress = mat_clean[1:2000, ] - qq %*% (t(qq) %*% mat_clean[1:2000, ])
  mat_regress_corr = cor(mat_regress)
  df_corr_tmp = mat_regress_corr %>% reshape2::melt() %>% 
    inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var1' = 'idp')) %>%
    inner_join(df_dmri %>% mutate(idp = FieldID) %>% select(idp, type, position, measure, lr), by = c('Var2' = 'idp'), suffix = c('.1', '.2'))
  pp = myplot_magic(df_corr_tmp %>% filter(type.1 == 'mean', type.2 == 'mean'))
  pp = pp + ggtitle(paste0('nFlash = ', nf))
  plist3[[length(plist3) + 1]] = pp
}
# do.call(gridExtra::grid.arrange, c(plist2, ncol = 2))
for(nf in ll) {
  ggsave(paste0('regress_out_corr_nf_', nf, '.png'), plist3[[which(ll == nf)]])
}
# ggsave('regress_out_corr.png', p0)


# take the residual after regressing out the top 50 PCs
# and do a region-by-region PCA analysis
region_list = unique(df_dmri$position[ df_dmri$type == 'mean' ])
for(rr in region_list) {
  df_dmri_sub = df_dmri[ df_dmri$position == rr & df_dmri$type == 'mean',  ]
  mat_regress_sub = mat_regress[, as.character(df_dmri_sub$FieldID) ]
  mat_regress_sub = apply(mat_regress_sub, 2, function(x) {( x - mean(x) ) / sd(x)})
  res_sub = svd(mat_regress_sub)
  vmat_sub = res_sub$v
  vmat_sub = as.data.frame(vmat_sub) 
  vmat_sub$idp = colnames(mat_regress_sub)
  df_v_sub = vmat_sub # %>% reshape2::melt(id.var = 'idp')
  df_v_sub = df_v_sub %>% left_join(df_dmri_sub %>% mutate(idp = as.character(FieldID)) %>% select(idp, type, position, measure, lr), by = 'idp')
  if(is.na(df_v_sub$lr[1])) {
    ppp = df_v_sub %>% ggplot() + geom_point(aes(x = V1, y = V2, color = measure))
  } else {
    ppp = df_v_sub %>% ggplot() + geom_point(aes(x = V1, y = V2, color = measure, shape = lr))
  }
  ggsave(paste0(stringr::str_replace_all(rr, ' ', '_'), '.png'), ppp + ggtitle(rr))
}
