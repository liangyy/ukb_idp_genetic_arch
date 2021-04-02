# setwd('misc_data/supplementary_materials_3rd/')

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

options(stringsAsFactors = F)

foldern = 'polygenicity'
dir.create(foldern)

color_code2 = c('T1' = 'orange', 'dMRI' = 'blue')

# load h2
{ 
df1 = read.table('~/Desktop/tmp/ukb_idp/heritability_3rd_round/third_round_dmri.tsv.gz', header = T)
df2 = read.table('~/Desktop/tmp/ukb_idp/heritability_3rd_round/third_round_t1.tsv.gz', header = T)
df1$is_pc = substr(df1$phenotype, 1, 2) == 'PC'
df2$is_pc = substr(df2$phenotype, 1, 2) == 'PC'
df1$pc = rep('Other brain IDP', nrow(df1))
df1$pc[df1$is_pc] = 'IDP PC'
df2$pc = rep('Other brain IDP', nrow(df2))
df2$pc[df2$is_pc] = 'IDP PC'
df_h2 = rbind(
  df1 %>% mutate(idp_type = factor('dMRI', levels = c('T1', 'dMRI'))), 
  df2 %>% mutate(idp_type = factor('T1', levels = c('T1', 'dMRI')))
)
df_h2 = df_h2 %>% rename(IDP = phenotype)
}

pcs = df_h2 %>% filter(is_pc) %>% select(IDP, idp_type) 

# load polygenicity
{
  idps = read.delim2('supp_table_1.tsv', header = T) %>% mutate(IDP = paste0('IDP-', ukb_field), idp_type = t1_or_dmri) %>% select(IDP, idp_type)
  type_list = list(dMRI = 'third_round_dmri', T1 = 'third_round_t1')
  idps = rbind(
    idps, 
    pcs
  )
  dd2 = list()
  for(i in 1 : nrow(idps)) {
    fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/', type_list[[idps$idp_type[i]]], '/', idps$IDP[i], '.sld4m_all.csv')
    if(file.exists(fn)) {
      tmp = read.csv(fn)
      dd2[[length(dd2) + 1]] = tmp %>% mutate(IDP = idps$IDP[i], idp_type = idps$idp_type[i])
    }
  }
  dd2 = do.call(rbind, dd2)
  df_poly = dd2 %>% filter(Var1 == 'Manual_aggregated')
  df_poly = df_poly %>% mutate(stable = Ma_est > 0 & Ma_est > 1.96 * Ma_err)
}

# load prediction performance
{
  df_perf = read.delim2('supp_table_2.tsv')
  df_perf$Spearman = as.numeric(df_perf$Spearman)
  df_perf = df_perf %>% rename(idp_type = IDP_type)
  tmp = inner_join(
    df_perf %>% filter(model_name == 'ridge'), 
    df_perf %>% filter(model_name == 'elastic net'), 
    by = c('IDP', 'idp_type'),
    suffix = c('.ridge', '.en')
  )
  df_perf = tmp %>% mutate(
    Spearman.diff = Spearman.ridge - Spearman.en
  )
}

# load external GWASs
{
  ind = read.table('../ld4m/ld4m_indep_gwas.txt')$V1
  dd3 = list()
  for(i in ind) {
    fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/indep_gwas/', i, '.sld4m_all.csv')
    if(!file.exists(fn)) {
      next
    }
    tmp = read.csv(fn)
    dd3[[length(dd3) + 1]] = tmp %>% mutate(phenotype = i, type = 'independent')
  }
  ext = read.table('../ld4m/ld4m_external_traits.txt')$V1
  for(i in ext) {
    fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/external_gwas/', i, '.sld4m_all.csv')
    if(!file.exists(fn)) {
      next
    }
    tmp = read.csv(fn)
    dd3[[length(dd3) + 1]] = tmp %>% mutate(phenotype = i, type = 'external')
  }
  dd3 = do.call(rbind, dd3)
  
  ss = read.table('../ld4m/sample_size.txt', header = T)
  dd3_sub = dd3 %>% filter(Var1 == 'Manual_aggregated')
  
  df_poly_other = left_join(dd3_sub, ss, by = c('phenotype' = 'trait')) %>% mutate(stable = Ma_est > 0 & Ma_est > 1.96 * Ma_err)
  
  # load paper results
  parse_number = function(xx) {
    tmp = stringr::str_match(xx, '([0-9\\.]+) \\(([0-9\\.]+)\\)')
    return(list(est = as.numeric(tmp[, 2]), se = as.numeric(tmp[, 3])))
  }
  dmap = readRDS('../ld4m/ld4m_map.rds')
  oo = parse_number(dmap$log10MeCommon)
  dmap$Ma_est = 10 ^ oo$est
  dmap$Ma_err = (10 ^ (oo$est + oo$se) - 10 ^ (oo$est - oo$se)) / 2
  df_poly_other = df_poly_other %>% left_join(dmap %>% select(Ma_est, Ma_err, FILE, Neff, Trait), by = c('phenotype' = 'FILE'), suffix = c('.myrun', '.paper'))
}

plot_myrun_vs_paper = T
if(isTRUE(plot_myrun_vs_paper)) {
  pp = df_poly_other %>% filter(stable, !is.na(Ma_est.paper), Ma_est.paper > Ma_err.paper * 1.96) %>% 
    ggplot() + 
    geom_point(aes(x = Ma_est.paper, y = Ma_est.myrun)) + 
    geom_errorbar(aes(x = Ma_est.paper, ymax = Ma_est.myrun + 1.96 * Ma_err.myrun, ymin = Ma_est.myrun - 1.96 * Ma_err.myrun), width = 0.025) + 
    geom_errorbarh(aes(y = Ma_est.myrun, xmax = Ma_est.paper + 1.96 * Ma_err.paper, xmin = Ma_est.paper - 1.96 * Ma_err.paper), height = 0.025) +
    scale_x_log10() + scale_y_log10() + coord_equal() + th + geom_abline(slope = 1, intercept = 0, color = 'blue') + 
    # ggtitle('Comparing Me reported in \n the paper and my runs') +
    ggrepel::geom_label_repel(
      aes(x = Ma_est.paper, y = Ma_est.myrun, label = Trait)
    ) + 
    xlab('Me reported in \n O\'Connor et al 2019') + 
    ylab('Me estimated from \n our pipeline')
  ggsave(paste0(foldern, '/our_runs_vs_paper.png'), pp, width = 6, height = 6)
}

plot_idp_vs_some_traits = T
if(isTRUE(plot_idp_vs_some_traits)) {
  selected_traits = c('Sunburn', 'RBC count', 'Type II diabetes', 'Platelet count', 'Eczema', 'Height', 'BMI', 'Neuroticism', 'Morning person', 'College', 'Smoking status', 'Schizophrenia')
  
  df_poly %>% filter(substr(IDP, 1, 2) != 'PC') %>% group_by(idp_type) %>% summarize(
    number_of_stable_Me = sum(stable, na.rm = T),
    total_number = n()
  ) %>% mutate(
    fraction_of_stable_Me = number_of_stable_Me / total_number
  ) %>% pander::pander(caption = 'Non PC IDPs')
  
  df_poly %>% group_by(idp_type) %>% summarize(
    number_of_stable_Me = sum(stable, na.rm = T),
    total_number = n(),
    median_Me = median(Ma_est, na.rm = T)
  ) %>% pander::pander(caption = 'Median Me')
  
  df_poly %>% summarize(
    number_of_stable_Me = sum(stable, na.rm = T),
    total_number = n(),
    median_Me = median(Ma_est, na.rm = T)
  ) %>% pander::pander(caption = 'Median Me (not stratified)')
  
  df_poly %>% filter(substr(IDP, 1, 2) == 'PC') %>% group_by(idp_type) %>% summarize(
    number_of_stable_Me = sum(stable, na.rm = T),
    total_number = n()
  ) %>% mutate(
    fraction_of_stable_Me = number_of_stable_Me / total_number
  ) %>% pander::pander(caption = 'PC IDPs')
  
  tmp2 = df_poly %>% filter(stable) %>% mutate(phenotype = paste(idp_type, IDP)) %>% 
    mutate(pheno = factor(phenotype, levels = phenotype[order(idp_type, Ma_est)]))
  labx = 550
  pp = tmp2 %>% ggplot() + 
    geom_errorbar(aes(x = pheno, ymax = Ma_est + 1.96 * Ma_err, ymin = Ma_est - 1.96 * Ma_err), width = 0.025, color = 'gray') + 
    geom_point(aes(x = pheno, y = Ma_est, color = idp_type))  + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    theme(legend.position = c(0.2, 0.95), legend.title = element_blank()) +
    guides(color = guide_legend(ncol = 2)) +
    th + 
    geom_segment(
      data = df_poly_other %>% filter(!is.na(Trait), stable, Trait %in% selected_traits), 
      aes(y = Ma_est.myrun, yend = Ma_est.myrun, x = 0, xend = labx),
      linetype = 2,
      alpha = 0.5
    ) + 
    ggrepel::geom_label_repel(
      data = df_poly_other %>% filter(!is.na(Trait), stable, Trait %in% selected_traits), 
      aes(x = labx, y = Ma_est.myrun, label = Trait),
      direction = "y",
      hjust = 0,
      nudge_x = 50,
      segment.size = 0.2,
      force = .5,
    ) +
    # geom_point(
      # data = df_poly_other %>% filter(!is.na(Trait), stable, Trait %in% selected_traits), 
      # aes(x = labx, y = Ma_est.myrun),
    # ) +
    scale_y_log10() + 
    coord_cartesian(ylim = c(100, 50000), xlim = c(0, 700)) +
    xlab('Brain IDPs') +
    ylab('Estimated Me') +
    scale_color_manual(values = color_code2) 
  ggsave(paste0(foldern, '/overview_of_brain_idp_Me.png'), pp, width = 8, height = 5)
}

plot_h2_vs_me = T
if(isTRUE(plot_h2_vs_me)) {
  tmp = inner_join(df_h2, df_poly, by = c('IDP', 'idp_type'))
  tmp = inner_join(tmp, df_perf, by = c('IDP', 'idp_type'))
  tmp %>% filter(stable) %>% 
    ggplot() + geom_point(aes(x = h2, y = Ma_est))
  pp = tmp %>% filter(stable) %>% 
    ggplot() + geom_point(aes(x = Ma_est, y = Spearman.diff, color = pc)) + 
    scale_x_log10() +
    xlab('Estimated Me') +
    ylab('Difference in Spearman correlation \n between ridge and elastic net') +
    th +
    theme(legend.position = c(0.8, 0.18), legend.title = element_blank()) +
    scale_color_manual(values = c('IDP PC' = 'red', 'Other brain IDP' = 'black'))
  ggsave(paste0(foldern, '/pred_perf_vs_Me.png'), pp, width = 5.5, height = 5)
}
