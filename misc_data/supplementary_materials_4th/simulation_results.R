# setwd('misc_data/supplementary_materials_4th/')
# NOTE: See simulation details/scripts at https://github.com/hakyimlab/yanyu-notebook/tree/master/submission/date_012923

library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
source('https://gist.githubusercontent.com/liangyy/af362f0a62aa440caae5b4423eb25bab/raw/8a042392fc359ac322af820bc8bc69d0265fa550/get_ci_from_pvalue.R')
outdir = 'simulation_results'
dir.create(outdir)

# look at results from all seeds
twas_permzs <- c()
corr_factors <- c(1.0, 1.05, 1.1)
for(ff in corr_factors) {
  twas_name <- glue::glue('twas_permz(f={sprintf("%.2f", ff)})')
  twas_permzs <- c(twas_permzs, twas_name)
}
new <- data.frame(
  method = c('indiv_twas', 'twas', 'twas_mixed', twas_permzs),
  nn = c('BrainXcan', 'S-BrainXcan', 'Mix-BrainXcan', twas_permzs)
) %>% mutate(nn = factor(nn, levels = nn))
seeds <- c(1, 10, 100, 1000, 10000)
dat.all <- list()
dat.all.alt <- list()
for(seed in seeds) {
  tmp <- readRDS(paste0('~/Desktop/cache_tmp_01292023.seed_', seed, '.rds'))
  {
    res <- tmp$res
    res_perm <- tmp$res_perm
    # some process
    res_perm_z_std <- res_perm %>% group_by(idx, gwas_n, model) %>%
      summarize(z_std = sd(z)) %>% ungroup()
    res2 <- res
    for(ff in c(1.0, 1.05, 1.1)) {
      twas_name <- glue::glue('twas_permz(f={sprintf("%.2f", ff)})')
      res.permz_adj <- res %>% filter(method == 'twas') %>%
        left_join(res_perm_z_std, by = c('idx', 'gwas_n', 'model')) %>%
        mutate(z = z / (z_std * ff)) %>%
        mutate(method = twas_name) %>%
        select(-z_std)
      res2 <- rbind(res2, res.permz_adj)
    }
    
    res22 <- res2 %>%
      left_join(new, by = 'method')
    # null
    dat.null <- res22 %>% 
      filter(method != 'marginal_true', model == 'null')  %>%
      filter(nn != 'Bacon') %>%
      mutate(pval = exp(pnorm(abs(z), log.p = T, lower.tail = F)) * 2)
    # alternative
    signal <- res22 %>% 
      filter(method == 'alt_true', model == 'alt') 
    dat.alt <- res22 %>% 
      filter(method != 'alt_true', model == 'alt') %>%
      left_join(signal %>% select(bhat, idx, midx), by = c('idx', 'midx'), suffix = c('', '.true')) %>%
      mutate(pval = exp(pnorm(abs(z), log.p = T, lower.tail = F)) * 2) 
    
    dat.all[[length(dat.all) + 1]] <- dat.null
    dat.all.alt[[length(dat.all.alt) + 1]] <- dat.alt
  }
}
dat.all <- do.call(rbind, dat.all)
dat.all.alt <- do.call(rbind, dat.all.alt)

p.min <- 1e-10
alphas <- c(0.01, 0.05, 0.1)
mtd_color <- c('Mix-BrainXcan' = 'gray', 'S-BrainXcan-Adj' = 'blue')
label_color <- c('null' = 'black', 'signal' = 'red')
methods2show <- c('BrainXcan', 'S-BrainXcan', 'Mix-BrainXcan', 'S-BrainXcan-Adj')
new2 <- new[, c('nn'), drop = FALSE]
new2$nn2 <- new2$nn
new2$nn2 <- as.character(new2$nn2)
new2$nn2[new2$nn == 'twas_permz(f=1.10)'] <- 'S-BrainXcan-Adj'
new2 <- new2 %>% filter(nn2 %in% methods2show) %>% 
  mutate(nn2 = factor(nn2, levels = methods2show))
# qq-plot null
{
  dat.null.plot <- dat.all %>%
    inner_join(new2, by = 'nn') %>%
    group_by(method, gwas_n) %>% 
    mutate(p_exp = rank(pval, ties.method = 'random') / (n() + 1)) %>%
    mutate(
      ci.high = get_ci_from_pvalue(pval)$ci.high,
      ci.low = get_ci_from_pvalue(pval)$ci.low) %>%
    ungroup() %>%
    mutate(pval = pmax(pval, p.min))
  p <- dat.null.plot %>% ggplot() +
    geom_point(aes(x = -log10(p_exp), -log10(pval))) +
    geom_line(aes(x = -log10(p_exp), -log10(ci.high)), color = 'gray') +
    geom_line(aes(x = -log10(p_exp), -log10(ci.low)), color = 'gray') +
    facet_grid(gwas_n ~ nn2, scales = 'free') +
    xlab(expression(paste(-log[10], p[expected]))) + 
    ylab(expression(paste(-log[10], p[observed]))) +
    geom_abline(slope = 1, intercept = 0) + th2
  ggsave(paste0(outdir, '/qqplot_null.png'), p, width = 7, height = 5)
}
# qq-plot alternative
{
  n <- dat.alt.plot %>% group_by(nn, gwas_n) %>% summarise(n = n())
  bon_pval_cutoff <- 0.05 / n$n[1] 
  dat.alt.plot <- dat.all.alt %>%
    inner_join(new2, by = 'nn') %>%
    mutate(label = ifelse(bhat.true == 0, 'null', 'signal')) %>%
    arrange(desc(label)) %>%
    group_by(method, gwas_n) %>% 
    mutate(p_exp = rank(pval, ties.method = 'random') / (n() + 1)) %>%
    ungroup() %>%
    mutate(pval = pmax(pval, p.min)) 
  p <- dat.alt.plot %>% ggplot() +
    geom_point(aes(x = -log10(p_exp), -log10(pval), color = label), alpha = 0.3) +
    facet_grid(gwas_n ~ nn2, scales = 'free') +
    geom_abline(slope = 1, intercept = 0) +
    geom_hline(yintercept = -log10(bon_pval_cutoff), linetype = 2) +
    xlab(expression(paste(-log[10], p[expected]))) + 
    ylab(expression(paste(-log[10], p[observed]))) +
    scale_color_manual(values = label_color) + th2 +
    theme(legend.title = element_blank())
  ggsave(paste0(outdir, '/qqplot_alt.png'), p, width = 7.5, height = 5)
}

# fpr 
{
  df.fp <- list()
  for(alpha in alphas) {
    fp <- dat.all %>%
      group_by(nn, gwas_n) %>% 
      summarize(fp_frac = mean(pval < alpha), n = n()) %>%
      ungroup() %>%
      mutate(se = sqrt(fp_frac * (1 - fp_frac) / n))
    df.fp[[length(df.fp) + 1]] <- fp %>% mutate(alpha = alpha)
  }
  df.fp <- do.call(rbind, df.fp)
  p <- df.fp %>% 
    inner_join(new2[-1:-2, ], by = 'nn') %>%
    ggplot() +
    geom_bar(aes(x = as.character(gwas_n), y = fp_frac, fill = nn2), stat = 'identity', position = position_dodge(1)) +
    facet_wrap(~alpha, labeller = label_both, scales = 'free_y') +
    geom_hline(data = df.fp %>% select(alpha) %>% distinct(), aes(yintercept = alpha), linetype = 2) +
    geom_errorbar(aes(x = as.character(gwas_n), ymin = fp_frac - 1.96 * se, ymax = fp_frac + 1.96 * se, group = nn2), position = position_dodge(1), width = 0.1) +
    xlab('GWAS sample size') +
    ylab('Fraction of false positives among nulls \n(false positive rate)') +
    scale_fill_manual(values = mtd_color) + th2 +
    theme(legend.title = element_blank())
  ggsave(paste0(outdir, '/fpr.png'), p, width = 9, height = 3.5)
}

# power 
{
  df.pow <- list()
  for(alpha in alphas) {
    pow <- dat.all.alt %>%
      group_by(nn, gwas_n) %>% 
      filter(bhat.true != 0) %>%
      mutate(positive = pval < alpha) %>%
      summarize(power = mean(positive), n = n()) %>%
      ungroup() %>%
      mutate(se = sqrt(power * (1 - power) / n))
    df.pow[[length(df.pow) + 1]] <- pow %>% mutate(alpha = alpha)
  }
  df.pow <- do.call(rbind, df.pow)
  p <- df.pow %>% 
    inner_join(new2[-1:-2, ], by = 'nn') %>%
    ggplot() +
    geom_bar(aes(x = as.character(gwas_n), y = power, fill = nn2), stat = 'identity', position = position_dodge(1)) +
    facet_wrap(~alpha, labeller = label_both, scales = 'free_y') +
    # geom_hline(data = df.fp %>% select(alpha) %>% distinct(), aes(yintercept = alpha), linetype = 2) +
    geom_errorbar(aes(x = as.character(gwas_n), ymin = power - 1.96 * se, ymax = power + 1.96 * se, group = nn2), position = position_dodge(1), width = 0.1) +
    xlab('GWAS sample size') +
    ylab('Fraction of positive among true signals \n(power)') +
    scale_fill_manual(values = mtd_color) + th2 +
    theme(legend.title = element_blank())
  ggsave(paste0(outdir, '/power.png'), p, width = 9, height = 3.5)
}

