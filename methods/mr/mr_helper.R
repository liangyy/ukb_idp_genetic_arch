load_idp_gwas = function(file_pattern) {
  idp_gwas = list()
  for(chr_num in 1 : 22) {
    idp_gwas[[length(idp_gwas) + 1]] = arrow::read_parquet(glue::glue(file_pattern, .open = '[', .close = ']'))
  }
  do.call(rbind, idp_gwas)
}

load_snp_meta = function(bim_pattern) {
  snp_meta = list()
  for(chr_num in 1 : 22) {
    snp_meta[[length(snp_meta) + 1]] = read.table(glue::glue(bim_pattern, .open = '[', .close = ']'), header = F) 
  }
  snp_meta = do.call(rbind, snp_meta)
  snp_meta = snp_meta %>% select(V2, V5, V6) %>% rename(rsid = V2, ref = V5, alt = V6)
}

ld_clump_local = function(dat, ld_clump_param, mode = 'idp2pheno') {
  d = data.frame(rsid = dat$SNP, pval = dat$pval.exposure, id = dat$id.exposure)
  out = tryCatch(
    {
      ieugwasr::ld_clump(
        d, 
        clump_kb = ld_clump_param[[mode]]$clump_kb, 
        clump_r2 = ld_clump_param[[mode]]$clump_r2, 
        clump_p = ld_clump_param[[mode]]$clump_p,
        plink_bin = ld_clump_param$plink_executable,
        bfile = ld_clump_param$bfile
      )
    },
    error = function(cond) {
      message('Cannot do LD clump with the current setting.', list2str(ld_clump_param[[mode]]))
      return(NULL)
    }
  )
  if(is.null(out)) {
    return(out)
  }
  keep = paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
  dat[keep, ]
}

list2str = function(ll) {
  str = c()
  for(n in names(ll)) {
      str = c(str, paste0(n, ': ', ll[[n]]))
  }
  paste0(str, collapse = ', ')
}

perf_mr = function(exp_dat, out_dat, ld_clump_param, ld_clump_mode) {
  clumped_exp_dat = ld_clump_local(idp_exp_dat, ld_clump_param, mode = ld_clump_mode)
  if(nrow(clumped_exp_dat) == 0) {
    return(list(data = NA, mr = NA))
  }
  extracted_out_dat = format_data(
    out_dat, type = 'outcome', snps = clumped_exp_dat$SNP
  )
  if(nrow(extracted_out_dat) == 0) {
    return(list(data = NA, mr = NA))
  }
  dat = harmonise_data(idp_exp_dat, gwas_dat)
  res = mr(dat)
  list(data = dat, mr = res)
}

load_pheno_gwas = function(filename, yaml_path) {
  col_yaml = yaml::read_yaml(yaml_path)
  df_gwas = data.table::fread(paste0('zcat ', filename), sep = '\t', data.table = F)
  df_gwas = df_gwas[, names(col_yaml)]
  colnames(df_gwas) = unlist(col_yaml)
  if(! 'effect_size' %in% colnames(df_gwas)) {
    if(sum(c('af', 'sample_size') %in% colnames(df_gwas)) < 2) {
      message('The effect_size column is missing. And cannot impute effect since zscore, af, and sample_size are also missing.')
      quit()
    }
    df_gwas$effect_size = impute_b_from_z(df_gwas$zscore, df_gwas$af, df_gwas$sample_size)
    df_gwas$se = df_gwas$effect_size / df_gwas$zscore
  }
  cols = c('variant_id', 'effect_allele', 'non_effect_allele', 'effect_size', 'se')
  df_gwas = df_gwas[, cols]
  df_gwas = df_gwas[ !is.na(df_gwas$variant_id), ]
  df_gwas = df_gwas[ !duplicated(df_gwas$variant_id), ]
  df_gwas
}