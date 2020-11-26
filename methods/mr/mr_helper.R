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
