load_gcta_hsq = function(fn) {
  cmd = paste0('cat ', fn, ' | grep "V(G)/Vp"')
  # system(cmd)
  tmp = tryCatch({
    data.table::fread(cmd = cmd, data.table = F, header = F)
  }, warning = function(error) {
    c(NA, NA, NA)
  }
  )
  # if(is.null(tmp)) {
  #   tmp = c(NA, NA, NA)
  # }
  data.frame(h2 = as.numeric(tmp[2]), h2_SE = as.numeric(tmp[3]))
}

# obtained from https://gist.github.com/hakyim/05ede4ba5e51d00196ee0f70e4cd8fa7
## this is robust to correlation between entries
## it uses the fact that the average of Cauchy r.v. is Cauchy regardless of correlation between them
## https://www.cell.com/ajhg/pdfExtended/S0002-9297(19)30002-3
acat = function(pvec) 
{
  TT = sum( tan( (0.5 - pvec) *pi ) )
  .5 - atan(TT / length(pvec)) / pi
}


p2z = function(p, b) {
  sign(b) * abs(qnorm(p / 2, lower.tail = T))
}

z2p <- function(z) {
  return(2 * exp(pnorm(abs(z), lower.tail = F, log.p = T)))
}

load_ldsc_rg = function(fn) {
  tmp = read.table(pipe(paste0('cat ', fn, ' |awk \'{if($1=="" && a==1){a=0};if($1=="p1" || a==1 ){a=1;print $0}}\'')), header = T, stringsAsFactors = F)
  p1 = unlist(lapply(strsplit(basename(tmp$p1), '\\.'), function(x) { 
    ind = (1 : length(x)) 
    ind = ind[ ! ind %in% c(1, 2, length(x) - 1, length(x))]
    x = x[ind]
    paste0(x, collapse = '.')
  }))
  p2 = stringr::str_remove(basename(tmp$p2), '.parquet')
  tmp$p1 = p1
  tmp$p2 = p2
  tmp
}
