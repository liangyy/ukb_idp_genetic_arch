library(optparse)

option_list <- list(
    make_option(c("-l", "--clump"), type="character", default=NULL,
                help="Output snp list from plink1.9 --clump",
                metavar="character"),
    make_option(c("-w", "--weight"), type="character", default=NULL,
                help="Input GWAS file for plink clumping",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output weight file",
                metavar="character"),
    make_option(c("-p", "--pval"), type="numeric", default=5e-5,
                help="P-value cutoff",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

library(data.table)

clumped <- fread(opt$clump, header = TRUE, data.table = FALSE)
clumped <- clumped[,-ncol(clumped)]
clumped <- na.omit(clumped[clumped$P <= opt$pval,])

tmpout = paste0(opt$output, '.rs_tmp')
write.table(clumped$SNP, tmpout, row.names = F, col.names = F, quote = F)

system(paste0("awk 'FNR==NR {a[$1]; next}; $2 in a' ", tmpout, " ", opt$weight, " > ", opt$output))
system(paste0('rm ', tmpout))
