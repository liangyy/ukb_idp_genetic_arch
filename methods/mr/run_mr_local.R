# This is a lightweight script to perform MR.
# using IDP GWAS and local phenotype GWAS. 
# Output the MR results in both directions:
# IDP -> Phenotype
# Phenotype -> IDP
# as RDS file.


library(optparse)

option_list <- list(
    make_option(c("-g", "--idp_gwas_pattern"), type="character", default=NULL,
                help="The IDP GWAS files (should contain [chr_num] as wildcards).",
                metavar="character"),
    make_option(c("-s", "--snp_meta"), type="character", default=NULL,
                help="SNP meta information in plink BIM files. (should contain [chr_num] as wildcards)",
                metavar="character"),
    make_option(c("-l", "--ld_clump_yaml"), type="character", default=NULL,
                help="LD clumping dependent files and parameters are specified in this YAML file.",
                metavar="character"),
    make_option(c("-w", "--gwas_file"), type="character", default=NULL,
                help="The path of GWAS summary statistics (assume TSV.GZ)",
                metavar="character"),
    make_option(c("-a", "--gwas_yaml"), type="character", default=NULL,
                help="YAML file to specify the column names of GWAS summary statistics",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output the Mendelian Randomization in two directions.",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

source('mr_helper.R')
library(dplyr)
library(TwoSampleMR)
# library(ieugwasr)

# logging config
logging::basicConfig(level = opt$log_level)

logging::loginfo('Loading LD clumping YAML.')
ld_clump_param = yaml::read_yaml(opt$ld_clump_yaml)

logging::loginfo('Loading phenotype GWAS.')
df_gwas = load_pheno_gwas(opt$gwas_file, opt$gwas_yaml)

logging::loginfo('Loading IDP GWAS.')
idp_gwas = load_idp_gwas(opt$idp_gwas_pattern)
snp_meta = load_snp_meta(opt$snp_meta)
idp_gwas = left_join(idp_gwas, snp_meta, by = c('variant_id' = 'rsid'))

logging::loginfo('Keeping the common SNPs between phenotype and IDP.')
common_snps = intersect(df_gwas$variant_id, idp_gwas$variant_id)
df_gwas = df_gwas[ df_gwas$variant_id %in% common_snps, ]
idp_gwas = idp_gwas[ idp_gwas$variant_id %in% common_snps, ]
logging::loginfo(paste0(nrow(idp_gwas), ' SNPs are kept.'))

logging::loginfo('Preparing IDP and phenotype data.frames.')
idp_df = data.frame(
  SNP = idp_gwas$variant_id,
  beta = idp_gwas$b,
  se = idp_gwas$b_se,
  effect_allele = idp_gwas$alt,
  other_allele = idp_gwas$ref
)
idp_dat = format_data(idp_df)
pheno_df = data.frame(
  SNP = df_gwas$variant_id,
  beta = df_gwas$effect_size,
  se = df_gwas$se,
  effect_allele = df_gwas$effect_allele,
  other_allele = df_gwas$non_effect_allele
)
pheno_dat = format_data(pheno_df)

logging::loginfo('Working on IDP -> Phenotype.')
idp2pheno = perf_mr(
  exposure = idp_dat, outcome = pheno_df, 
  ld_clump_param = ld_clump_param, ld_clump_mode = 'idp2pheno'
)

logging::loginfo('Working on Phenotype -> IDP.')
pheno2idp = perf_mr(
  exposure = pheno_dat, outcome = idp_df, 
  ld_clump_param = ld_clump_param, ld_clump_mode = 'pheno2idp'
)

logging::loginfo('Saving results.')
saveRDS(
  list(idp2pheno = idp2pheno, pheno2idp = pheno2idp),
  opt$output
)

logging::loginfo('Done.')
