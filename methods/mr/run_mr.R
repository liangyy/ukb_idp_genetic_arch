# This is a lightweight script to perform MR.
# using IDP GWAS and publicly available GWASs in open GWAS database. 
# Output the MR results in both directions:
# IDP -> Phenotype
# Phenotype -> IDP
# as RDS file.


library(optparse)

option_list <- list(
    make_option(c("-g", "--idp_gwas_pattern"), type="character", default=NULL,
                help="The IDP GWAS files (should contain {chr_num} as wildcards).",
                metavar="character"),
    make_option(c("-s", "--snp_meta"), type="character", default=NULL,
                help="SNP meta information in plink BIM files.",
                metavar="character"),
    make_option(c("-l", "--ld_clump_yaml"), type="character", default=NULL,
                help="LD clumping dependent files and parameters are specified in this YAML file.",
                metavar="character"),
    make_option(c("-f", "--open_gwas_id"), type="character", default=NULL,
                help="The GWAS ID in open GWAS database.",
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

gwas_code = opt$open_gwas_id
logging::loginfo(paste0('Phenotype code is ', gwas_code))

logging::loginfo('Loading IDP GWAS.')
idp_gwas = load_idp_gwas(opt$idp_gwas_pattern)
snp_meta = load_snp_meta(opt$snp_meta)
idp_gwas = left_join(idp_gwas, snp_meta, by = c('variant_id' = 'rsid'))
idp_exp_dat = format_data(
  data.frame(
    SNP = idp_gwas$variant_id, 
    beta = idp_gwas$b, 
    se = idp_gwas$b_se, 
    effect_allele = idp_gwas$alt, 
    other_allele = idp_gwas$ref
  )
)

logging::loginfo('Working on IDP -> Phenotype.')
logging::loginfo('** IDP -> Phenotype: LD clumping')
idp_exp_dat = ld_clump_local(idp_exp_dat, ld_clump_param, mode = 'idp2pheno')

logging::loginfo('** IDP -> Phenotype: Extracting outcome GWAS.')
outcome_dat = extract_outcome_data(
	snps = idp_exp_dat$SNP,
	outcomes = gwas_code
)

logging::loginfo('** IDP -> Phenotype: Running MR.')
dat_forward = harmonise_data(idp_exp_dat, outcome_dat)
res_forward = mr(dat_forward)
res_forward %>% pander::pander(caption = 'IDP -> Phenotype')

logging::loginfo('Working on Phenotype -> IDP.')
logging::loginfo('** Phenotype -> IDP: Loading instrument GWAS.')
exp_dat2 = extract_instruments(
  outcomes = gwas_code,
  clump = T,
  p1 = ld_clump_param$pheno2idp$clump_p,
  r2 = ld_clump_param$pheno2idp$clump_r2,
  kb = ld_clump_param$pheno2idp$clump_kb
)

logging::loginfo('** Phenotype -> IDP: Loading IDP as outcome.')
idp_dat = format_data(
  data.frame(
    SNP = idp_gwas$variant_id, 
    beta = idp_gwas$b, 
    se = idp_gwas$b_se, 
    effect_allele = idp_gwas$alt, 
    other_allele = idp_gwas$ref
  ),
  type = 'outcome',
  snps = exp_dat2$SNP
)

logging::loginfo('** Phenotype -> IDP: Running MR.')
dat_backward = harmonise_data(exp_dat2, idp_dat)
res_backward = mr(dat_backward)
res_backward %>% pander::pander(caption = 'Phenotype -> IDP')

logging::loginfo('Saving results.')
saveRDS(
  list(
    idp2pheno = list(data = dat_forward, mr = res_forward),
    pheno2idp = list(data = dat_backward, mr = res_backward),
    gwas_code = gwas_code
  ),
  opt$output
)

logging::loginfo('Done.')
