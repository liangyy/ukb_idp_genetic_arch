Here is the one-time-use script to prepare data for brainxcan ([link](https://github.com/liangyy/brainxcan)) software.

# Files

* Genotype covariances: 
    - Input: `/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/geno_covar/*`
    - Output Pattern: 
        + `$OUTDIR/geno_cov/chr{chr_num}.banded.npz`
        + `$OUTDIR/geno_cov/chr{chr_num}.banded.snp_meta.parquet`
* IDP weights:
    - Input: `/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/*`
    - Output Pattern:
        + `$OUTDIR/idp_weights/{model_type}/{idp_type}.{idp_modality}.parquet`
        + About `idp_type`: Label `no_pc` as `original` and `w_pc` as `residual`.
        + Need to copy PCs to `no_pc` ones.
* IDP GWASs:
    - Input: `/scratch/t.cri.yliang/ukb_idp/idp_gwas_4th/*` (also zipped at `/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_gwas_4th.zip`)
    - Output pattern: 
        + `$OUTDIR/idp_gwas/{idp_type}.{idp_modality}.chr{chr_num}/{idp_code}.parquet` 
        + About `idp_type`: Label `no_pc` as `original` and `w_pc` as `residual`.
        + Need to copy PCs to `no_pc` ones.  
* IDP GWAS SNPs:
    - Input: `/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/*.bim`
    - Output pattern:
        + `$OUTDIR/idp_gwas/snp_bim/chr{chr_num}.bim`
* MR data:
    - Input:
        + LD clump default YAML: `/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/mr/ld_clump_another.yaml`
        + ieugwasr reference LD panel: `/gpfs/data/im-lab/nas40t2/yanyul/data/ieugwasr/`
    - Output pattern:
        + YAML: `$OUTDIR/mr/ld_clump.yaml`
        + LD panel: `$OUTDIR/mr/ieugwasr/{pop}.[bed|bim|fam]`
* IDP model performance:
    - Input: `/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/*`
        