OUTDIR=/gpfs/data/im-lab/nas40t2/yanyul/data/brainxcan_data

# conda activate ukb_idp

echo ---------------------- geno covar  ------------------------
# /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/geno_covar/ukb_idp.geno_cov.chr1.banded.npz
# /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/geno_covar/ukb_idp.geno_cov.chr1.banded.snp_meta.parquet
echo Working on geno covar
subdir=geno_covar
mkdir -p $OUTDIR/$subdir

for chr_num in `seq 1 22`
do 
  echo chr_num = $chr_num
  cp /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/geno_covar/ukb_idp.geno_cov.chr${chr_num}.banded.npz $OUTDIR/$subdir/chr${chr_num}.banded.npz
  cp /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/geno_covar/ukb_idp.geno_cov.chr${chr_num}.banded.snp_meta.parquet $OUTDIR/$subdir/chr${chr_num}.banded.snp_meta.parquet
done

echo -----------------------------------------------------------

echo ---------------------- idp weights ------------------------
# /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/fourth_round.dmri_no_pc.gw_elastic_net_beta.parquet 
subdir=idp_weights
mkdir -p $OUTDIR/$subdir
modalities="t1 dmri"
models="ridge elastic_net"

# do w_pc
idp_type="residual"
idp_type_tag="w_pc"
for model_type in $models
do
  mkdir -p $OUTDIR/$subdir/$model_type
  for idp_modality in $modalities
  do
    echo model_type = $model_type, idp_modality = $idp_modality, idp_type=$idp_type
    cp /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/fourth_round.${idp_modality}_${idp_type_tag}.gw_${model_type}_beta.parquet $OUTDIR/$subdir/$model_type/${idp_type}.${idp_modality}.parquet
  done
done

# do no_pc
idp_type="original"
idp_type_tag="no_pc"
idp_type_tag2="w_pc"
for model_type in $models
do
  mkdir -p $OUTDIR/$subdir/$model_type
  for idp_modality in $modalities
  do
    echo model_type = $model_type, idp_modality = $idp_modality, idp_type=$idp_type
    out=$OUTDIR/$subdir/$model_type/${idp_type}.${idp_modality}.parquet
    input=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/fourth_round.${idp_modality}_${idp_type_tag}.gw_${model_type}_beta.parquet 
    inputpc=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/fourth_round.${idp_modality}_${idp_type_tag2}.gw_${model_type}_beta.parquet 
    python add_pc.py \
      --out $out \
      --input $input \
      --input_pc $inputpc
  done
done

echo -----------------------------------------------------------


echo ---------------------- idp gwas  --------------------------
# /scratch/t.cri.yliang/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.t1_w_pc.chr22/IDP-25900.parquet
subdir=idp_gwas
mkdir -p $OUTDIR/$subdir
modalities="t1 dmri"

# do w_pc
idp_type="residual"
idp_type_tag="w_pc"
for chr_num in `seq 1 22`
do 
  for idp_modality in $modalities
  do 
    echo chr_num = $chr_num, idp_modality = $idp_modality, idp_type=$idp_type
    curdir=$OUTDIR/$subdir/${idp_type}.${idp_modality}.chr${chr_num}
    mkdir -p $curdir
    cp /scratch/t.cri.yliang/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.${idp_modality}_${idp_type_tag}.chr${chr_num}/*.parquet $curdir/
  done
done

# do no_pc
idp_type="original"
idp_type_tag="no_pc"
idp_type_tag2="w_pc"
for chr_num in `seq 1 22`
do 
  for idp_modality in $modalities
  do 
    echo chr_num = $chr_num, idp_modality = $idp_modality, idp_type=$idp_type
    curdir=$OUTDIR/$subdir/${idp_type}.${idp_modality}.chr${chr_num}
    mkdir -p $curdir
    cp /scratch/t.cri.yliang/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.${idp_modality}_${idp_type_tag}.chr${chr_num}/*.parquet $curdir/
    cp /scratch/t.cri.yliang/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.${idp_modality}_${idp_type_tag2}.chr${chr_num}/PC-*.parquet $curdir/
  done
done

echo -----------------------------------------------------------

echo ---------------------- idp weights ------------------------
# /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr1.bim
subdir=idp_gwas/snp_bim
mkdir -p $OUTDIR/$subdir

for chr_num in `seq 1 22`
do 
  echo chr_num = $chr_num
  cp /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr${chr_num}.bim $OUTDIR/$subdir/chr${chr_num}.bim
done
echo -----------------------------------------------------------

echo ---------------------- mr data  ---------------------------
# /gpfs/data/im-lab/nas40t2/yanyul/data/ieugwasr/EUR.bed
# /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/mr/ld_clump_another.yaml
subdir=mr
mkdir -p $OUTDIR/$subdir

pops="EUR SAS AMR EAS AFR"

subsubdir=ieugwasr
mkdir -p $OUTDIR/$subdir/$subsubdir
suffixes="bed fam bim"
for pop in $pops
do 
  for suffix in $suffixes
  do
    echo pop = $pop, suffix = $suffix
    cp /gpfs/data/im-lab/nas40t2/yanyul/data/ieugwasr/$pop.$suffix $OUTDIR/$subdir/$subsubdir/$pop.$suffix
  done
done

cp /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/mr/ld_clump_another.yaml $OUTDIR/$subdir/ld_clump.yaml
echo -----------------------------------------------------------
