# environment setup
source /vol/bmd/yanyul/miniconda3/etc/profile.d/conda.sh
conda activate py3
source /vol/bmd/yanyul/GitHub/ptrs-ukb/scripts/source_snakemake_on_nucleus.sh

pippath=/vol/bmd/yanyul/GitHub/misc-tools/hapmap3_snps

wkdir=`pwd`
cd $pippath

$MYSNMK -s hapmap.snmk --configfile $wkdir/config.eur.yaml -p > $wkdir/run.log 2>&1
