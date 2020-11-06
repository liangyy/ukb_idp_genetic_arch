# args1: phenotype parquet file
# args2: grm prefix 
# args3: output prefix

source /vol/bmd/yanyul/miniconda3/etc/profile.d/conda.sh
conda activate ukb_idp

export PYTHONPATH=/vol/bmd/yanyul/GitHub/misc-tools/pyemma:/vol/bmd/yanyul/GitHub/misc-tools/pyutil

python /vol/bmd/yanyul/GitHub/misc-tools/pyemma/run_pyemma.py \
  --grm $2 \
  --grm_cache $2.pyemma_cache \
  --y_table $1 individual \
  --reml \
  --output $3.tsv.gz > $3.log 2>&1
# --reml \

