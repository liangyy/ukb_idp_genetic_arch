#!/bin/sh

#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

module load matlab

source /home/yanyul/conda_init.sh
source /home/yanyul/snmk_init.sh

conda activate py38


cd /project2/haky/yanyul/GitHub/misc-tools/est_polygenicity

$SNMK -s sld4m.snmk --configfile $SLURM_SUBMIT_DIR/config.$MYJOB_NAME.yaml \
  --config mode=all gwas_list=$NAME_LIST \
  -j1 -p \
  > $SLURM_SUBMIT_DIR/logs/sld4m.$MYJOB_NAME.$SLURM_JOBID.log 2>&1
