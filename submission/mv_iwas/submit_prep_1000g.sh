mkdir -p logs
for i in `seq 1 22`
do
  # export CHR=$i 
  # export DOWNLOAD=Yes
  # screen -dmS chr$i bash prep_1000g_genotype.qsub
  # sbatch --export=CHR=$i prep_1000g_genotype.sbatch
  qsub -v CHR=$i prep_1000g_genotype.qsub
done
