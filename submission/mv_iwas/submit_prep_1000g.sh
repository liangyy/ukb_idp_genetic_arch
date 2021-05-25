mkdir -p logs
for i in `seq 1 22`
do
  # export CHR=$i 
  # screen -dmS chr$i bash prep_1000g_genotype.sbatch
  # sbatch --export=CHR=$i prep_1000g_genotype.sbatch
  qsub -v CHR=$i prep_1000g_genotype.qsub
done
