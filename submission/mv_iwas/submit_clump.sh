for i in `seq 1 22`
do
  qsub -v CHR=$i run_clump_by_chr.qsub
done