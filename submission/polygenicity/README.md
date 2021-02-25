# Pre-compute LD matrix for S-BayesS

```
for i in `seq 1 5`
do
  qsub -v CHRNUM=$i -l mem=64gb run_sbayess_precomp_ld.qsub
done
for i in `seq 6 16`
do
  qsub -v CHRNUM=$i -l mem=32gb run_sbayess_precomp_ld.qsub
done
for i in `seq 17 22`
do
  qsub -v CHRNUM=$i -l mem=16gb run_sbayess_precomp_ld.qsub
done
``` 