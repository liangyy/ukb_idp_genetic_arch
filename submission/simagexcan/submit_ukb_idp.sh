# ARGS: middle name of the config (tailored to ukb idp cohort)

for i in `seq 1 22`
do 
  if (( $i < 7)); then nbatch=200; fi
  if (( $i < 17 & $i > 6 )); then nbatch=100; fi
  if (( $i < 23 & $i > 16 )); then nbatch=50; fi
  filename=logs/run_geno_cov_${1}_$i.log
  if [[ -f $filename ]]
  then
    ff=`cat logs/run_geno_cov_${1}_$i.log |tail -n 1|grep shadow`
    if [[ ! -z $ff ]]
    then
      echo qsub -v NBATCH=$nbatch,CHR=$i,NAME=$1 run_geno_cov.qsub
      qsub -v NBATCH=$nbatch,CHR=$i,NAME=$1 run_geno_cov.qsub
    fi
  else
    echo qsub -v NBATCH=$nbatch,CHR=$i,NAME=$1 run_geno_cov.qsub
    qsub -v NBATCH=$nbatch,CHR=$i,NAME=$1 run_geno_cov.qsub
  fi
done
