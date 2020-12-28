mkdir -p logs

for i in `cat trait_list.txt`
do
  fn=logs/harmonize_$i.log
  if [[ -f $fn ]]
  then
    e=`cat $fn|tail -n 1 | grep 'unlock\|shadow'`
    if [[ ! -z $e ]]
    then
      qsub -v CONFIG=$i run.qsub 
    fi
  else
    qsub -v CONFIG=$i run.qsub
  fi
done
