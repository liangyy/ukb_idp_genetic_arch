phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/idp_list.txt
batchsize=8

mkdir -p logs
mkdir -p batch_list
cd batch_list

split -d -a 3 -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list/`
do
  if [[ -f logs/run_$i.out ]]
  then
    tmp=`cat logs/run_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      echo qsub -v BATCH=$i -N gw_lasso_$i run.qsub
    fi
  fi
  # qsub -v BATCH=$i -N gw_lasso_$i run.qsub
done
