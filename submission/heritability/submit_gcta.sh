phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/idp_list.txt
batchsize=2

mkdir -p logs_gcta
mkdir -p batch_list_gcta
cd batch_list_gcta

split -d -a 3 -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list_gcta/`
do
  if [[ -f logs_gcta/run_$i.out ]]
  then
    tmp=`cat logs_gcta/run_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      echo qsub -v BATCH=$i -N gcta_$i run_gcta.qsub
    fi
  else
     :
    # qsub -v BATCH=$i -N gw_lasso_$i run_lasso.qsub
  fi
done
