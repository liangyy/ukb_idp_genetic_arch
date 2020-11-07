phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.txt
batchsize=2

mkdir -p logs_lasso
mkdir -p batch_list_lasso
cd batch_list_lasso

split -d -a 3 -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list_lasso/`
do
  if [[ -f logs_lasso/run_$i.out ]]
  then
    tmp=`cat logs_lasso/run_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      echo qsub -v BATCH=$i -N gw_lasso_$i run_lasso.qsub
    fi
  else
    # :
    qsub -v BATCH=$i -N gw_lasso_$i run_lasso.qsub
  fi
done
