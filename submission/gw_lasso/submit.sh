phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/idp_list.txt
batchsize=8

mkdir -p logs
mkdir -p batch_list
cd batch_list

split -d -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list/`
do
  echo qsub -v BATCH=$i -N gw_lasso_$i run.qsub
done
