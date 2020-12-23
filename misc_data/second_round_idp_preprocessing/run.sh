# * `orignal_t1_all_covar_no_pc`
# * `orignal_t1_non_idp_covar_no_pc`
# * `orignal_t1_all_covar_w_pc`
# * `orignal_t1_non_idp_covar_w_pc`
# * `scaled_t1_all_covar_no_pc`
# * `scaled_t1_non_idp_covar_no_pc`
# * `scaled_t1_all_covar_w_pc`
# * `scaled_t1_non_idp_covar_w_pc`
# * `regress_t1_all_covar_no_pc`
# * `regress_t1_non_idp_covar_no_pc`
# * `regress_t1_all_covar_w_pc`
# * `regress_t1_non_idp_covar_w_pc`
# * `orignal_dmri_all_covar_no_pc`
# * `orignal_dmri_non_idp_covar_no_pc`
# * `orignal_dmri_all_covar_w_pc`
# * `orignal_dmri_non_idp_covar_w_pc`

mkdir -p output
outdir=output

# idptype=t1
# pve=0.3
# mode1s="scaled regress original"
# mode2s="non_idp_covar all_covar"
# mode3s="w_pc no_pc"
# for mode1 in $mode1s
# do
#   # T1: step1
#   Rscript scripts/scale_idp.R \
#     --mode $mode1 \
#     --type $idptype \
#     --output $outdir/$idptype.$mode1.parquet 
#   for mode2 in $mode2s
#   do 
#     # T1: step2
#     Rscript scripts/regress_out_covar.R \
#       --mode $mode2 \
#       --input_matrix $outdir/$idptype.$mode1.parquet \
#       --output $outdir/$idptype.$mode1.$mode2.parquet 
#     for mode3 in $mode3s
#     do
#       # T1: step3 
#       echo $idptype $mode1 $mode2 $mode3 
#       Rscript scripts/adjust_w_pc.R \
#         --mode $mode3 \
#         --input_matrix $outdir/$idptype.$mode1.$mode2.parquet \
#         --output_prefix $outdir/$idptype.$mode1.$mode2.$mode3 \
#         --pve_cutoff $pve \
#         --annot_type $idptype
#     done
#   done
# done
# 
  

idptype=dmri
pve=0.5
mode1s="regress original"
mode2s="non_idp_covar all_covar"
mode3s="w_pc no_pc"
for mode1 in $mode1s
do
  # dMRI: step1
  Rscript scripts/scale_idp.R \
    --mode $mode1 \
    --type $idptype \
    --output $outdir/$idptype.$mode1.parquet 
  for mode2 in $mode2s
  do 
    # dMRI: step2
    Rscript scripts/regress_out_covar.R \
      --mode $mode2 \
      --input_matrix $outdir/$idptype.$mode1.parquet \
      --output $outdir/$idptype.$mode1.$mode2.parquet 
    for mode3 in $mode3s
    do
      # dMRI: step3 
      echo $idptype $mode1 $mode2 $mode3 
      Rscript scripts/adjust_w_pc.R \
        --mode $mode3 \
        --input_matrix $outdir/$idptype.$mode1.$mode2.parquet \
        --output_prefix $outdir/$idptype.$mode1.$mode2.$mode3 \
        --pve_cutoff $pve \
        --annot_type $idptype
    done
  done
done