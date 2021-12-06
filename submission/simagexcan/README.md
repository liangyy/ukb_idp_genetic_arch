The list of name tags for the 2nd round.

```
dmri.original.all_covar.no_pc.gw_ridge_beta
dmri.original.all_covar.w_pc.gw_ridge_beta
dmri.original.non_idp_covar.no_pc.gw_ridge_beta
dmri.original.non_idp_covar.w_pc.gw_ridge_beta
dmri.regress.all_covar.no_pc.gw_ridge_beta
dmri.regress.all_covar.w_pc.gw_ridge_beta
dmri.regress.non_idp_covar.no_pc.gw_ridge_beta
dmri.regress.non_idp_covar.w_pc.gw_ridge_beta
t1.original.all_covar.no_pc.gw_elastic_net_beta
t1.original.all_covar.no_pc.gw_ridge_beta
t1.original.all_covar.w_pc.gw_elastic_net_beta
t1.original.all_covar.w_pc.gw_ridge_beta
t1.original.non_idp_covar.no_pc.gw_ridge_beta
t1.original.non_idp_covar.w_pc.gw_ridge_beta
t1.regress.all_covar.no_pc.gw_elastic_net_beta
t1.regress.all_covar.no_pc.gw_ridge_beta
t1.regress.all_covar.w_pc.gw_elastic_net_beta
t1.regress.all_covar.w_pc.gw_ridge_beta
t1.regress.non_idp_covar.no_pc.gw_ridge_beta
t1.regress.non_idp_covar.w_pc.gw_ridge_beta
t1.scaled.all_covar.no_pc.gw_elastic_net_beta
t1.scaled.all_covar.no_pc.gw_ridge_beta
t1.scaled.all_covar.w_pc.gw_elastic_net_beta
t1.scaled.all_covar.w_pc.gw_ridge_beta
t1.scaled.non_idp_covar.no_pc.gw_ridge_beta
t1.scaled.non_idp_covar.w_pc.gw_ridge_beta
[AND MORE ...]
```

The list of name tags for the 3rd round.

```
third_round_dmri.gw_elastic_net_beta
third_round_dmri.gw_ridge_beta
third_round_t1.gw_elastic_net_beta
third_round_t1.gw_ridge_beta
```

The list of name tags for the 4th round.

```
fourth_round.dmri_no_pc.gw_elastic_net_beta
fourth_round.dmri_w_pc.gw_elastic_net_beta
fourth_round.dmri_no_pc.gw_ridge_beta
fourth_round.dmri_w_pc.gw_ridge_beta
fourth_round.t1_no_pc.gw_elastic_net_beta
fourth_round.t1_w_pc.gw_elastic_net_beta
fourth_round.t1_no_pc.gw_ridge_beta
fourth_round.t1_w_pc.gw_ridge_beta
```

New 4th round (wPCadj)

```
bash submit_4th.sh dmri psychiatric_4th_en_wPCadj 2
bash submit_4th.sh t1 psychiatric_4th_en_wPCadj 12
bash submit_4th.sh dmri psychiatric_4th_ridge_wPCadj 2
bash submit_4th.sh t1 psychiatric_4th_ridge_wPCadj 12
bash submit_4th.sh dmri gtex_gwas_4th_en_wPCadj 2
bash submit_4th.sh t1 gtex_gwas_4th_en_wPCadj 12
bash submit_4th.sh dmri gtex_gwas_4th_ridge_wPCadj 2
bash submit_4th.sh t1 gtex_gwas_4th_ridge_wPCadj 12
```

New new 4th round (residual only)

```
bash submit_4th.sh dmri psychiatric_4th_en_residual 2
bash submit_4th.sh t1 psychiatric_4th_en_residual 12
bash submit_4th.sh dmri psychiatric_4th_ridge_residual 2
bash submit_4th.sh t1 psychiatric_4th_ridge_residual 12
bash submit_4th.sh dmri gtex_gwas_4th_en_residual 2
bash submit_4th.sh t1 gtex_gwas_4th_en_residual 12
bash submit_4th.sh dmri gtex_gwas_4th_ridge_residual 2
bash submit_4th.sh t1 gtex_gwas_4th_ridge_residual 12
```

New new 4th round (residual only) with permutation-adjusted z-score

```
bash submit_4th_permz.sh dmri psychiatric_4th_en_residual 2
bash submit_4th_permz.sh t1 psychiatric_4th_en_residual 12
bash submit_4th_permz.sh dmri psychiatric_4th_ridge_residual 2
bash submit_4th_permz.sh t1 psychiatric_4th_ridge_residual 12
```

