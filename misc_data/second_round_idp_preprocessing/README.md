This is the second round of IDP preprocessing. 
We try out various approaches to pre-processing our T1 and dMRI IDPs.

The strategies to try out are:

* The way to scale T1 IDPs using head size factor (25000): `orignal` vs `scaled` vs `regress`
* To include IDP technical covariates or not: `all_covar` vs `non_idp_covar`
* To have PCA-based adjustment or not: `no_pc` vs `w_pc`

With this considerations, we have the following matrices:

* `orignal_t1_all_covar_no_pc`
* `orignal_t1_non_idp_covar_no_pc`
* `orignal_t1_all_covar_w_pc`
* `orignal_t1_non_idp_covar_w_pc`
* `scaled_t1_all_covar_no_pc`
* `scaled_t1_non_idp_covar_no_pc`
* `scaled_t1_all_covar_w_pc`
* `scaled_t1_non_idp_covar_w_pc`
* `regress_t1_all_covar_no_pc`
* `regress_t1_non_idp_covar_no_pc`
* `regress_t1_all_covar_w_pc`
* `regress_t1_non_idp_covar_w_pc`
* `orignal_dmri_all_covar_no_pc`
* `orignal_dmri_non_idp_covar_no_pc`
* `orignal_dmri_all_covar_w_pc`
* `orignal_dmri_non_idp_covar_w_pc`
* `regress_dmri_all_covar_no_pc`
* `regress_dmri_non_idp_covar_no_pc`
* `regress_dmri_all_covar_w_pc`
* `regress_dmri_non_idp_covar_w_pc`

To do this systematically, we have three scripts for each of the steps:

* `scale_idp.R`: do scaling
* `regress_out_covar.R`: regressing out covariates
* `adjust_w_pc.R`: adjusting for PCs followed by inverse normalization


The IDP covariates: 

* IDP-25756
* IDP-25757
* IDP-25758
* IDP-25759

Ideally, as suggested by [UKB IDP documentation](https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/brain_mri.pdf) (page 18 and 19), more covariates should be included:

```
Code Short name Modality Description
25921 NewEddy dMRI Whether increased search space in eddy current estimation was used for dMRI
25922 YTranslation dMRI Standard deviation of apparent translation in the Y axis as measured by eddy
25923 TErfMRI rfMRI Echo Time for the rfMRI
25924 TEtfMRI tfMRI Echo Time for the tfMRI
25925 T1Scaling T1 Intensity scaling for T1
25926 T2FLAIRScaling T2 Intensity scaling for T2_FLAIR
25927 SWIScaling SWI Intensity scaling for SWI
25928 dMRIScaling dMRI Intensity scaling for dMRI
25929 rfMRIScaling rfMRI Intensity scaling for rfMRI
25930 tfMRIScaling tfMRI Intensity scaling for tfMRI
```

But we leave it for the future but NOT now.
