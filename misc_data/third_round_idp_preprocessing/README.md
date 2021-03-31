Here is the third round of IDP preprocessing which is for dMRI and T1.
We start with `dmri.original.all_covar` and `t1.scaled.all_covar`.

For dMRI IDPs, we focus on ICVF, OD, ISOVF, FA only.

We perform PC adjustment for each measurement separately.

For T1 IDPs, we perform PC adjustment for cortical, subcortical, 
and cerebellum regions separately. 
We keep the total and brainstem volume ones **asis**.
 