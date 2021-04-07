This is the fourth round of preprocessing. 
We start with `dmri.original.all_covar` and `t1.scaled.all_covar`.

They are made by the following preprocessing steps:

* For T1: 
    1. Scaled score = Original score x IDP-25000
    2. Residual from regressing out pc1 ~ pc10, age_recruitment, sex, age_recruitment ^ 2, age x sex, age ^2 x sex, IDP-25756 ~ IDP-25759.
* For dMRI: Same as T1 but skip step 1.

In the fourth round of preprocessing, we do the following after the steps described above:

* For T1:
    - Regress out the top-1 PC for each group.
    - Keep the PC. Check if the PC loadings are all positive or negative. If not, pulse (it does not happen in practice). If they are all negative, flip the sign of the PC.
    - Group definition: Total, Subcortical-vol, Subcortical-GMvol, Brainstem, Cortical, Cerebellum
    - Inverse normalization
* For dMRI:
    - Regress out the top-1 PC for each group.
    - Keep the PC. Check if the PC loadings are all positive or negative. If not, pulse (it does not happen in practice). If they are all negative, flip the sign of the PC.
    - Group definition: ICVF, ISOVF, OD, FA for TBSS and ProbTrack respectively.
    - Inverse normalization
Also, keep a copy of T1 and dMRI in which we skip the PC adjustment but still do inverse normalization.
    