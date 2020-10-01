Here I implement the ridge regression for genome-wide SNPs.

The script is designed to use the genotype in PLINK binary PED format.
And it outputs the cross-validated prediction performance of the ridge regression, where the hyper-parameter is determined by a nested round of cross-validation within training folds.

Here we follow the BLUP based formula and we assume that the GRM can be fit into memory. 
If this is not the case, consider reduce the precision of the floating number in GRM. 
For instance, reducing from `float64` to `float32` will reduce the memory need by a half. 

The script to test the performance on a grid of `theta_g` where we assign `theta_g` weight to GRM and `1 - theta_g` weight to the identity component. 
`theta_g` varies between 0 and 1 by construction.
The relation between `theta_g` and `sigma2` (in vanilla ridge regression formula) is `sigma2 = (1 - theta_g) M / theta_g` where `M` is the number of SNPs.  