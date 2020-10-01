Here we rely on the BSAIL algorithm implemented at [here](https://github.com/junyangq/snpnet) and the details of the algorithm is described in [this paper](https://www.biorxiv.org/content/10.1101/630079v3).

We first need to format genotype data. Specifically, we need to merge the genotypes of 22 chromosomes (in plink BED format) into one genotype file in PGEN format (the new format for `plink2`). 
This step is done at `../../misc_data/merge_plink_bed.sh` following the post [here](https://www.biostars.org/p/148657/).

To obtain the prediction performance, we follow the procedure similar to `../gw_ridge/`. 
Essentially, we want to obtain the prediction on held-out data.
To do so, we split the data into K fold and at each fold we use the rest K-1 fold to training a lasso/elastic net model where the hyperparameter is determined by an inner round of cross-validation (among the K-1 folds). 
