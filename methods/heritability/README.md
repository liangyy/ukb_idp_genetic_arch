Here I provide two approaches to estimate heritability. 

1. Based on R package `EMMA` which perform matrix factorization of GRM once and apply to all phenotypes. We reimplemented the algorithm in Python, see [here](https://github.com/liangyy/misc-tools/tree/master/pyemma). See `run_pyemma.py` script for details. The calculation is done on Washington.
2. Based on `gcta` (for sanity check). 

