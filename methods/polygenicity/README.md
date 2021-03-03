We use the pipeline at [here](https://github.com/liangyy/misc-tools/tree/master/est_polygenicity)

# S-BayesS

* Step 1: Pre-calculte LD matrix (rule `all_ldm_chr` for `chr_num = 1 .. 22`).
* Step 2: Run through all IDPs in parallel.

# SLD4M

* Step 1: Format GWAS into MATLAB MAT file (chisq and rsid).
* Step 2: Run SLD4M.

We also want to run some other complex trait GWASs to use as reference. 
See config `config.test_sld4m_txt.yaml` for external GWASs.
