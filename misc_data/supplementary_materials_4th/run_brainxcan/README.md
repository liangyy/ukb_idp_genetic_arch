This is an example run of brainxcan pipeline.

To submit:

```
# psychiatric
pheno=SCZ_PGC_2020
qsub -v PHENOTYPE=$pheno,CONFIGNAME=psychiatric -N $pheno run.qsub

# gtex-gwas
pheno=pgc.scz2
qsub -v PHENOTYPE=$pheno,CONFIGNAME=gtex_gwas -N $pheno run.qsub
```

