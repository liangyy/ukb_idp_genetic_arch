Here we re-use the pipeline built [here](https://github.com/liangyy/misc-tools/tree/master/gw_qtl).

Example run.

```
SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/mixqtl/bin/snakemake
pipe=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/run.snmk
$SNMK -s $pipe --configfile config.test.yaml -pn
```