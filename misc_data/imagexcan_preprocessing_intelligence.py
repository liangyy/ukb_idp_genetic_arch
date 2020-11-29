import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats

def plot_binary(f, fn):
    tmp = sns.countplot(y=f) 
    tmp.figure.savefig(fn)
    plt.clf()

def plot_quant(f, fn):
    plt.xlim(np.nanquantile(f, 0), np.nanquantile(f, 0.95))
    tmp = sns.histplot(x=f, binwidth=1) 
    tmp.figure.savefig(fn)
    plt.clf()


if __name__ == '__main__':
    import pandas as pd
    import matplotlib.pyplot as plt

    # input files
    phenotype = '/gpfs/data/im-lab/nas40t2/Data/UKB/ukbrest-queries/2020-11-24-fluid-intelligence/fluid_intelligence.csv'

    # output files
    fig_outdir = 'imagexcan_preprocessing_output'
    pheno_out = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/imagexcan_fluid_intelligence.parquet'
    

    # TASK
    # construct a cleaner version of fluid intelligence phenotype

    print('TASK: Process fluid intelligence phenotype.')
    df = pd.read_csv(phenotype)
    df.drop_duplicates(subset='eid', inplace=True)
    cols = [ 'instance0', 'instance1', 'instance2', 'instance3', 'field20191' ]
    fintel = None
    for cc in cols:
        if fintel is None:
            fintel = df[cc].values
        else:
            fintel[ np.isnan(fintel) ] = df[cc].values[ np.isnan(fintel) ]
    df['fluid_intelligence'] = fintel
    df = df[ ~ df.handedness.isna() ].reset_index(drop=True)
    
    for col in df.columns:
        if col == 'eid':
            continue
        tmp = df[col].values
        num_na = np.isnan(tmp).sum()
        if len(np.unique(tmp)) == 2:
            plot_binary(tmp, f'{fig_outdir}/pheno_{col}.png')
            ncase = np.nansum(tmp)
            print(f'Binary phenotype {col}, # NA = {num_na}, # case = {ncase}')
        else:
            plot_quant(tmp, f'{fig_outdir}/pheno_{col}.png')
            max_, min_ = np.nanmax(tmp), np.nanmin(tmp)
            print(f'Quantitative phenotype {col}, # NA = {num_na}, min = {min_}, max = {max_}')
   
    
    df.to_parquet(pheno_out, index=False)

    
