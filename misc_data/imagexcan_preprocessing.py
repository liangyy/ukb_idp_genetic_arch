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
    covariate = '/vol/bmd/data/ukbiobank/psychiatric_traits/2019-12-17_psychiatric-trait_covariates.csv'
    covariate_pc = '/vol/bmd/yanyul/GitHub/ptrs-ukb/output/query_phenotypes_output.csv'
    phenotype = '/vol/bmd/meliao/data/psychiatric_trait_phenotypes/2020-04-10_collected-phenotypes.txt'
    phenotype_handedness = '/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/data/handedness_query.csv'
    phenotype_ptrs = '/vol/bmd/yanyul/GitHub/ptrs-ukb/output/query_phenotypes_cleaned_up.csv'
    idp_train = '/vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet'

    # output files
    fig_outdir = 'imagexcan_preprocessing_output'
    covar_out = '/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/data/imagexcan_covariate_round_1.parquet'
    pheno_out = '/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/data/imagexcan_phenotype_round_1.parquet'
    idp_indiv_out = '/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/data/idp_cohort.txt'
    
    # TASK 1
    # construct a cleaner version of covariates:
    # age, sex, age ^ 2, sex * age, sex * age ^ 2
    
    print('TASK 1: Covariates.')
    df = pd.read_csv(covariate)
    df = df[['eid', 'age_recruitment', 'sex']].rename(columns={'age_recruitment': 'age'})
    df['age_squared'] = df['age'] ** 2
    df['sex_x_age'] = df['age'] * df['sex']
    df['sex_x_age_squared'] = df['age_squared'] * df['sex']
    dfp = pd.read_csv(covariate_pc)
    dfp.rename(columns={'individual': 'eid'}, inplace=True)
    dfp = dfp[['eid'] + [ f'pc{i}' for i in range(1, 11)] ]
    df = pd.merge(df, dfp, on='eid')
    df.to_parquet(covar_out, index=False)

    # TASK 2
    # construct a cleaner version of phenotypes
    # make binary trait really binary trait.
    # make sure one individual per row without duplication
    # and add handedness
    # and also add the 17 quantitative traits composed for UKB PTRS project

    print('TASK 2: Phenotypes.')
    df = pd.read_csv(phenotype, sep='\t')
    df.drop_duplicates(subset='eid', inplace=True)
    dfh = pd.read_csv(phenotype_handedness)
    cols = ['handedness_x_instance_0', 'handedness_x_instance_1', 'handedness_x_instance_2']
    handedness = None
    for cc in cols:
        if handedness is None:
            handedness = dfh[cc].values
        else:
            handedness[ np.isnan(handedness) ] = dfh[cc].values[ np.isnan(handedness) ]
    hdd_binary = handedness
    hdd_binary[handedness == 1] = 0
    hdd_binary[handedness == 2] = 1
    hdd_binary[np.logical_or(handedness == 3, handedness == -3)] = np.nan
    dfh['handedness'] = handedness
    dfh = dfh[ ~ dfh.handedness.isna() ].reset_index(drop=True)
    df = pd.merge(df, dfh[['eid', 'handedness']], on='eid')
    ad = df['parent_AD_score'].values
    ad[ad == 2] = 1
    df['parent_AD_score'] = ad
    dp = df['parent_depression_score'].values
    dp[dp == 2] = 1
    df['parent_depression_score'] = dp
    df.rename(columns={
        'parent_AD_score': 'parent_AD',
        'parent_depression_score': 'parent_depression'
    }, inplace=True)
    
    df_ptrs = pd.read_csv(phenotype_ptrs)
    cols = ['height', 'dbp', 'sbp', 'bmi', 'wbc', 'rbc', 'hb', 'ht', 'mcv', 'mch', 'mchc', 'platelet', 'lymphocyte', 'monocyte', 'neutrophil', 'eosinophil', 'basophil']
    df_ptrs = df_ptrs[['eid'] + cols]
    df = pd.merge(df, df_ptrs, on='eid')
    
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

    # TASK 3
    # extract the IDP cohort for prediction model training
    # they will be excluded for the downstream association test
    
    print('TASK 3: IDP individual list.')
    df = pd.read_parquet(idp_train)
    df[['individual']].to_csv(idp_indiv_out, header=False, index=False)

    # TASK 4
    # construct the list of non related European individuals as the target cohort 
    # to be included for analysis. 
    # NOT for now since I will re-use the British individual list constructed for 
    # PTRS project.


## meta
# cat /vol/bmd/meliao/data/psychiatric_trait_phenotypes/2020-04-10_collected-phenotypes.txt | for i in `seq 2 13`; do cat /vol/bmd/meliao/data/psychiatric_trait_phenotypes/2020-04-10_collected-phenotypes.txt | cut -f $i | sort|uniq -c; done
 # 158159 0.0
 #   5665 1.0
 #  19209 10.0
 #     43 100.0
 #      1 101.0
 #      8 102.0
 #     10 103.0
 #      4 104.0
 #      7 105.0
 #      9 106.0
 #      2 107.0
 #      8 108.0
 #      5 109.0
 #  10652 11.0
 #      4 110.0
 #      3 111.0
 #      5 112.0
 #      1 113.0
 #      5 114.0
 #      2 115.0
 #      2 116.0
 #      3 117.0
 #      2 118.0
 #      1 119.0
 #  18337 12.0
 #      4 120.0
 #      1 121.0
 #      1 122.0
 #      5 123.0
 #      1 124.0
 #      9 125.0
 #      1 126.0
 #      1 129.0
 #   8772 13.0
 #      2 130.0
 #      3 131.0
 #      1 132.0
 #      5 133.0
 #      1 136.0
 #      1 138.0
 #  13404 14.0
 #      3 140.0
 #      2 141.0
 #      2 144.0
 #      1 145.0
 #      1 146.0
 #      1 148.0
 #      1 149.0
 #   9300 15.0
 #      3 150.0
 #   9157 16.0
 #      1 160.0
 #      1 165.0
 #      1 168.0
 #   5269 17.0
 #      1 170.0
 #      1 174.0
 #      6 175.0
 #      1 176.0
 #      1 178.0
 #      1 179.0
 #   8345 18.0
 #      2 182.0
 #      1 183.0
 #      1 184.0
 #      1 185.0
 #   3936 19.0
 #      2 190.0
 #      1 195.0
 #  18018 2.0
 #   7287 20.0
 #      2 200.0
 #      1 203.0
 #      1 205.0
 #   5064 21.0
 #      2 210.0
 #      1 214.0
 #      1 216.0
 #   4362 22.0
 #      1 220.0
 #      2 225.0
 #   2612 23.0
 #      1 235.0
 #   4325 24.0
 #      1 246.0
 #   3174 25.0
 #      1 252.0
 #      1 256.0
 #      1 259.0
 #   2557 26.0
 #   1918 27.0
 #   2480 28.0
 #   1291 29.0
 #  21144 3.0
 #   2957 30.0
 #   1190 31.0
 #   1420 32.0
 #      1 320.0
 #      1 327.0
 #    939 33.0
 #    996 34.0
 #   1120 35.0
 #   1088 36.0
 #    649 37.0
 #    629 38.0
 #      1 380.0
 #    446 39.0
 #  26667 4.0
 #    938 40.0
 #    336 41.0
 #    914 42.0
 #    350 43.0
 #    391 44.0
 #    367 45.0
 #    370 46.0
 #    214 47.0
 #    292 48.0
 #      1 483.0
 #    227 49.0
 #  21499 5.0
 #    450 50.0
 #    154 51.0
 #    231 52.0
 #    134 53.0
 #    174 54.0
 #    145 55.0
 #    211 56.0
 #     83 57.0
 #     90 58.0
 #     65 59.0
 #  30347 6.0
 #    165 60.0
 #     65 61.0
 #     91 62.0
 #     71 63.0
 #     89 64.0
 #     68 65.0
 #     50 66.0
 #     40 67.0
 #     61 68.0
 #     26 69.0
 #  21592 7.0
 #     99 70.0
 #     27 71.0
 #     52 72.0
 #     19 73.0
 #     30 74.0
 #     49 75.0
 #     30 76.0
 #     28 77.0
 #     24 78.0
 #     17 79.0
 #  22839 8.0
 #     57 80.0
 #     20 81.0
 #     17 82.0
 #     16 83.0
 #     37 84.0
 #     16 85.0
 #     18 86.0
 #     18 87.0
 #     20 88.0
 #     13 89.0
 #  16267 9.0
 #     26 90.0
 #      6 91.0
 #      8 92.0
 #      9 93.0
 #      8 94.0
 #     13 95.0
 #     11 96.0
 #      3 97.0
 #      8 98.0
 #      5 99.0
 #      1 weekly_alcohol
 # 111308 0.0
 # 100373 1.0
 #   3140 10.0
 #     29 11.0
 #    537 12.0
 #     12 13.0
 #    140 14.0
 #    244 15.0
 #     34 16.0
 #      3 17.0
 #     23 18.0
 #      1 19.0
 #  93899 2.0
 #    175 20.0
 #     49 21.0
 #      5 22.0
 #      2 23.0
 #      4 24.0
 #     17 25.0
 #      1 26.0
 #     24 28.0
 #  61123 3.0
 #     33 30.0
 #      1 32.0
 #      9 35.0
 #      3 36.0
 #  41988 4.0
 #     15 40.0
 #      1 41.0
 #      4 42.0
 #      1 44.0
 #      1 49.0
 #  24837 5.0
 #      4 50.0
 #      1 55.0
 #  16474 6.0
 #      1 60.0
 #   4177 7.0
 #   4894 8.0
 #      1 80.0
 #    677 9.0
 #      1 daily_coffee
 #  38600 NA
 # 466669 0.0
 #    169 1.0
 #   5745 10.0
 #      3 100.0
 #     96 11.0
 #   1668 12.0
 #      2 120.0
 #    183 13.0
 #    321 14.0
 #      2 140.0
 #   6480 15.0
 #    305 16.0
 #    185 17.0
 #    456 18.0
 #     17 19.0
 #    427 2.0
 #   8295 20.0
 #     13 21.0
 #     94 22.0
 #     50 23.0
 #     35 24.0
 #   2056 25.0
 #     16 26.0
 #     15 27.0
 #     25 28.0
 #    644 3.0
 #   1854 30.0
 #      3 32.0
 #      5 33.0
 #      2 34.0
 #    212 35.0
 #      3 36.0
 #      1 37.0
 #      2 38.0
 #    675 4.0
 #    622 40.0
 #      1 42.0
 #      1 44.0
 #     44 45.0
 #   1967 5.0
 #    102 50.0
 #      1 51.0
 #      5 55.0
 #      1 56.0
 #   1063 6.0
 #     54 60.0
 #      1 65.0
 #      1 66.0
 #    818 7.0
 #      2 70.0
 #   1242 8.0
 #      9 80.0
 #    202 9.0
 #      1 90.0
 #      1 daily_cigarettes
 # 471489 0
 #  31376 1
 #      1 recurrent_depressive_disorder
 # 502283 0
 #    582 1
 #      1 ptsd
 # 502661 0
 #    204 1
 #      1 alzheimers
 # 500762 0
 #   2103 1
 #      1 alcohol_abuse
 # 502792 0
 #     73 1
 #      1 opioid_abuse
 # 502135 0
 #    730 1
 #      1 bipolar_disorder
 # 502604 0
 #    261 1
 #      1 schizophrenia
 # 446001 0
 #  54853 1
 #   2011 2
 #      1 parent_AD_score
 # 458072 0
 # 42865 1
 #  1928 2
 #      1 parent_depression_score
