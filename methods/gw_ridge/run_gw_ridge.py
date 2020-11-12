import pandas as pd
import numpy as np
import scipy.stats
from pandas_plink import read_plink1_bin 

# for debugging
import pdb

def load_genotype_from_bedfile(bedfile, indiv_list, snplist_to_exclude, chromosome=None, load_first_n_samples=None, 
    missing_rate_cutoff=0.5, return_snp=False):
    G = read_plink1_bin(bedfile, verbose=False)
    
    if chromosome is not None:
        chr_str = G.chrom[0].values.tolist()
        if 'chr' in chr_str:
            chromosome = 'chr' + str(chromosome)
        else:
            chromosome = str(chromosome)
        G = G.where(G.chrom == chromosome, drop=True)
    
    df_geno_indiv = pd.DataFrame({'indiv': G.sample.to_series().tolist()})
    df_geno_indiv['idx'] = [ i for i in range(df_geno_indiv.shape[0]) ]
    
    if indiv_list is None:
        indiv_list = G.sample.to_series().tolist()
        if load_first_n_samples is not None:
            indiv_list = indiv_list[:load_first_n_samples]
    df_target_indiv = pd.DataFrame({'indiv': indiv_list})
    df_geno_indiv = pd.merge(df_geno_indiv, df_target_indiv, on='indiv').sort_values(by=['idx'])
    if df_geno_indiv.shape[0] != len(indiv_list):
        raise ValueError('There are input individuals that do not appear in BED file.')
    query_indiv_list = df_geno_indiv.indiv.tolist()
    
    
    snpid = G.variant.variant.to_series().to_list()
    snpid = np.array([ s.split('_')[1] for s in snpid ])
    if return_snp is True:
        a0 = G.variant.a0.to_series().to_numpy()
        a1 = G.variant.a1.to_series().to_numpy()       
        chrom = G.variant.chrom.to_series().to_numpy()    
    
    geno = G.sel(sample=query_indiv_list).values

    # re-order to target indiv_list
    geno = geno[match_y_to_x(np.array(query_indiv_list), np.array(indiv_list)), :]
    
    # filter out unwanted snps
    geno = geno[:, ~np.isin(snpid, snplist_to_exclude)]
    if return_snp is True:
        a0 = a0[~np.isin(snpid, snplist_to_exclude)]
        a1 = a1[~np.isin(snpid, snplist_to_exclude)]
        chrom = chrom[~np.isin(snpid, snplist_to_exclude)]
        
    snpid = snpid[~np.isin(snpid, snplist_to_exclude)]
   
    # filter out genotypes with high missing rate
    missing_rate = np.isnan(geno).mean(axis=0)
    geno = geno[:, missing_rate < missing_rate_cutoff]
    if return_snp is True:
        snpid = snpid[missing_rate < missing_rate_cutoff]
        a0 = a0[missing_rate < missing_rate_cutoff]
        a1 = a1[missing_rate < missing_rate_cutoff]
        chrom = chrom[missing_rate < missing_rate_cutoff]
        
    maf = np.nanmean(geno, axis=0) / 2
    
    # impute genotype missing value
    miss_x, miss_y = np.where(np.isnan(geno))
    geno[(miss_x, miss_y)] = maf[miss_y] * 2
    var_geno = 2 * maf * (1 - maf)
    
    # keep only genotypes with variance != 0
    to_keep = var_geno != 0
    geno = geno[:, to_keep]
    if return_snp is True:
        snpid = snpid[to_keep]
        a0 = a0[to_keep]
        a1 = a1[to_keep]
        chrom = chrom[to_keep]
        
    maf = maf[to_keep]
    var_geno = var_geno[to_keep]
    geno = (geno - 2 * maf) / np.sqrt(var_geno)
    
    if return_snp is True:
        return geno, indiv_list, np.sqrt(var_geno), (snpid.tolist(), a0.tolist(), a1.tolist(), chrom.tolist())
    else:
        return geno, indiv_list, np.sqrt(var_geno)

def compute_grm_from_bed(bedfile_pattern, snplist_to_exclude=None, load_first_n_samples=None, missing_rate_cutoff=0.5):
    
    indiv_list = None
    nsnp = 0
    grm = None
    if snplist_to_exclude is None:
        snplist_to_exclude = set([])
    
    read_by_file = False
    if '{chr_num}' in bedfile_pattern:
        read_by_file = True
    
    for i in range(1, 23):
        
        if read_by_file is True:
            geno, indiv_list, _ = load_genotype_from_bedfile(
                bedfile_pattern.format(chr_num=i),
                indiv_list,
                snplist_to_exclude,
                load_first_n_samples=load_first_n_samples,
                missing_rate_cutoff=missing_rate_cutoff
            )
        else:
            geno, indiv_list, _ = load_genotype_from_bedfile(
                bedfile_pattern,
                indiv_list,
                snplist_to_exclude,
                chromosome=i,
                load_first_n_samples=load_first_n_samples,
                missing_rate_cutoff=missing_rate_cutoff
            )
        
        # calc sub-GRM
        M = geno.shape[1]
        grm_now = geno @ (geno.T / M)
        
        # update GRM
        if grm is None:
            grm = np.zeros((len(indiv_list), len(indiv_list)))
        w1 = nsnp / (M + nsnp)
        w2 = 1 - w1
        grm = grm * w1 + grm_now * w2
        nsnp += M
    
    return grm, indiv_list

def load_phenotype_parquet(filename):
    df = pd.read_parquet(filename)
    indiv_list = df.individual.tolist()
    pheno_info = df.columns[1:].tolist()
    mat = df.iloc[:, 1:].values
    del df
    return mat, indiv_list, pheno_info
    
def intersection(l1, l2):
    s1 = set(l1)
    return sorted(list(s1.intersection(set(l2))))

def match_y_to_x(x, y):
    '''
    x, y are 1d np.array 
    y is a subset of x.
    return idx such that x[idx] = y
    '''
    return np.where(y.reshape(y.size, 1) == x)[1]
    
def rearrange(mat, reference_label, target_label, mode):
    rearrange_idx = match_y_to_x(np.array(reference_label), np.array(target_label))
    if mode == 'row':
        return mat[rearrange_idx, :]
    elif mode == 'both':
        return mat[rearrange_idx, :][:, rearrange_idx]
    else:
        raise ValueError('Mode = {} has not been implemented yet.'.format(mode))

def gen_partition(nsample, nfold, seed):
    np.random.seed(seed)
    size = nsample // nfold 
    rest = nsample - nfold * size
    out = []
    for i in range(nfold):
        out += [i] * size
        if i < rest:
            out += [i]
    return np.random.permutation(np.array(out))

def pearson(y1, y2):
    '''
    y1, y2 are 1d np.array
    '''
    return np.corrcoef(y1, y2)[1, 0]

def spearman(y1, y2):
    '''
    y1, y2 are 1d np.array
    '''
    return scipy.stats.spearmanr(y1, y2).correlation

def r2(ypred, yobs):
    '''
    ypred, yobs are n x k np.array with n being sample size
    '''
    mean_total_error = np.power(yobs - yobs.mean(axis=0), 2).mean(axis=0)
    mean_resid_error = np.power(yobs - ypred, 2).mean(axis=0)
    return 1 - mean_resid_error / mean_total_error

def evaluate_performance(ypred, yobs):
    k = ypred.shape[1]
    pearson_col = []
    spearman_col = []
    for i in range(k):
        pearson_col.append(pearson(ypred[:, i], yobs[:, i]))
        spearman_col.append(spearman(ypred[:, i], yobs[:, i]))
    r2_ = r2(ypred, yobs)
    return pd.DataFrame({'R2': r2_, 'Pearson': pearson_col, 'Spearman': spearman_col})

def obtain_bhat_from_bed(bedfile_pattern, beta_partial, theta_g, indiv_list=None, snplist_to_exclude=None, 
    load_first_n_samples=None, missing_rate_cutoff=0.5):
    nsnp = 0
    beta_hat = []
    snpid = []
    a0 = []
    a1 = []
    chrom = []
    if snplist_to_exclude is None:
        snplist_to_exclude = set([])
    
    
    read_by_file = False
    if '{chr_num}' in bedfile_pattern:
        read_by_file = True

    for i in range(1, 23):
        
        if read_by_file is True:
            geno, indiv_list, geno_sd, snp_info = load_genotype_from_bedfile(
                bedfile_pattern.format(chr_num=i),
                indiv_list,
                snplist_to_exclude,
                chromosome=i,
                load_first_n_samples=load_first_n_samples,
                missing_rate_cutoff=missing_rate_cutoff, 
                return_snp=True
            )
        else:
            geno, indiv_list, geno_sd, snp_info = load_genotype_from_bedfile(
                bedfile_pattern,
                indiv_list,
                snplist_to_exclude,
                chromosome=i,
                load_first_n_samples=load_first_n_samples,
                missing_rate_cutoff=missing_rate_cutoff, 
                return_snp=True
            )
        
        beta_unscaled_i = geno.T @ beta_partial / geno_sd[:, np.newaxis]
        nsnp += geno.shape[1]
        
        beta_hat.append(beta_unscaled_i)
        snpid += snp_info[0]
        a0 += snp_info[1]
        a1 += snp_info[2]
        chrom += snp_info[3]
        del beta_unscaled_i
        del snp_info
    
    beta_hat = np.concatenate(beta_hat, axis=0)
    beta_hat = theta_g / beta_hat.shape[0] * beta_hat
    
    return beta_hat, snpid, a0, a1, chrom
    
def load_list(filename):
    o = []
    with open(filename, 'r') as f:
        for i in f:
            o.append(i.strip())
    return list(set(o))    


def load_grm_id(grm_id):
    o = []
    with open(grm_id, 'r') as f:
        for i in f:
            i = i.strip()
            o.append(i.split('\t')[1])
    return o

def load_grm(grm_txt_gz, grm_id):
    grm_indiv = load_grm_id(grm_id)
    nindiv = len(grm_indiv)
    grm_mat = np.zeros((nindiv, nindiv))
    with gzip.open(grm_txt_gz, 'rt') as f:
        for i in f:
            i = i.strip()
            x, y, _, val = i.split('\t')
            x, y, val = int(x), int(y), float(val)
            grm_mat[x - 1, y - 1] = val
            grm_mat[y - 1, x - 1] = val
    return grm_mat, grm_indiv

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_gw_ridge.py', description='''
        Run genome-wide ridge regression.
        We assume the GRM can be fit into memory.
        It outputs either CV-based prediction performance or prediction model.
    ''')
    parser.add_argument('--geno_bed_pattern', help='''
        Genotype file in binary PED format (plink).
        It takes {chr_num} as wildcard.
        If you have all chromosomes in one bed file, no wildcard is needed and
        the script will load one chromosome at a time assuming the genotype file
        has 1 .. 22 chromosomes. 
    ''')
    parser.add_argument('--gcta_grm_prefix', default=None, help='''
        Optional. If it is specified, the script will use
        GCTA GRM format. 
        It assumes there are 
        [gcta_grm_prefix].grm.gz and [gcta_grm_prefix].grm.id files.
        CAUTION: the GRM should be constructed using EXACTLY the same SNP set of
        the --geno_bed_pattern file (no extra SNP filters applied). 
        We highly discourage user using this option if they are unsure about 
        how the GRM is constructed. 
    ''')
    parser.add_argument('--phenotype_parquet', help='''
        Phenotype in parquet format.
    ''')
    parser.add_argument('--rand_seed', type=int, default=1, help='''
        Random seed for numpy 
        (to determine the cross-validation partition). 
    ''')
    parser.add_argument('--nfold', nargs='+', default=None, type=int, help='''
        Set the cross-validation fold for 
        outer CV and nested CV respectively.
    ''')
    parser.add_argument('--theta_g_grid', nargs='+', type=float, default=None, 
        help='''
        Set the cross-validation fold for 
        outer CV and nested CV respectively.
    ''')
    parser.add_argument('--first_n_indiv', type=int, default=None, help='''
        For debugging purpose, , it run with the first N individuals.
        And the rest of the samples will be discarded.
    ''')
    parser.add_argument('--output', help='''
        Report the cross-validated R2, 
        Pearson's correlation and Spearman's correlation.
    ''')
    parser.add_argument('--snplist_to_exclude', default=None, help='''
        the list of SNP to be excluded in analysis (e.g. duplicated SNPs). 
    ''')
    parser.add_argument('--train_full_model', action='store_true', help='''
        If specified, the script will generate training model 
        without inner loop of cross-validation 
        (args.nfold[1] will not be used).
        The output betahat is in original scale.
        Doing `plink2 --score` with 'cols=scoresums' gives the desired 
        PRS with a mean shift (the mean shift is due to 2 * MAF * betahat).
    ''')
    parser.add_argument('--grm_cache', default=None, help='''
        If specified, will try to load [grm_cache]. If the file does not exist,
        will load genotype and cache the GRM to this path.
        Otherwise, will cache the GRM to [output].grm_cache.pkl.gz
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    import os.path
    import gzip, pickle
    from tqdm import tqdm 
    import gw_ridge 
    
    outer_nfold, inner_nfold = args.nfold
    
    if args.snplist_to_exclude is None:
        snplist_to_exclude = set([])
    else:
        snplist_to_exclude = load_list(args.snplist_to_exclude)
        logging.info('To exclude {} SNPs'.format(len(snplist_to_exclude)))
       
    logging.info('Loading phenotypes.')
    pheno_mat, pheno_indiv_info, pheno_col_info = load_phenotype_parquet(args.phenotype_parquet)
    
    logging.info('Computing GRM.')
    if args.gcta_grm_prefix is None:
        if args.grm_cache is not None:
            grm_cache_file = args.grm_cache
        else:
            grm_cache_file = args.output + '.grm_cache.pkl.gz'
        if os.path.isfile(grm_cache_file):
            with gzip.open(grm_cache_file, 'rb') as f:
                tmp = pickle.load(f)
                grm = tmp['grm']
                grm_indiv_info = tmp['grm_indiv_info']
        else:
            grm, grm_indiv_info = compute_grm_from_bed(
                args.geno_bed_pattern, snplist_to_exclude,
                load_first_n_samples=args.first_n_indiv)
            with gzip.open(grm_cache_file, 'wb') as f:
                tmp = {
                    'grm': grm,
                    'grm_indiv_info': grm_indiv_info
                }
                pickle.dump(tmp, f, protocol=4)
    else:
        grm, grm_indiv_info = load_grm(
            args.gcta_grm_prefix + '.grm.gz', 
            args.gcta_grm_prefix + '.grm.id'
        )
        if args.first_n_indiv is not None:
            grm = grm[:, :500][:500, :]
            grm_indiv_info = grm_indiv_info[:500]
    
    logging.info('Finalizing GRM and phenotype matrices.')
    indiv_info = intersection(pheno_indiv_info, grm_indiv_info)
    pheno_mat = rearrange(pheno_mat, pheno_indiv_info, indiv_info, mode='row')
    grm = rearrange(grm, grm_indiv_info, indiv_info, mode='both')
    nsample = grm.shape[0]
    logging.info('-> Filtering out phenotypes with constant values.')
    # add a checker to filter out phenotypes with constant values
    npheno0 = pheno_mat.shape[1]
    pheno_std = pheno_mat.std(axis=0)
    pheno_non_constant_idx = np.where(pheno_std != 0)[0]
    pheno_mat = pheno_mat[:, pheno_non_constant_idx]
    pheno_col_info = [ pheno_col_info[i] for i in pheno_non_constant_idx ] # pheno_col_info[pheno_non_constant_idx]
    npheno = pheno_mat.shape[1]
    logging.info('-> {} phenotypes have been removed since they have constant value.'.format(npheno0 - npheno))
    logging.info(
        'There are {} individuals and {} phenotypes in the final evaluation.'.format(
            nsample,
            npheno
        )
    )
    
    if args.train_full_model is False:
        logging.info('{}-fold CV (outer loop).'.format(outer_nfold))
        solver = gw_ridge.blupRidgeSolver(
            grm=grm, y=pheno_mat, 
            theta_g_grid=args.theta_g_grid, 
            inner_cv_fold=inner_nfold
        )
        sample_partitions = gen_partition(nsample, outer_nfold, args.rand_seed)
        Ypred_collector = []
        Yobs_collector = []
        for i in tqdm(range(outer_nfold)):
            train_idx = np.where(sample_partitions != i)[0]
            test_idx = np.where(sample_partitions == i)[0]
            Ypred, Yobs = solver.cv_train(train_idx, test_idx, rand_seed=args.rand_seed + i + 1)
            Ypred_collector.append(Ypred)
            Yobs_collector.append(Yobs)
        Ypred = np.concatenate(Ypred_collector, axis=0)
        Yobs = np.concatenate(Yobs_collector, axis=0)
        # return R2, pearson coef, spearman coef
        df_res = evaluate_performance(Ypred, Yobs)  
        df_res['phenotype'] = pheno_col_info
        df_res.to_csv(args.output, compression='gzip', sep='\t', index=False)
    else:
        logging.info('Training with {}-fold CV.'.format(outer_nfold))
        solver = gw_ridge.blupRidgeSolver(
            grm=grm, y=pheno_mat, 
            theta_g_grid=args.theta_g_grid, 
            inner_cv_fold=outer_nfold
        )
        beta_partial, best_theta_g = solver.cv_train(rand_seed=args.rand_seed)
        logging.info('Obtaining best betahat from beta_partial, best_theta_g, and genotypes.')
        betahat, snpid, ref, alt, chrom = obtain_bhat_from_bed(
            args.geno_bed_pattern, snplist_to_exclude=snplist_to_exclude,
            indiv_list=indiv_info,
            load_first_n_samples=args.first_n_indiv,
            beta_partial=beta_partial, theta_g=best_theta_g
        )
        df_meta = pd.DataFrame({'snpid': snpid, 'a0': ref, 'a1':alt, 'chr': chrom})
        df_beta = pd.DataFrame(betahat, columns=pheno_col_info)
        df_beta = pd.concat([df_meta, df_beta], axis=1)
        del df_meta
        
        _, filetype = os.path.splitext(args.output)
        if filetype == '.gz':
            df_beta.to_csv(
                args.output, compression='gzip', 
                sep='\t', index=False
            )
        elif filetype == '.parquet':
            df_beta.to_parquet(args.output, index=False)
        
    logging.info('Done.')
        
    
        
