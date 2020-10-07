import pandas as pd
import numpy as np
import scipy.stats
from pandas_plink import read_plink1_bin 

# for debugging
import pdb

def load_genotype_from_bedfile(bedfile, indiv_list, load_first_n_samples=None, missing_rate_cutoff=0.5, return_snp=False):
    G = read_plink1_bin(bedfile, verbose=False)
    if indiv_list is None:
        indiv_list = G.sample.to_series().tolist()
        if load_first_n_samples is not None:
            indiv_list = indiv_list[:load_first_n_samples]
            
    geno = G.sel(sample=indiv_list).values
    # filter out genotypes with high missing rate
    missing_rate = np.isnan(geno).mean(axis=0)
    geno = geno[:, missing_rate < missing_rate_cutoff]
    maf = np.nanmean(geno, axis=0) / 2
    # impute genotype missing value
    miss_x, miss_y = np.where(np.isnan(geno))
    geno[(miss_x, miss_y)] = maf[miss_y]
    var_geno = 2 * maf * (1 - maf)
    # keep only genotypes with variance != 0
    to_keep = var_geno != 0
    geno = geno[:, to_keep]
    maf = maf[to_keep]
    var_geno = var_geno[to_keep]
    geno = (geno - 2 * maf) / np.sqrt(var_geno)
    
    if return_snp is True:
        snpid = G.variant.variant.to_series().to_list()
        a0 = G.variant.a0.to_series().to_list()
        a1 = G.variant.a1.to_series().to_list()
        return geno, indiv_list, (snpid, a0, a1)
    return geno, indiv_list

def compute_grm_from_bed(bedfile_pattern, load_first_n_samples=None, missing_rate_cutoff=0.5):
    
    indiv_list = None
    nsnp = 0
    grm = None
    
    for i in range(1, 23):
        geno, indiv_list = load_genotype_from_bedfile(
            bedfile_pattern.format(chr_num=i),
            indiv_list,
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

def obtain_bhat_from_bed(bedfile_pattern, beta_partial, theta_g, 
    load_first_n_samples=None, missing_rate_cutoff=0.5):
    indiv_list = None
    nsnp = 0
    beta_hat = []
    snpid = []
    a0 = []
    a1 = []
    
    for i in range(1, 23):
        geno, indiv_list, snp_info = load_genotype_from_bedfile(
            bedfile_pattern.format(chr_num=i),
            indiv_list,
            load_first_n_samples=load_first_n_samples,
            missing_rate_cutoff=missing_rate_cutoff, 
            return_snp=True
        )
        if i == 22:
            breakpoint() 
        beta_unscaled_i = geno.T @ beta_partial
        nsnp += geno.shape[1]
        
        beta_hat.append(beta_unscaled_i)
        snpid += snp_info[0]
        a0 += snp_info[1]
        a1 += snp_info[2]
        del beta_unscaled_i
        del snp_info
    
    beta_hat = np.concatenate(beta_hat, axis=0)
    beta_hat = theta_g / beta_hat.shape[0] * beta_hat
    
    return beta_hat, snpid, a0, a1
    
        

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
    parser.add_argument('--train_full_model', action='store_true', help='''
        If specified, the script will generate training model 
        without inner loop of cross-validation 
        (args.nfold[1] will not be used).
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
       
    logging.info('Loading phenotypes.')
    pheno_mat, pheno_indiv_info, pheno_col_info = load_phenotype_parquet(args.phenotype_parquet)
    
    logging.info('Computing GRM.')
    grm_cache_file = args.output + '.grm_cache.pkl.gz'
    if os.path.isfile(grm_cache_file):
        with gzip.open(grm_cache_file, 'rb') as f:
            tmp = pickle.load(f)
            grm = tmp['grm']
            grm_indiv_info = tmp['grm_indiv_info']
    else:
        grm, grm_indiv_info = compute_grm_from_bed(args.geno_bed_pattern, load_first_n_samples=args.first_n_indiv)
        with gzip.open(grm_cache_file, 'wb') as f:
            tmp = {
                'grm': grm,
                'grm_indiv_info': grm_indiv_info
            }
            pickle.dump(tmp, f, protocol=4)
    
    logging.info('Finalizing GRM and phenotype matrices.')
    indiv_info = intersection(pheno_indiv_info, grm_indiv_info)
    pheno_mat = rearrange(pheno_mat, pheno_indiv_info, indiv_info, mode='row')
    grm = rearrange(grm, grm_indiv_info, indiv_info, mode='both')
    nsample = grm.shape[0]
    npheno = pheno_mat.shape[1]
    logging.info(
        'There are {} individuals and {} phenotypes in the final evaluation.'.format(
            grm.shape[0],
            pheno_mat.shape[1]
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
        betahat, snpid, ref, alt = obtain_bhat_from_bed(
            args.geno_bed_pattern, load_first_n_samples=args.first_n_indiv,
            beta_partial=beta_partial, theta_g=best_theta_g
        )
        df_meta = pd.DataFrame({'snpid': snpid, 'a0': ref, 'a1':alt})
        df_beta = pd.DataFrame(betahat, columns=pheno_col_info)
        df_beta = pd.concat([df_meta, df_beta], axis=1)
        del df_meta
        
        _, filetype = os.path.splitext(args.output)
        print(filetype)
        if filetype == '.gz':
            df_beta.to_csv(
                args.output, compression='gzip', 
                sep='\t', index=False
            )
        elif filetype == '.parquet':
            df_beta.to_parquet(args.output, index=False)
        
    logging.info('Done.')
        
    
        
