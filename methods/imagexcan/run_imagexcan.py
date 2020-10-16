from util import *
from solver import *

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_imagexcan.py', description='''
        Run association test between observed phenotype and predicted IDPs.
        There are linear regression and logistic regression available.
        The script will go through one phenotype at a time.
    ''')
    parser.add_argument('--phenotype_table', nargs='+', help='''
        In parquet or csv format.
        Specify the filename followed by the column of individual ID.
    ''')
    parser.add_argument('--covariate_table', nargs='+', help='''
        In parquet or csv format.
        Specify the filename followed by the column of individual ID.
    ''')
    parser.add_argument('--phenotype_yaml', help='''
        Specify the list of phenotype to test and the test to use 
        (logistic regression or linear regression).
    ''')
    parser.add_argument('--individual_list', help='''
        The list of individuals to include in the analysis.
    ''')
    parser.add_argument('--individual_list_exclude', default=None, help='''
        The list of individuals to exclude from the analysis.
    ''')
    parser.add_argument('--covariate_yaml', default=None, help='''
        If it is not specified, all covariates will be used as quantitative measure.
        Otherwise, specify the list of covariate to use and 
        the type of the covariate (continuous or categorical).
    ''')
    parser.add_argument('--idp_table', nargs='+', help='''
        In parquet or csv format.
        Specify the filename followed by the column of individual ID.
    ''')
    parser.add_argument('--output', help='''
        A table in csv format where each row is the test result of one
        observed phenotype and predicted IDP pair.
    ''')
    parser.add_argument('--first_30_idp', action='store_true', help='''
        For debugging purpose, if set, only the first 30 IDPs will be tested.
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
    from tqdm import tqdm
    
    logging.info('Loading tables.')
    df_covar = load_covariate(args.covariate_table, args.covariate_yaml)
    df_pheno, list_pheno_info = load_phenotype(
        args.phenotype_table, args.phenotype_yaml
    )
    df_idp = load_idp(args.idp_table)
    indiv_lists = [df_covar.indiv.to_list(), df_pheno.indiv.to_list(), df_idp.indiv.to_list()]
    if args.individual_list is not None:
        indiv_lists.append(load_list(args.individual_list))
    indiv_list = take_intersect(indiv_lists)
    if args.individual_list_exclude is not None:
        indiv_list = exclude_b_from_a(
            a=indiv_list, 
            b=load_list(args.individual_list_exclude)
        )
    indiv_list = sorted(indiv_list)
   
    df_covar = rearrange_rows(df_covar, indiv_list)
    df_pheno = rearrange_rows(df_pheno, indiv_list)
    df_idp = rearrange_rows(df_idp, indiv_list)
    logging.info('There are {} individauls being included.'.format(df_covar.shape[0]))
    
    if args.first_30_idp:
        df_idp = df_idp.iloc[:, :31]

    Covar = get_matrix(df_covar)
    # add intercept to covariate
    Covar = np.concatenate((np.ones((Covar.shape[0], 1)), Covar), axis=1)
    Idp = get_matrix(df_idp)
    idp_cols = df_idp.columns[1:].to_list()

    df_list = []
    for pheno_col, test_type in list_pheno_info.items():
        logging.info(f'Working on {pheno_col}')
        y = df_pheno[pheno_col].to_numpy()
        not_nan = np.logical_not(np.isnan(y))
        if test_type == 'logistic_regression':
            bhat, pval = logistic_regression(y[not_nan], X=Idp[not_nan, :], C=Covar[not_nan, :])
        elif test_type == 'linear_regression':
            bhat, pval = linear_regression(y[not_nan], X=Idp[not_nan, :], C=Covar[not_nan, :])
        df = pd.DataFrame({
            'IDP': idp_cols, 
            'phenotype': pheno_col, 
            'bhat': bhat, 
            'pval': pval
        })
        df_list.append(df)
    df_list = pd.concat(df_list, axis=0)
    
    df_list.to_csv(args.output, index=False)
    
    
    