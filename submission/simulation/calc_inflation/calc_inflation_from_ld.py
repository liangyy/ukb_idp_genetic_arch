if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='calc_inflation_from_ld.py', description='''
        Calculate the terms in the inflation formula, Tr(R'R) and Tr^2(R), 
        using genotype covariance matrix (sparse or dense). 
    ''')
    parser.add_argument('--genotype_covariance', help='''
        The genotype covariance computed in build_genotype_covariance.py
        Accept wildcard {chr_num}.
        Will automatically search for the corresponding meta SNP file.
    ''')
    parser.add_argument('--correlation', action='store_true', help='''
        If specified, it will use genotype correlation instead of genotype 
        covariance.
    ''')
    parser.add_argument('--output', help='''
        Output file (TAB-delimited) including Tr(R'R) and Tr^2(R) for each 
        chromosome. 
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
    import pandas as pd
    from CovConstructor import CovMatrix
    
    logging.info('Calculating by chromosome')
    res = []
    for i in range(1, 23):
        logging.info(f'Working on Chromosome {i}')
        # df_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
        cov_mat = CovMatrix(args.genotype_covariance.format(chr_num=i))
        TrR = cov_mat.eval_trace(param=5000, cor=args.correlation)
        TrRtR = cov_mat.eval_sum_of_squares(param=5000, cor=args.correlation)
        res.append(pd.DataFrame({'chr': [ i ], 'TrR': [ TrR ], 'TrRtR': [ TrRtR ]}))
    res = pd.concat(res, axis=0)
    
    logging.info('Saving as data frames')
    res.to_csv(args.output, index=False, sep='\t')
    
