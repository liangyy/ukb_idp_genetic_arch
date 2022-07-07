if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='calc_xxt.py', description='''
        Calculate XX' / M
    ''')
    parser.add_argument('--geno_bed_pattern', help='''
        Genotype file in binary PED format (plink).
        It takes {chr_num} as wildcard.
        If you have all chromosomes in one bed file, no wildcard is needed and
        the script will load one chromosome at a time assuming the genotype file
        has 1 .. 22 chromosomes. 
    ''')
    parser.add_argument('--output_prefix', help='''
        Cache the XX'/M to [output_prefix].xxt.pkl.gz
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
    
    import pickle, gzip
    import numpy as np
    from run_gw_ridge import load_genotype_from_bedfile
    
    xxt = None
    indiv_list = None
    nsnp = 0
    for i in range(1, 23):
        logging.info(f'-> Working on chromosome{i}')
        geno_prefix = args.geno_bed_pattern.format(chr_num=i)
        geno_i, indiv_list, _ = load_genotype_from_bedfile(
            f'{geno_prefix}.bed', indiv_list, snplist_to_exclude=set([]), 
            return_snp=False, standardize=False)
        geno_i_mean = np.nanmean(geno_i, axis=0)
        geno_i -= geno_i_mean
        logging.info(f'-> {geno_i.mean().sum()}')
        
        M = geno_i.shape[1]
        xxt_now = geno_i @ (geno_i.T / M)
        
        if xxt is None:
            xxt = np.zeros((len(indiv_list), len(indiv_list)))
        w1 = nsnp / (M + nsnp)
        w2 = 1 - w1
        xxt = xxt * w1 + xxt_now * w2
        nsnp += M
        logging.info(f'M = {M}, nsnp = {nsnp}')
    
    logging.info('Saving')
    xxt_file = args.output_prefix + '.xxt.pkl.gz'
    with gzip.open(xxt_file, 'wb') as f:
        tmp = {
            'xxt': xxt,
            'indiv_list': indiv_list
        }
        pickle.dump(tmp, f, protocol=4)
    
