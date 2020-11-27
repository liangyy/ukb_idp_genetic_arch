from run_gw_ridge import load_genotype_from_bedfile
from CovConstructor import CovConstructor

def load_mode(mode_list):
    if mode_list[0] not in ['naive', 'cap', 'banded']:
        raise ValueError('Wrong mode.')
    elif mode_list[0] in ['cap', 'banded']:
        if len(mode_list) != 2:
            raise ValueError('Wrong number of parameters for mode.')
        if mode_list[0] == 'cap':
            param = int(mode_list[1])
        else:
            param = float(mode_list[1])
    else:
        param = None
    return mode_list[0], param

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='build_genotype_covariance.py', description='''
        Build genotype covariance matrix.
        It supports pushing small values to zero.
        Or construct banded matrix.
        Need to export PYTHONPATH=path-to-gw_ridge
    ''')
    parser.add_argument('--genotype_bed', help='''
        Genotype in binary PED format.
    ''')
    parser.add_argument('--mode', nargs='+', help='''
        Support mode: 1. naive; 2. cap; 3. banded.
        For naive, it will construct a dense matrix and save it as HDF5.
        For cap, need to set the threshold to push to zero and the result will 
        be saved as 
    ''')
    parser.add_argument('--nbatch', type=int, help='''
        Number of batch to process SNPs.
    ''')
    parser.add_argument('--output_prefix', help='''
        The prefix of output.
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
    
    # sanity check on mode
    mode, param = load_mode(args.mode)
    
    logging.info('Loading genotype matrix.')
    geno_mat, _, _, snp_info = load_genotype_from_bedfile(
        args.genotype_bed,
        indiv_list=None,
        snplist_to_exclude=set([]),
        return_snp=True
    )
    
    logging.info('Computing genotype covariance.')
    constructor = CovConstructor(
        data=geno_mat,
        nbatch=args.nbatch
    )
    constructor.compute_to_disk(
        mode=mode
        param=param,
        output_prefix=args.output_prefix
    )
    
    logging.info('Done.')
    
    
