import numpy as np
import pandas as pd

BASE_PAIR = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def get_partitions(total, npart):
    '''
    Return the partition encoded by [0, 0, 1, 1, 2, 2] for instance.
    And the size of the largest partition.
    '''
    size = total // npart
    remain = total - npart * size
    o = []
    for i in range(npart):
        o += [ i ] * size
        if remain >= i + 1:
            o += [ i ]
    return o, size + 1 if remain > 0 else 0

def prepare_args_for_worker(df_info, partitions, bgen, bgi, chunk_size=20):
    '''
    According to partition, split df_info into len(partitions) part.
    In each part, prepare the argument list for calc_prs_at_worker:
    bgen, bgi, and a data frame with 
    prs_score, snpid, a0, a1.
    '''
    partitions = np.array(partitions)
    df_list = []
    for i in range(partitions.max() + 1):
        df_sub = df_info[ partitions == i ].reset_index(drop=True)
        df_sub.drop(columns='chr', inplace=True)
        df_list.append((bgen, bgi, df_sub, i, chunk_size))
    return df_list

def worker_run_wrapper(args):
    return calc_prs_at_worker(*args)

def calc_prs_at_worker(bgen, bgi, df_info, worker_idx, chunk_size=20):
    '''
    df_info should first 3 columns being snpid, a0, a1
    and the rest columns being prs. 
    '''
    reader = ukb_imp_reader.UKBReader(
        bgen=bgen, 
        bgi=bgi
    )
    snplist = df_info.snpid.tolist()
    nchunk = reader.get_nchunk(snplist, chunk_size=chunk_size)
    step_size = max(1, nchunk // 10)
    out = None
    counter = 0
    for res, _, idx in reader.dosage_generator_by_chunk(snplist, chunk_size=chunk_size):
        dosage, a0, a1, _ = res
        prs_effect_size = df_info.iloc[idx, 3:].values
        check_direction_and_flip(
            prs_effect_size, 
            target=(a0, a1), 
            current=(df_info.iloc[idx, :].a0.tolist(), df_info.iloc[idx, :].a1.tolist()) 
        )
        # dosage: nsnp x nindiv; prs_effect_size: nsnp x ntrait
        out_i = np.einsum('ij,jk->ik', np.rint(dosage.T).astype(np.float32), prs_effect_size).astype(np.float32)
        if out is None:
            out = out_i
        else:
            out += out_i
        counter += 1
        if counter % step_size == 0:
            print('Worker {} has done {} / {}'.format(worker_idx, counter, nchunk))
    return out

def check_direction_and_flip(beta_mat, target, current):
    n = beta_mat.shape[0]
    for i in range(n):
        valid_ret = is_valid(target[0][i], target[1][i], current[0][i], current[1][i])
        if valid_ret == 1:
            beta_mat[i, :] = -beta_mat[i, :]
        elif valid_ret == -1:
            beta_mat[i, :] = 0

def is_valid(a0, a1, b0, b1):
    '''
    Return 0 if exactly match.
    Return 1 if flip.
    Return -1 if not valid not ambiguious.
    '''
    # remove ambiguious first.
    if a0 == BASE_PAIR[a1] or b0 == BASE_PAIR[b1]:
        return -1
    # exact match
    if a0 == b0 and a1 == b1:
        return 0
    # flip
    if a0 == b1 and a1 == b0:
        return 1    
    # compliment match
    if a0 == BASE_PAIR[b0] and a1 == BASE_PAIR[b1]:
        return 0
    # compliment flip
    if a0 == BASE_PAIR[b1] and a1 == BASE_PAIR[b0]:
        return 1    
    # if all above does not return, it has to be invalid.
    return -1

def read_sample_file_as_list(fn):
    '''
    Load bgen sample file as list of string.
    So, will discard the first 2 rows.
    '''
    o = []    
    with open(fn, 'r') as f:
        next(f)
        next(f)
        for i in f:
            o.append(i.strip().split(' ')[0])
    return o

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_prs.py', description='''
        Given UKB BGEN and a PRS parquet file, output the 
        polygenic risk score of all UKB individuals.
    ''')
    parser.add_argument('--ukb_bgen_pattern', help='''
        UKB imputed genotype (v3).
        It takes {chr_num} as wildcard.
        We load genotype dosage as float32 to save memory.
    ''')
    parser.add_argument('--ukb_bgi_pattern', help='''
        The bgi files that are associated with bgen (for rbgen and bgenix).
        It takes {chr_num} as wildcard.
    ''')
    parser.add_argument('--ukb_sample_file', help='''
        Sample file.
    ''')
    parser.add_argument('--prs_parquet', help='''
        PRS parquet (nsnp x ntrait).
        We load PRS score as float32 to save memory.
    ''')
    parser.add_argument('--snp_batch_size', type=int, default=20, help='''
        How many variants to retrieve each time.  
        Recommended to use 20.
    ''')
    parser.add_argument('--nthread', type=int, default=2, help='''
        We use multithreading to speed up the process (mostly I/O).
        Specify the number of processors here.
    ''')
    parser.add_argument('--output', help='''
        The output is saved as parquet.
    ''')
    parser.add_argument('--ukb_imp_reader_path', help='''
        Path to ukb_imp_reader script.
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
    
    from multiprocessing import Pool
    
    sys.path.insert(0, args.ukb_imp_reader_path)
    import ukb_imp_reader
    
    nthread = args.nthread
    
    logging.info('Loading PRS scores.')
    df_var = pd.read_parquet(args.prs_parquet)
    trait_list = df_var.columns[4:].tolist()
    df_var = pd.concat([df_var.iloc[:, :4], df_var.iloc[:, 4:].astype(np.float32)], axis=1)
    
    logging.info('Looping over chromosomes.')
    prs = None
    for i in range(1, 23):
        logging.info(f'Working on chromosome {i}')
        
        logging.info(f'Chr {i}: Preparing data.')
        df_sub = df_var[ df_var.chr == str(i) ].reset_index(drop=True)
        if df_sub.shape[0] == 0:
            logging.info(f'Chr {i}: No SNP. Skip.')
            continue
        # chr_list = df_sub.chr.tolist()
        # snp_list = df_sub.snpid.tolist()
        # effect_alleles = df_sub.a1.toilst()
        # score_mat = df_sub.iloc[4:].astype(np.float32).values
        
        logging.info(f'Chr {i}: Dividing the job.')
        nsnp = df_sub.shape[0]
        partitions, task_size = get_partitions(total=nsnp, npart=nthread)
        logging.info(f'Chr {i}: Each thread needs to process {task_size} SNPs.')
        args_by_worker = prepare_args_for_worker(
            df_sub, partitions,
            args.ukb_bgen_pattern.format(chr_num=i),
            args.ukb_bgi_pattern.format(chr_num=i),
            chunk_size=args.snp_batch_size
        )
        with Pool(nthread) as pool:
            res = pool.map(
                worker_run_wrapper, args_by_worker
            )
        # we expect res as a list of nindiv x ntrait arrays
        res = np.array(res) # nwork x nindiv x ntrait
        res = res.sum(axis=0) 
        if prs is None:
            prs = res
        else:
            prs += res
    
    df_prs = pd.DataFrame(prs, columns=trait_list)
    df_indiv = pd.DataFrame(
        { 'indiv': read_sample_file_as_list(args.ukb_sample_file) }, 
    ) 
    df_prs = pd.concat([df_indiv, df_prs], axis=1)
    df_prs.to_parquet(args.output)   
    
    logging.info('Done.')
    
    
        
