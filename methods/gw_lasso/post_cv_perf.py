
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='post_cv_perf.py', description='''
        Perform post-process on the models trained by 
        run_gw_lasso.R --mode cv_performance. The output is 
        the table of cross-validated performance (TSV.GZ table). 
        Need export PYTHONPATH=path-to-misc-tools/pyutil.
    ''')
    parser.add_argument('--model_list', help='''
        The list of input models (file path).
    ''')
    parser.add_argument('--rename_yaml', default=None, help='''
        (Optional) Sometimes, we may want to rename the phenotype.
        Specify the rename rule in the yaml (curr_name: new_name).
    ''')
    parser.add_argument('--output', help='''
        Output.
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
    from pyutil import load_list, intersection
    
    if args.rename_yaml is not None:
        from pyutil import read_yaml
        rename_dict = read_yaml(args.rename_yaml)
    else:
        rename_dict = None
    
    out = []
    model_files = load_list(args.model_list)
    logging.info('There are {} model files to load.'.format(len(model_files)))
    for model_path in tqdm(model_files):
        model = load_perf(model_path)
        out.append(model)
    
    # clean rename dict so that it only contains keys occurring in out
    rename_dict_new = {}
    for k in rename_dict.keys():
        if k in out.columns:
            rename_dict_new[k] = rename_dict[k]
    rename_dict = rename_dict_new
    
    for i in range(out.shape[0]):
        out.phenotype[i] = rename_dict[out.phenotype[i]]
 
    out.to_csv(args.output, compression='gzip', sep='\t', index=False)
    