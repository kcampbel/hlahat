""" TME Reporting """
import pandas as pd
import argparse
import os
import subprocess as sb
import yaml
import logging
from importlib.resources import files
from tme import data, nextflow, R
from tme.fileio import check_paths_exist

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_dir', help='EPIC pipeline input path')
    parser.add_argument('-s', '--nextflow-script', default=package_file_path(nextflow, 'main.nf'), help='Nextflow pipeline script')
    parser.add_argument('-tsv', '--nextflow-tsv', default='input.tsv', help='Nextflow input tsv to write')
    parser.add_argument('-p', '--nextflow-params', help='Nextflow parameters file')
    parser.add_argument('-e', '--extra', help='Comma separated list of extra args')
    parser.add_argument('-t', '--tracing', action='store_true', help='Enable Nextflow report file generation')
    parser.add_argument('-n', '--dryrun', action='store_true', help='Dry run')
    
    return parser.parse_args(args)

def package_file_path(package, filename):
    filepath = files(package).joinpath(filename)
    if filepath.exists():
        return str(filepath.absolute())
    else:
        return None
    

def nextflow_cmd(script:str, input_tsv:str, params_file:str, tracing:bool, extra:bool):
    cmd = [
        'nextflow', 'run', script, 
        '--input', input_tsv,
    ]
    if params_file:
        cmd.extend([ '-params-file', params_file])
    if tracing:
        cmd.extend([
            '-with-report',
            '-with-timeline',
            '-with-dag'
        ])
    if extra:
        args = extra.strip("\"'").split(',')
        cmd.extend(args)

    return cmd

def main():
#    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
#    args = parse_args([
#        './test_data/PACT056_T_196454/manifest.yml',
#        './test_data',
#        '--nextflow-tsv', '/tmp/input.tsv',
#        '--tracing',
#        '-e', "'-bg,--resume'",
#        '-n'
#    ])
    # Check inputs exist
    exists = check_paths_exist([args.manifest, args.input_dir])
    for k,v in exists.items():
        if not v:
            raise FileNotFoundError(k)

    # Write nextflow tsv
    mf_d = yaml.safe_load(open(args.manifest))
    mf_d = mf_d['pipeline']
    mf_df = pd.json_normalize(mf_d, sep='_')
    mf_df.insert(0, 'sample', mf_d['pact_id'])
    mf_df['manifest'] = os.path.abspath(args.manifest)
    mf_df['input_folder'] = os.path.abspath(args.input_dir)
    mf_df.columns = mf_df.columns.str.replace('.', '_', regex=False)
    mf_df.to_csv(args.nextflow_tsv, sep='\t', index=False)

    # Run nextflow
    cmd = nextflow_cmd(args.nextflow_script, args.nextflow_tsv, args.nextflow_params, args.tracing, args.extra)
    tme_rmd = R.__path__[0]
    #tme_rmd = package_file_path(R, 'main.Rmd')
    cmd.extend(['--rmd', tme_rmd])
    if args.dryrun:
        print(' '.join(cmd))
    else:
        job = sb.run(cmd)
#        if job.return_code == 0:
#            logging.info('Pipeline completed')
#        else:
#            logging.info(f'Pipeline failed: {job.return_code}')

if __name__ == "__main__":
    main()
