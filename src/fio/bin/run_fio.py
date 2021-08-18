""" FIO Pipeline """
import pandas as pd
import argparse
import os
import subprocess as sb
import yaml
import logging
from importlib.resources import files
from fio import nextflow
from process_exprs import data as pe_data
from tme import R as tme_R
from commonLib.lib.fileio import check_paths_exist, package_file_path
from commonLib.lib.search import locate

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('pipeline', type=str, default=package_file_path(nextflow, 'main.nf'),
        choices=['fio', 'tme', 'process_exprs'], help='Pipeline to run')
    parser.add_argument('input_dir', help='EPIC pipeline input path')
    parser.add_argument('-m' ,'--manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('-s', '--nextflow-script', default=package_file_path(nextflow, 'main.nf'), help='Nextflow pipeline script')
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('-tsv', '--nextflow-tsv', default='input.tsv', help='Nextflow input tsv to write')
    parser.add_argument('-p', '--nextflow-params', help='Nextflow parameters file')
    parser.add_argument('-e', '--extra', help='Comma separated list of extra args')
    parser.add_argument('-t', '--tracing', action='store_true', help='Enable Nextflow report file generation')
    parser.add_argument('-n', '--dryrun', action='store_true', help='Dry run')
    parser.add_argument('--email', help='Email address for pipeline messaging')
    parser.add_argument('--output_dir', default='output', help='Output directory')
    
    return parser.parse_args(args)

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

def fio_config(input_folder:str):
    def _find_file(input_folder, fn, pattern:str = None):
        hits = list(locate(fn, input_folder))
        if pattern:
            hits = [x for x in hits if re.search(pattern, x)]
        if len(hits) != 1:
            raise ValueError(f'!= 1 input file:\n{fn} {hits}')
        else:
            return hits[0]

    # Manifest
    fn = f'manifest*.yml'
    manifest_f = _find_file(input_folder, fn)
 
    config = {
        'manifest_f': manifest_f
    }
    return config

def main():
    args = parse_args()
#    args = parse_args([
#        'fio',
#        './test_data/PACT056_T_196454',
#        '--nextflow-tsv', '/tmp/input.tsv',
#        '--tracing',
#        '-e', "'-bg,--resume'",
#        '-n',
#        '--manifest', './test_data/manifest_for_testing.yml',
#    ])
    # Check inputs exist
    exists = check_paths_exist([args.input_dir])
    for k,v in exists.items():
        if not v:
            raise FileNotFoundError(k)

    ## Write nextflow input tsv
    # Manifest
    config = fio_config(args.input_dir)
    if not args.manifest:
        args.manifest = config['manifest_f']
    mf_d = yaml.safe_load(open(args.manifest))
    sample = mf_d['pipeline']['pact_id']

    row = {
        'sample': mf_d['pipeline']['pact_id'],
        'manifest': os.path.abspath(args.manifest),
        'input_folder': os.path.abspath(args.input_dir)
    }
    df = pd.DataFrame([row])
    df.to_csv(args.nextflow_tsv, sep='\t', index=False)

    ## Nextflow command
    if args.pipeline == 'fio':
        nextflow_script = package_file_path(nextflow, 'main.nf')
        cmd = nextflow_cmd(nextflow_script, args.nextflow_tsv, args.nextflow_params, args.tracing, args.extra)
        # process_exprs
        cmd.extend(['--tcga_gtex_map', package_file_path(pe_data, 'tcga_gtex.tsv')])
        # tme
        cmd.extend(['--rmd', package_file_path(tme_R, '')])
    if args.pipeline == 'tme':
        nextflow_script = package_file_path(nextflow, 'tme/main.nf')
        cmd = nextflow_cmd(nextflow_script, args.nextflow_tsv, args.nextflow_params, args.tracing, args.extra)
        cmd.extend(['--rmd', package_file_path(tme_R, '')])
    if args.pipeline == 'process_exprs':
        nextflow_script = package_file_path(nextflow, 'process_exprs/main.nf')
        cmd = nextflow_cmd(nextflow_script, args.nextflow_tsv, args.nextflow_params, args.tracing, args.extra)
        cmd.extend(['--tcga_gtex_map', package_file_path(pe_data, 'tcga_gtex.tsv')])
    if args.email:
        cmd.extend(['--email', args.email])
    cmd.extend(['--output', args.output_dir])

    if args.dryrun:
        print(' '.join(cmd))
    else:
        job = sb.run(cmd)

if __name__ == "__main__":
    main()
