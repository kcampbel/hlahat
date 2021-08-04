""" TME Reporting """
import pandas as pd
import argparse
import os
import subprocess as sb
import yaml
import logging
from importlib.resources import files
from process_exprs import data, nextflow
from commonLib.lib.fileio import check_paths_exist, package_file_path
from commonLib.lib.search import locate

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

def pipeline_config(specimen_id:str, input_folder:str):
    def _find_file(input_folder, fn, pattern:str = None):
        hits = list(locate(fn, input_folder))
        if pattern:
            hits = [x for x in hits if re.search(pattern, x)]
        if len(hits) != 1:
            raise ValueError(f'{specimen_id} has != 1 input file:\n{fn} {hits}')
        else:
            return hits[0]

    # Pipeline datafiles and parameters
    config = {
    'tcga_gtex_map': package_file_path(data, 'tcga_gtex.tsv'),
    }
    
    # EPIC inputs
    # Counts for sample
    fn = f'*{specimen_id}*.genes.tsv'
    counts_f = _find_file(input_folder, fn)
    config.update(
        {
            'counts': counts_f,
        }
    )
    return config

def main():
    args = parse_args()
#    args = parse_args([
#        './test_data/PACT056_T_196454/manifest.yml',
#        './test_data/PACT056_T_196454',
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

    ## Write nextflow tsv
    # Add manifest fields
    mf_d = yaml.safe_load(open(args.manifest))
    mf_d = mf_d['pipeline']
    mf_df = pd.json_normalize(mf_d, sep='_')
    mf_df.insert(0, 'sample', mf_d['pact_id'])
    mf_df['manifest'] = os.path.abspath(args.manifest)
    config = pipeline_config(mf_d['pact_id'], args.input_dir)

    # Add RNAseq counts
    mf_df['gene_counts'] = config['counts']

    # Map TCGA to Xena primary site
    tg = pd.read_csv(config['tcga_gtex_map'], sep='\t')
    tcga_study_code = mf_d['patient_info']['patient.tumorType']
    primary_site = tg[tg.tcga_study_code==tcga_study_code].primary_site.unique()[0]
    mf_df['primary_site'] = primary_site

    mf_df.columns = mf_df.columns.str.replace('.', '_', regex=False)
    cols = ['sample', 'gene_counts', 'manifest', 'primary_site']
    out = mf_df[cols]
    out.to_csv(args.nextflow_tsv, sep='\t', index=False)

    # Parameters
    if not args.nextflow_params:
        args.nextflow_params = package_file_path(data, 'pipeline.yml')
    
    # Run nextflow
    cmd = nextflow_cmd(args.nextflow_script, args.nextflow_tsv, args.nextflow_params, args.tracing, args.extra)
    cmd.extend(['--tcga_gtex_map', f'{config["tcga_gtex_map"]}'])
    if args.dryrun:
        print(' '.join(cmd))
    else:
        job = sb.run(cmd)

if __name__ == "__main__":
    main()
