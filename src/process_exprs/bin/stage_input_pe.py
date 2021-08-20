""" Nextflow stage process_exprs 
Generates a tsv for input into the process_exprs nextflow pipeline
"""
import pandas as pd
import argparse
import os
import subprocess as sb
import yaml
import logging
from importlib.resources import files
from process_exprs import data 
from commonLib.lib.fileio import check_paths_exist, package_file_path
from commonLib.lib.search import locate

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_dir', help='EPIC pipeline input path')
    parser.add_argument('-o', '--outfile', default='process_exprs.tsv', help='Pipeline tsv to write')
    
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

def process_exprs_config(specimen_id:str, tcga_study_code:str, input_folder:str):
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

    # Xena primary site
    tg = pd.read_csv(config['tcga_gtex_map'], sep='\t')
    primary_site = tg[tg.tcga_study_code==tcga_study_code].primary_site.unique()[0]

    # Multiqc
    fn = 'multiqc_data.json'
    multiqc_f = _find_file(input_folder, fn)

    config.update(
        {
            'counts': counts_f,
            'primary_site': primary_site,
            'multiqc': multiqc_f,
        })
    return config

def manifest2tsv(manifest:dict):
    cols = ['pact_id', 'patient_info_patient_id', 'patient_info_study_id', 
            'patient_info_dob', 'patient_info_patient_tumorType']
    mf_d = manifest['pipeline']
    df = pd.json_normalize(mf_d, sep='_')
    df.columns = df.columns.str.replace('.', '_', regex=False)
    df = df[cols]
    df.insert(0, 'sample', mf_d['pact_id'])
    return(df)

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
#    args = parse_args([
#        './test_data/PACT056_T_196454/manifest.yml',
#        './test_data/PACT056_T_196454',
#        '-o', '/tmp/process_exprs.tsv'
#    ])
    # Check inputs exist
    exists = check_paths_exist([args.manifest, args.input_dir])
    for k,v in exists.items():
        if not v:
            raise FileNotFoundError(k)

    ## Write nextflow input tsv
    # Manifest
    mf_d = yaml.safe_load(open(args.manifest))
    df = manifest2tsv(mf_d)
    df['manifest'] = os.path.abspath(args.manifest)
    df['input_folder'] = os.path.abspath(args.input_dir)
    tcga_study_code = mf_d['pipeline']['patient_info']['patient.tumorType']

    ## process_exprs
    pe_cfg = process_exprs_config(df['sample'].iloc[0], tcga_study_code, args.input_dir)
    df['primary_site'] = pe_cfg['primary_site']
    df['gene_counts'] = pe_cfg['counts']
    df['multiqc'] = pe_cfg['multiqc']

    logging.info(f'Writing {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
