""" TME Input stager
Writes a pipeline yaml from epic manifest and pipeline data
"""
import argparse
import os
import re
import yaml
import logging
import pandas as pd
from importlib.resources import files
from process_exprs import data
from commonLib.fileio import check_paths_exist, package_file_path
from commonLib.search import locate

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_folder', help='EPIC pipeline input folder')
    parser.add_argument('-o', '--outfile', default='process_exprs.tsv', help='Pipeline tsv to write')
    
    return parser.parse_args(args)

def pipeline_config(specimen_id:str, tcga_study_code, input_folder:str):
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
    config.update(
        {
            'counts': counts_f,
            'primary_site': primary_site,
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

    ## Write nextflow input tsv
    # Manifest
    mf_d = yaml.safe_load(open(args.manifest))
    df = manifest2tsv(mf_d)
    df['manifest'] = os.path.abspath(args.manifest)
    df['input_folder'] = os.path.abspath(args.input_dir)

    ## Pipeline config
    config = pipeline_config(df['sample'], args.input_folder)
    df['gene_counts'] = config['counts']
    df['primary_site'] = config('primary_site')

    df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    main()
