""" Nextflow HLACN input stager
Generates a tsv for input into the HLACN nextflow pipeline
"""
import pandas as pd
import argparse
import os
import subprocess as sb
import yaml
import logging
from collections import OrderedDict
from importlib.resources import files
from pactescape.hlacn import data 
from commonLib.lib.fileio import check_paths_exist, package_file_path, find_file
from commonLib.lib.search import locate
from commonLib.lib.munge import get_timestamp

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_dir', help='EPIC pipeline input path')
    parser.add_argument('-o', '--outfile', help='Pipeline tsv to write')
    
    return parser.parse_args(args)

def pipeline_config(specimen_id:str, input_folder:str):
    # Pipeline datafiles and parameters
    config = OrderedDict()
    
    # EPIC inputs
    pact_id = specimen_id.split('_')[0]
    config.update(
        {
            'normal_bam': find_file(input_folder, f'*{pact_id}*normal_dna-final.bam'),
            'normal_bai': find_file(input_folder, f'*{pact_id}*normal_dna-final*bai'),
            'tumor_bam': find_file(input_folder, f'*{specimen_id}*tumor_dna-final.bam'),
            'tumor_bai': find_file(input_folder, f'*{pact_id}*tumor_dna-final*bai'),
            'sequenzaModelRData': find_file(input_folder, 'sequenzaModel.RData'),
            'epic_hlatypes': find_file(input_folder, 'hla-final.tsv'),
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
#        './test_data/PACT291_T_840213-001T/manifest.yml',
#        './test_data/PACT291_T_840213-001T',
#        '-o', '/tmp/hlacn.tsv'
#    ])
    timestamp = get_timestamp().split('+')[0]
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

    # Pipeline config
    cfg = pipeline_config(df['sample'].iloc[0], args.input_dir)
    df = pd.concat([df, pd.DataFrame([cfg])], axis=1)
    
    logging.info(f'Writing {args.outfile}')
    if not args.outfile:
        args.outfile = f'hlacn_{timestamp}.tsv'
    df.to_csv(args.outfile, sep='\t', index=False)
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
