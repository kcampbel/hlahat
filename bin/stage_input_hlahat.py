""" Nextflow HLAHAT input stager
Generates a tsv for input into the HLAHAT nextflow pipeline
"""
import pandas as pd
import argparse
import os
import yaml
import logging
from collections import OrderedDict
from importlib.resources import files
from commonLib.lib.fileio import check_paths_exist, package_file_path, find_file
from commonLib.lib.search import locate
from commonLib.lib.munge import get_timestamp

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_dir', help='EPIC pipeline input path')
    parser.add_argument('-o', '--outfile', help='Pipeline tsv to write')
    
    return parser.parse_args(args)

def pipeline_config(specimen_id:str, input_folder:str, seqtype:str):
    # Pipeline datafiles and parameters
    config = OrderedDict()
    
    pact_id = specimen_id.split('_')[0]
    config.update(
        {
            'fastq_1': find_file(input_folder, f'*{pact_id}*{seqtype}*reads1.fastq.gz'),
            'fastq_2': find_file(input_folder, f'*{pact_id}*{seqtype}*reads2.fastq.gz'),
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
#        '-o', '/tmp/hlahat.tsv'
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
    meta = manifest2tsv(mf_d)
    specimen_id = mf_d["pipeline"]["pact_id"]

    # Pipeline config
    df = pd.DataFrame()
    for seqtype in ['normal_dna', 'tumor_dna', 'tumor_rna']:
        cfg = pipeline_config(specimen_id, args.input_dir, seqtype)
        row = meta.copy()
        row.iloc[0]['sample'] = f'{specimen_id}_{seqtype}'
        tmp = pd.concat([row, pd.DataFrame([cfg])], axis=1)
        df = pd.concat([df, tmp])
    
    if not args.outfile:
        args.outfile = f'{specimen_id}_hlahat_{timestamp}.tsv'
    logging.info(f'Writing {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)

    logging.info(f'{os.path.basename(__file__)} finished')
if __name__ == "__main__":
    main()
