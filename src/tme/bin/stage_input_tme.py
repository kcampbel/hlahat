""" TME Input stager
Writes a pipeline yaml from epic manifest and pipeline data
"""
import argparse
import os
import re
import pandas as pd
import yaml
import logging
from importlib.resources import files
from tme import R, data
from commonLib.lib.fileio import package_file_path, find_file
from commonLib.lib.search import locate
from commonLib.lib.munge import get_timestamp

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_folder', help='EPIC pipeline input folder')
    parser.add_argument('-t' ,'--threads', default=8, help='Number of threads')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Output directory')
    
    return parser.parse_args(args)

def tme_config(specimen_id:str, input_folder:str, threads:int):
    # Pipeline datafiles and parameters
    config = {
    'goi_f': package_file_path(data, 'pact_goi.tsv'),
    'geneset_f': package_file_path(data, 'genesets.txt'),
    'geneset_meta_f': package_file_path(data, 'genesets_meta.tsv'),
    'xcell_celltypes_f': package_file_path(data, 'xcell_types.txt'),
    'hotspots_f': package_file_path(data, 'chang2017_hotspots.tsv'),
    'gtf_tsv_f': package_file_path(data, 'Homo_sapiens.GRCh37.87.tsv.gz'),
    'hgnc_f': package_file_path(data, 'hgnc_complete_set_2021-07-01.txt.gz'),
    'threads': threads,
    }
    params = yaml.safe_load(files(data).joinpath('tme_report_params.yml').read_text())
    config.update(params)
    
    # EPIC inputs
    # Sequenza copy number segments file
    fn = f'{specimen_id}_segments.txt'
    # Exclude segments files in the alt directories
    cn_f = find_file(input_folder, fn, pattern='^((?!alt).)*$')

    # Epipope annotator
    fn = f'{specimen_id}*_oncotator_results_Annotated.tsv'
    vars_f = find_file(input_folder, fn)

    config.update(
        {
            'cn_f': cn_f,
            'vars_f': vars_f
        }
    )
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
#        './test_data',
#		'-o', '/tmp',
#    ])
    logging.info(args)
    timestamp = get_timestamp().split('+')[0]

    mf_d = yaml.safe_load(open(args.manifest))
    specimen_id = mf_d['pipeline']['pact_id']
    tme_tsv = f'{args.outdir}/{specimen_id}_tme_{timestamp}.tsv'
    tme_yml = f'{args.outdir}/{specimen_id}_tme_{timestamp}.yml'

    # Write tsv
    df = manifest2tsv(mf_d)
    df['config'] = tme_yml
    logging.info(f'Writing {tme_tsv}')
    df.to_csv(tme_tsv, sep='\t', index=False)

    # Write yml
    mf_d['manifest_f'] = os.path.abspath(args.manifest)
    config = tme_config(specimen_id, args.input_folder, args.threads)
    config.update(mf_d)
    
    logging.info(f'Writing {tme_yml}')
    stream = open(tme_yml, 'w')
    yaml.dump(config, stream)
    stream.close()

    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
