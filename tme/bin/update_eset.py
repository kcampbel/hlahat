#!/usr/bin/env python
""" Update expression set, metadata, and eset log
""" 
import os
import sys
import pandas as pd 
import argparse
import logging
from datetime import datetime
import yaml
import re
import subprocess as sb
from process import read_exprs, counts2mat, get_extension

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    #parser.add_argument('-c', '--config', required=True, default='batch_correct.yml', help = 'YAML config file')
    #parser.add_argument('metric', choices=['FPKM', 'TPM', 'counts'], help='Count metric')
    #parser.add_argument('counts', nargs='+', help='count files to process (e.g. *genes.tsv)')
    parser.add_argument('counts', nargs=1, help='count file to process (e.g. *genes.tsv)')
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--metadata_tsv', required=True)
    parser.add_argument('--eset', required=True)
    parser.add_argument('--tcga_gtex_map', required=True)
    parser.add_argument('--pact_eset_log', required=True)
    parser.add_argument('-o', '--outpath', default='rev', help='Output path')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of samples with the same name')

    return parser.parse_args(args)

def manifest_to_meta(manifest:dict):
    """ Converts EPIC manifest fields to a dict
    """
    pi = manifest['pipeline']['patient_info']
    row = {
        'specimen_id': pi['patient'],
        'pact_patient_id': pi['patient'].split('_')[0],
        'patient_id': pi['patient.id'],
        'sample_name': pi['sample.id'],
        'date_of_birth': pi['dob'],
        'study_id': pi['study.id'],
        'tcga_study_code': pi['patient.tumorType'],
    }
    return(row)

def update_metadata_tsv(manifest:dict, meta, tcga_gtex, force:bool = False):
    """ Update metadata flat file 
        Args:
            manifest(dict): EPIC manifest
            meta(pandas): metadata df
            tcga_gtex(pandas): TCGA to GTEX mapping df
            force(bool): Force update if True

        Returns:
            pandas or None if specimen_id is found and force=False
    """ 
    row = manifest_to_meta(manifest)
    specimen_id = row['specimen_id']
    if meta.index.isin([specimen_id]).any():
        if force:
            logging.info(f'--force enabled. Dropping {sampleId} from metadata')
            meta = eset.drop(specimen_id, axis=0)
        else:
            logging.warning(f'{specimen_id} already exists in metadata. Set --force to overwrite.')
            return None
    primary_site = tcga_gtex[['tcga_study_code', 'primary_site']]\
                    [tcga_gtex.tcga_study_code == row['tcga_study_code']]\
                    .drop_duplicates().primary_site.unique()[0]
    row.update({'primary_site': primary_site})
    tmp = pd.DataFrame([row]).set_index('specimen_id')
    out = pd.concat([meta, tmp])
    return(out)

def update_eset_log(manifest:dict, eset_log, eset_path:str, timestamp:str):
    """ Update eset log
        Args:
            manifest(dict): EPIC manifest
            eset_log(pandas): Eset log df
            eset_path(str): Path to updated eset
            timestamp(str): timestamp

        Returns:
            pandas 
    """ 
    mf = manifest_to_meta(manifest)
    specimen_id = mf['specimen_id']
    row = {
        'date': timestamp,
        'specimen_id': specimen_id,
        'eset_path': eset_path,
    }
    tmp = pd.DataFrame([row])
    out = pd.concat([eset_log, tmp], ignore_index=True)
    return(out)
    
def main():
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.DEBUG)
    logging.info(f'Starting {os.path.basename(__file__)}')
    now = datetime.now()
    timestamp = now.strftime("%Y%m%dT%H%M%S")

#    args = parse_args()
#    args = parse_args(
#        [
#        '--manifest', '/home/csmith/git/bioinfo-fio/tme/test/manifest.PACT004_T_560351F.yml',
#        '--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset.tsv',
#        #'--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset_20210625T165418.tsv', # Has PACT004
#        #'--eset', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        '--eset', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding.tsv',
#        '--tcga_gtex_map', '/home/csmith/git/bioinfo-fio/tme/data/tcga_gtex.tsv',
#        '--pact_eset_log', '/home/csmith/git/bioinfo-fio/tme/test/pact_eset_log.tsv',
#        '-o', '/tmp/stage_esets',
#        #'-o', '/tmp/stage_esets',
#        '/home/csmith/git/bioinfo-fio/tme/test/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
#        #'-f', 
#        ]
#    )

    # S3
    args = parse_args(
        [
        '--manifest', '/home/csmith/git/bioinfo-fio/tme/test/manifest.PACT004_T_560351F.yml',
        '--metadata_tsv', 's3://pact-research/csmith/tme/test/meta_eset.tsv',
        #'--metadata_tsv', 's3://pact-research/csmith/tme/test/meta_eset_20210625T165418.tsv', # Has PACT004
        #'--eset', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
        #'--eset', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding.tsv',
        '--eset', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding.parquet.gz',
        '--tcga_gtex_map', 's3://pact-research/csmith/tme/data/tcga_gtex.tsv',
        '--pact_eset_log', 's3://pact-research/csmith/tme/test/pact_eset_log.tsv',
        #'-o', 's3://pact-research/csmith/tme/test/output',
        #'-o', '/tmp/stage_esets',
        '-o', '/tmp/rev',
        '/home/csmith/git/bioinfo-fio/tme/test/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
        #'-f', 
        ]
    )
    #cfg = yaml.safe_load(open(args.config, 'r'))
    manifest = yaml.safe_load(open(args.manifest))
    mf = manifest_to_meta(manifest)
    specimen_id = mf['specimen_id']

    meta = pd.read_csv(args.metadata_tsv, sep='\t', index_col=0)
    eset = read_exprs(args.eset, index_col=0)
    eset_log = pd.read_csv(args.pact_eset_log, sep='\t')
    tcga_gtex = pd.read_csv(args.tcga_gtex_map, sep='\t')

    # Check eset and meta for updates
    eset_new = counts2mat(counts=args.counts, metric='TPM', gene_name=eset.index.name, eset=eset,
        force=args.force)
    meta_new = update_metadata_tsv(manifest, meta, tcga_gtex, args.force)
    
    # Write if eset and manifest have been updated
    if eset_new is not None and meta_new is not None:
        logging.info(f'Updating eset and metadata for {specimen_id}')
        # Write outputs
        if not os.path.exists(args.outpath):
            os.makedirs(args.outpath, exist_ok=True)

        # Metadata
        prefix, ext = get_extension(args.metadata_tsv)
        fp = f'{args.outpath}/{prefix}_{timestamp}{ext}'
        fp_orig = f'{args.outpath}/{prefix}{ext}'
        logging.info(f'Writing {fp} and {fp_orig}')
        meta_new.to_csv(fp, sep='\t')
        meta_new.to_csv(fp_orig, sep='\t')

        # Eset
        prefix, ext = get_extension(args.eset)
        #fn = f'{prefix}_{timestamp}{ext}'
        fp = f'{args.outpath}/{prefix}_{timestamp}{ext}'
        fp_orig = f'{args.outpath}/{prefix}{ext}'
        logging.info(f'Writing {fp} and {fp_orig}')
        if 'parquet' in fp:
            eset_new.to_parquet(fp, compression='gzip')
            eset_new.to_parquet(fp_orig, compression='gzip')
        else:
            eset_new.to_csv(fp, sep='\t')
            eset_new.to_csv(fp_orig, sep='\t')

        # Log
        eset_log_new = update_eset_log(manifest, eset_log, fp, timestamp)
        prefix, ext = get_extension(args.pact_eset_log)
        fp = f'{args.outpath}/{prefix}_{timestamp}{ext}'
        fp_orig = f'{args.outpath}/{prefix}{ext}'
        logging.info(f'Writing {fp}')
        eset_log_new.to_csv(fp, sep='\t', index=False)
        eset_log_new.to_csv(fp_orig, sep='\t', index=False)
    else:
        logging.info(f'{specimen_id} found in eset or metadata and --force not specified. Exiting...')
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == '__main__':
    main()