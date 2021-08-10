#!/usr/bin/env python
""" Update expression matrix.
Updates expression matrix with counts from long format tsvs (e.g. RSEM) and updates
expression set log and metadata using inputs from the sample's EPIC pipeline manifest. 
""" 
import os
import sys
import pandas as pd 
import argparse
import logging
from datetime import datetime
import yaml
import re
from process_exprs.process import counts2mat, manifest_to_meta
from commonLib.lib.fileio import get_extension, file_time, read_exprs, write_exprs

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('counts', nargs=1, help='count file to process (e.g. *genes.tsv)')
    parser.add_argument('--manifest', required=True, help='EPIC pipeline manifest file')
    parser.add_argument('--exprs', required=True, help='Gene expression matrix to update')
    parser.add_argument('--exprs_log', help='Update log file')
    parser.add_argument('--metadata_tsv', help='Specimen metadata')
    parser.add_argument('--tcga_gtex_map', help='TCGA to GTEX metadata mapping file')
    parser.add_argument('-o', '--outpath', default='rev', help='Output path')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of samples with the same name')

    return parser.parse_args(args)

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
            meta = exprs.drop(specimen_id, axis=0)
        else:
            logging.warning(f'{specimen_id} already exists in metadata. Set --force to overwrite.')
            return None
    primary_site = tcga_gtex[['tcga_study_code', 'primary_site']]\
                    [tcga_gtex.tcga_study_code == row['tcga_study_code']]\
                    .drop_duplicates().primary_site.unique()[0]
    row.update({
        'primary_site': primary_site,
        'study': 'PACT',
        'tumor_normal': 'Tumor',
        })
    tmp = pd.DataFrame([row]).set_index('specimen_id')
    out = pd.concat([meta, tmp])
    return(out)

def main():
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.INFO)
    logging.info(f'Starting {os.path.basename(__file__)}')

    args = parse_args()
#    args = parse_args(
#        [
#        '--manifest', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/manifest.PACT004_T_560351F.yml',
#        '--metadata_tsv', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/meta_eset.tsv',
#        #'--metadata_tsv', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/meta_eset_20210625T165418.tsv', # Has PACT004
#        '--tcga_gtex_map', '/home/csmith/git/bioinfo-fio/src/process_exprs/data/tcga_gtex.tsv',
#        #'--exprs', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        '--exprs', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/eset_pact_geneid_proteincoding.tsv',
#        '--exprs_log', '/home/csmith/git/bioinfo-fio/process_exprs/test_data/pact_eset_log.tsv',
#        '-o', '/tmp/updat_exprs',
#        #'-o', '/tmp/stage_esets',
#        '/home/csmith/git/bioinfo-fio/process_exprs/test_data/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
#        #'-f', 
#        ]
#    )


    # S3
#    args = parse_args(
#        [
#        '--manifest', '/home/csmith/git/bioinfo-fio/tme/test/manifest.PACT004_T_560351F.yml',
#        #'--metadata_tsv', 's3://pact-research/csmith/tme/test/meta_eset.tsv',
#        #'--metadata_tsv', 's3://pact-research/csmith/tme/test/meta_eset_20210625T165418.tsv', # Has PACT004
#        #'--tcga_gtex_map', 's3://pact-research/csmith/tme/data/tcga_gtex.tsv',
#        #'--exprs', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        #'--exprs', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding.tsv',
#        '--exprs', 's3://pact-research/csmith/tme/test/eset_pact_geneid_proteincoding.parquet.gz',
#        '--exprs_log', 's3://pact-research/csmith/tme/test/pact_eset_log.tsv',
#        #'-o', 's3://pact-research/csmith/tme/test/output',
#        #'-o', '/tmp/stage_esets',
#        '-o', '/tmp/rev',
#        '/home/csmith/git/bioinfo-fio/tme/test/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
#        #'-f', 
#        ]
#    )
    manifest = yaml.safe_load(open(args.manifest))
    mf = manifest_to_meta(manifest)
    specimen_id = mf['specimen_id']

    # Check meta for updates
    if args.metadata_tsv:
        if not args.tcga_gtex_map:
            raise Exception(f'--tcga_gtex_map required for metadata update. Exiting...')
        meta = pd.read_csv(args.metadata_tsv, sep='\t', index_col=0)
        tcga_gtex = pd.read_csv(args.tcga_gtex_map, sep='\t')
        meta_new = update_metadata_tsv(manifest, meta, tcga_gtex, args.force)
    else:
        meta_new = 'namespace' 

    # Check exprs for updates
    exprs = read_exprs(args.exprs, index_col=0)
    exprs_new = counts2mat(counts=args.counts, metric='TPM', gene_name=exprs.index.name, exprs=exprs,
        force=args.force)
    
    # Write if exprs and manifest have been updated
    if exprs_new is not None and meta_new is not None:
        logging.info(f'Updating exprs for {specimen_id}')
        # Write outputs
        if not os.path.exists(args.outpath):
            os.makedirs(args.outpath, exist_ok=True)

        if args.metadata_tsv:
            # Metadata
            prefix, ext = get_extension(args.metadata_tsv)
            fp = f'{args.outpath}/{os.path.basename(args.metadata_tsv)}'
            logging.info(f'Writing {fp}')
            meta_new.to_csv(fp, sep='\t')

        # Expression 
        prefix, ext = get_extension(args.exprs)
        fp = f'{args.outpath}/{os.path.basename(args.exprs)}'
        logging.info(f'Writing {fp}')
        write_exprs(exprs_new, fp, compression='gzip')
        timestamp = file_time(fp)

        # Log
        if args.exprs_log:
            exprs_log = pd.read_csv(args.exprs_log, sep='\t')
            row = {
                'date': timestamp,
                'specimen_id': specimen_id,
                'exprs_path': fp
            }
            tmp = pd.DataFrame([row])
            exprs_log_new = pd.concat([exprs_log, tmp], ignore_index=True)
            prefix, ext = get_extension(args.exprs_log)
            fp = f'{args.outpath}/{os.path.basename(args.exprs_log)}'
            logging.info(f'Writing {fp}')
            exprs_log_new.to_csv(fp, sep='\t', index=False)
    else:
        logging.info(f'{specimen_id} found in exprs or metadata and --force not specified. Exiting...')
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == '__main__':
    main()
