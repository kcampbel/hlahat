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
from process_exprs import data
from commonLib.lib.fileio import get_extension, file_time, read_exprs, write_exprs, package_file_path

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('counts', nargs=1, help='count file to process (e.g. *genes.tsv)')
    parser.add_argument('--manifest', required=True, help='EPIC pipeline manifest file')
    parser.add_argument('--exprs', required=True, help='Gene expression matrix to update')
    parser.add_argument('--exprs_log', help='Update log file')
    parser.add_argument('--metadata_tsv', help='Specimen metadata')
    parser.add_argument('--multiqc', help='Multiqc json for QC check')
    parser.add_argument('--blacklist', help='Blacklist file to update if fail')
    parser.add_argument('-o', '--outpath', default='rev', help='Output path')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of samples with the same name are present')
    parser.add_argument('-p', '--passthrough', action='store_true', help='Output without overwrite if samples with the same name are present')
    parser.add_argument('--config', default=package_file_path(data, 'pipeline.yml'), help='YAML with QC thresholds')
    parser.add_argument('--tcga_gtex_map', help='TCGA to GTEX metadata mapping file')

    return parser.parse_args(args)

def update_metadata_tsv(manifest:dict, meta, tcga_gtex, force:bool = False, passthrough:bool = False):
    """ Update metadata flat file 
        Args:
            manifest(dict): EPIC manifest
            meta(pandas): metadata df
            tcga_gtex(pandas): TCGA to GTEX mapping df
            force(bool): Force update if True
            passthrough(bool): leave existing entry in meta if True

        Returns:
            pandas or None if specimen_id is found and force=False
    """ 
    row = manifest_to_meta(manifest)
    specimen_id = row['specimen_id']
    if meta.index.isin([specimen_id]).any():
        if force and not passthrough:
            logging.info(f'--force enabled. Overwriting {specimen_id} in metadata')
            meta = meta.drop(specimen_id, axis=0)
        elif passthrough:
            logging.info(f'--passthrough enabled. Retaining existing data for {specimen_id} in metadata')
            return(meta)
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

def multiqc_json_to_df(sample_id:str, multiqc:dict, cutoffs:dict):
    dat = multiqc['report_saved_raw_data']
    yaml_map = {
        'pe_antisense': dat['multiqc_rseqc_infer_experiment'][f'{sample_id}.rseqc_infer_experiment.txt'],
        'pe_sense': dat['multiqc_rseqc_infer_experiment'][f'{sample_id}.rseqc_infer_experiment.txt'],
        #'failed': dat['multiqc_rseqc_infer_experiment'][f'{sample_id}.rseqc_infer_experiment.txt'],
        'PCT_PF_READS_ALIGNED': dat['multiqc_picard_AlignmentSummaryMetrics'][sample_id],
        'PCT_CODING_BASES': dat['multiqc_picard_RnaSeqMetrics'][sample_id],
        'PCT_RIBOSOMAL_BASES': dat['multiqc_picard_RnaSeqMetrics'][sample_id],
    }

    out = pd.DataFrame()
    for k, v in cutoffs.items():
        observed = yaml_map[k][k]
        direction = v[0]
        threshold = float(v[1::])
        if direction == '>':
            pass_ = observed > threshold
        elif direction == '>=':
            pass_ = observed >= threshold
        elif direction == '<':
            pass_ = observed < threshold
        elif direction == '<=':
            pass_ = observed <= threshold
        elif direction == '==':
            pass_ = observed == threshold
        else:
            raise ValueError(f'No directionality specified for {k}:{v}. Format must be metric: <>=value')

        row = {
            'metric': k,
            'threshold': v,
            'observed': observed,
            'pass': pass_,
        }
        tmp = pd.DataFrame([row])
        out = pd.concat([out, tmp]).round(3)

    return(out)

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')

    args = parse_args()
#    args = parse_args(
#        [
#        '--exprs', '/home/csmith/git/bioinfo-fio/test_data/eset_pact_geneid_proteincoding.tsv', # Has PACT004 but not PACT056
##        #'/home/csmith/git/bioinfo-fio/test_data/PACT004_T_560351F/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
##        '/home/csmith/git/bioinfo-fio/test_data/PACT004_T_560351F/RNA_PACT004_T_560351F_tumor_rna.genes.onegene.tsv',
##        '--manifest', '/home/csmith/git/bioinfo-fio/test_data/PACT004_T_560351F/manifest.PACT004_T_560351F.yml',
#        '/home/csmith/git/bioinfo-fio/test_data/PACT056_T_196454/RNA_PACT056_T_196454_tumor_rna.genes.tsv',
#        '--manifest', '/home/csmith/git/bioinfo-fio/test_data/PACT056_T_196454/manifest.yml',
#        '--metadata_tsv', '/home/csmith/git/bioinfo-fio/test_data/meta_pact_xena.tsv',
#        '--blacklist', '/home/csmith/git/bioinfo-fio/test_data/blacklist.txt',
#        #'--multiqc', '/home/csmith/git/bioinfo-fio/test_data/PACT056_T_196454/multiqc_data.json',
#        '--multiqc', '/home/csmith/git/bioinfo-fio/test_data/PACT056_T_196454/multiqc_data_fail.json',
#        '--tcga_gtex_map', '/home/csmith/git/bioinfo-fio/src/process_exprs/data/tcga_gtex.tsv',
#        '--exprs_log', '/home/csmith/git/bioinfo-fio/test_data/pact_eset_log.tsv',
#        '-o', '/tmp/updat_exprs',
#        #'--passthrough',
#        '-f', 
#        ]
#    )

    # S3
#    args = parse_args(
#        [
#        '--manifest', '/home/csmith/git/bioinfo-fio/test/manifest.PACT004_T_560351F.yml',
#        #'--metadata_tsv', 's3://pact-research/csmith/test/meta_eset.tsv',
#        #'--metadata_tsv', 's3://pact-research/csmith/test/meta_eset_20210625T165418.tsv', # Has PACT004
#        #'--tcga_gtex_map', 's3://pact-research/csmith/data/tcga_gtex.tsv',
#        #'--exprs', 's3://pact-research/csmith/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        #'--exprs', 's3://pact-research/csmith/test/eset_pact_geneid_proteincoding.tsv',
#        '--exprs', 's3://pact-research/csmith/test/eset_pact_geneid_proteincoding.parquet.gz',
#        '--exprs_log', 's3://pact-research/csmith/test/pact_eset_log.tsv',
#        #'-o', 's3://pact-research/csmith/test/output',
#        #'-o', '/tmp/stage_esets',
#        '-o', '/tmp/rev',
#        '/home/csmith/git/bioinfo-fio/test/RNA_PACT004_T_560351F_tumor_rna.genes.tsv',
#        #'-f', 
#        ]
#    )
    logging.info(args)
    if args.force and args.passthrough:
        raise Exception('Cannot pass --force and --passthrough')

    manifest = yaml.safe_load(open(args.manifest))
    mf = manifest_to_meta(manifest)
    specimen_id = mf['specimen_id']

    # Check meta for updates
    if args.metadata_tsv:
        if not args.tcga_gtex_map:
            raise Exception(f'--tcga_gtex_map required for metadata update. Exiting...')
        meta = pd.read_csv(args.metadata_tsv, sep='\t', index_col=0)
        tcga_gtex = pd.read_csv(args.tcga_gtex_map, sep='\t')
        meta_new = update_metadata_tsv(manifest, meta, tcga_gtex, args.force, args.passthrough)
    else:
        meta_new = 'namespace' 

    # Check exprs for updates
    exprs = read_exprs(args.exprs, index_col=0)
    exprs_new = counts2mat(counts=args.counts, metric='TPM', gene_name=exprs.index.name, exprs=exprs,
        force=args.force, passthrough=args.passthrough)
    
    # Write if exprs and manifest have been updated
    if exprs_new is not None and meta_new is not None:
        # Write outputs
        if not os.path.exists(args.outpath):
            os.makedirs(args.outpath, exist_ok=True)
        # Check QC
        if args.multiqc:
            # Load YAML and compare to QC results
            mqc = yaml.safe_load(open(args.multiqc))
            config = yaml.safe_load(open(args.config))
            config = config['exprs_qc']
            rna_tumor = manifest['pipeline']['rna']['tumor']
            qc_df = multiqc_json_to_df(rna_tumor, mqc, config)
            blacklist = pd.read_csv(args.blacklist)
            logging.info(f'Sample QC:\n{qc_df}')
            if (qc_df['pass'] == False).any():
                logging.warning(f'{specimen_id} failed QC')
                # Add to blacklist
                if args.blacklist:
                    if not blacklist.specimen_id.isin([specimen_id]).any():
                        logging.info(f'Adding {specimen_id} to blacklist')
                        tmp = pd.DataFrame([{'specimen_id': specimen_id}])
                        blacklist = pd.concat([blacklist, tmp])
                    else:
                        logging.info(f'{specimen_id} already in {args.blacklist}. Skipping')
                else:
                    logging.warning(f'Skipping blacklist update as --blacklist is not set.')
            fp = f'{args.outpath}/{os.path.basename(args.blacklist)}'
            blacklist.to_csv(fp, index=False)

            # Write QC tsv
            fp = f'{args.outpath}/exprs_qc.tsv'
            logging.info(f'Writing {fp}')
            qc_df.to_csv(fp, sep='\t', index=False)

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
