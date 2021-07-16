#!/usr/bin/env python
""" Combat batch correction.
Performs batch correction after splitting expression data by primary tissue site.
"""

import sys
import os
import pandas as pd
import patsy
import numpy as np
import yaml
import argparse
import logging
from datetime import datetime
from process_exprs.combat import combat
from process_exprs.fileio import read_exprs, write_exprs, get_extension

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('exprs', nargs='+', help='Gene expression matrices to compile')
    parser.add_argument('--batch_name', required=True, help='Batch column in metadata')
    parser.add_argument('--primary_site', required=True, nargs=1, help='Primary tissue site(s), comma-delimited')
    parser.add_argument('--metadata_tsv', required=True, help='Specimen metadata')
    parser.add_argument('--exprs_log', help='Update log filename')
    parser.add_argument('--formula', help='Model covariate formula (e.g. ~ factor)')
    parser.add_argument('--blacklist', help='Specimen blacklist')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of samples with the same name')
    parser.add_argument('-o', '--outfile', required=True, help='Output path')

    return parser.parse_args(args)
    
def run_combat(exprs, meta, batch_name:str, formula:str, count_min=np.log2(0.001)):
    """ Combat batch correction 
    
    Args:
        exprs(pandas): expression matrix, genes in rows and samples in columns
        meta(pandas}: metadata, sample names as index, column with batch_name
        batch_name(str): batch column in meta
        formula(str): formula for preserving group differences (e.g. "~ tumor_normal"), where tumor_normal
            is a column in meta
        count_min(float): minimum expression value for a gene to be retained, where every sample 
            must exceed this value to pass

    Returns:
        pandas 
    """
    batch = meta.loc[:, batch_name]
    samples = meta.index

    # Remove genes where all samples have < count_min
    exprs_s = exprs[samples]
    genes_pass = exprs_s.apply(lambda x: (x > count_min).any(), axis=1)
    exprs_w = exprs_s[genes_pass]

    logging.info(f'Retained {genes_pass.sum()} genes with counts greater than count_min={count_min}')
    mod = patsy.dmatrix(formula, meta, return_type="dataframe")
    edat = combat(exprs_w, batch, mod)
        
    # Combat error handling
    if (edat.isna().apply(sum) != 0).any():
        for count_min in [0.005, 0.05, 0.5]:
            genes_pass = exprs_w.apply(lambda x: (x > np.log2(count_min)).any(), axis=1)
            exprs_pass = exprs_w.loc[genes_pass]
            logging.warning(f'Combat failed. Removing genes where all samples have a count <' +
                f'{count_min}:{genes_pass.sum()}/{exprs_w.shape[0]} genes retained.')
            edat = combat(exprs_pass, batch, mod)
            if edat.isna().any().any():
                continue
        if edat.isna().any().any():
            raise Exception('Combat failed. Aborting...')
    return edat

def combat_split(exprs, meta, split_name:str, batch_name:str, formula:str = None, count_min = np.log2(0.001)):
    """ Combat batch correction split by split_name 
    Input expression matrices are expected to be log2(x + 0.001) transformed

    Args:
        exprs(pandas): expression matrix, genes in rows and samples in columns
        meta(pandas}: metadata, sample names as index, column with batch_name
        batch_name(str): batch column in meta
        split_name(str): column to split matrix on in meta
        formula(str): formula for preserving group differences (e.g. "~ tumor_normal"), where tumor_normal
            is a column in meta
        count_min(float): minimum expression value for a gene to be retained, where every sample 
            must exceed this value to pass

    Returns:
        pandas 
    """
    out = pd.DataFrame()
    for ii in meta[split_name].unique():
        meta_split = meta[meta[split_name] == ii]
        exprs_split = exprs.loc[:, meta_split.index]
        if not formula:
            if meta_split.tumor_normal.str.contains('TCGANormal').any():
                formula = '~ tumor_normal'
            else:
                formula = '~ 1'
        logging.info(f'Combat split_name:{ii} batch_name:{batch_name} formula:{formula}')
        edat_combat = run_combat(exprs_split, meta_split, batch_name, formula)
        out = pd.concat([out, edat_combat], axis=1).fillna(np.log2(0.001))
    return(out)

def main():
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.INFO)
    logging.info(f'Starting {os.path.basename(__file__)}')
    now = datetime.now()
    timestamp = now.strftime("%Y%m%dT%H%M%S")

    args = parse_args()
#    args = parse_args(
#        [
#            './test/eset_xena_hgnc_proteincoding.tsv',
#            './test/eset_pact_hgnc_proteincoding.tsv',
##            './test/eset_xena_geneid_proteincoding.tsv',
##            './test/eset_pact_geneid_proteincoding.tsv',
#            '--batch_name', 'study',
#            '--metadata_tsv', './test/meta_eset.tsv',
#            '--exprs_log', './test/bc_exprs_log.tsv',
#            '--blacklist', './test/blacklist.txt',
#            '--primary_site', 'Ovary',
#            '-o', '/tmp/batch_correct/bc.tsv',
#        ])

    # Metadata
    primary_site = args.primary_site[0].split(',')
    meta = pd.read_table(args.metadata_tsv, index_col='specimen_id')
    meta = meta[
        (meta.primary_site.isin(args.primary_site)) &
        (meta.sample_type != 'Cell Line') &
        (meta.study.isin(['TCGA', 'GTEX', 'PACT']))
    ]
    if not meta.primary_site.isin(primary_site).all():
        diffs = set(primary_site) - set(meta.primary_site.unique())
        raise Exception(f'{",".join(diffs)} missing from metadata')

    # Expression matrices
    exprs = pd.Series(dtype=str)
    gene_names = pd.Series()
    for ii in args.exprs:
        tmp = read_exprs(ii, index_col=0)
        exprs = pd.concat([exprs, tmp], axis=1).fillna(0)
        if gene_names.empty:
            gene_names = tmp.index
        else:
            if not gene_names.isin(tmp.index.to_list()).any():
                raise Exception(f'Gene matrices do not share any genes')

    exprs_m = exprs.loc[:, exprs.columns.isin(meta.index)]
    exprs_m.index.name = tmp.index.name
    meta_m = meta.loc[exprs_m.columns]

    if args.blacklist:
        logging.info(f'Checking for blacklisted samples from {args.blacklist}')
        blacklist = [x.strip() for x in open(args.blacklist)]
        dropme = exprs_m.filter(items = blacklist, axis=1).columns.to_list()
        if dropme:
            logging.info(f'Blacklisting {dropme}')
            exprs_m.drop(dropme, axis=1, inplace=True)
            meta_m.drop(dropme, axis=0, inplace=True)

    # Merge metadata
    sample_table = meta_m.groupby(['primary_site', 'study', 'tumor_normal']).primary_site.value_counts()
    logging.info(f'Merged sample counts:\n{sample_table}')
    logging.info(f'Merged genes : {exprs_m.shape[0]} samples: {exprs_m.shape[1]}')

    # Combat
    exprs_m = np.log2(exprs_m + 0.001)
    out = combat_split(exprs_m, meta=meta_m, split_name='primary_site', batch_name=args.batch_name, formula=args.formula)

    dn = os.path.dirname(args.outfile)
    if not os.path.exists(dn):
        os.makedirs(dn, exist_ok=True)
    fp = args.outfile
    logging.info(f'Writing {fp}')
    write_exprs(exprs_m, fp, compression='gzip')

    # Log file
    if args.exprs_log:
        exprs_log = pd.read_csv(args.exprs_log, sep='\t')
        row = {
            'date': timestamp,
            'method': 'combat',
            'exprs_path': args.outfile
        }
        tmp = pd.DataFrame([row])
        exprs_log_new = pd.concat([exprs_log, tmp], ignore_index=True)
        fp = f'{dn}/{os.path.basename(args.exprs_log)}'
        logging.info(f'Writing {fp}')
        exprs_log_new.to_csv(fp, sep='\t', index=False)

    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
