#!/usr/bin/env python
""" Combat batch correction
"""
#    Expression matrix must have Ensembl gene ids in rows, samples in columns. A config yml is required with the following:
#    xena_tpm_raw, pact_tpm_raw, meta, hgnc, blacklist 

import sys
import pandas as pd
import patsy
import numpy as np
import yaml
import argparse
import logging
from combat import combat
from process import read_exprs

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('primary_site', nargs='+', help='Primary tissue site(s)')
    parser.add_argument('-c', '--config', required=True, default='batch_correct.yml', help = 'YAML config file')
    parser.add_argument('-o', '--outfile', required=True, help='Output expression matrix file')
    parser.add_argument('-b', '--batch_name', required=True, help='Batch column in metadata')
    parser.add_argument('-f', '--formula', help='Model covariate formula (e.g. ~ factor)')
    parser.add_argument('-g', '--gid', default='hgnc', choices=['hgnc', 'gene_id'], help='Gene identifier to return')
    #parser.add_argument('-b', '--blacklist', default='/home/csmith/git/bioinfo-tme/pact/blacklist.txt',
    #parser.add_argument('-x', '--xena_tpm_raw', help='Xena expression matrix')
    #parser.add_argument('-p', '--pact_tpm_raw', help='PACT expression matrix')
    #parser.add_argument('-m' ,'--meta', help='Metadata tsv')
    #    help = 'Sample blacklist')
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
    logging.basicConfig(format=formatter, level=logging.DEBUG)
    logging.info('Starting batch correction')

    args = parse_args()
#    args = parse_args(
#        [
#            '-c', '/home/csmith/git/bioinfo-tme/test/batch_correct.yml',
#            '-o', '/tmp/bc.parquet.gz',
#            '-b', '_study',
#            'Ovary'
#        ])

    cfg = yaml.safe_load(open(args.config, 'r'))
    # Metadata
    meta = pd.read_table(cfg['meta'], index_col='specimen_id')
    meta = meta[
        (meta.primary_site.isin(args.primary_site)) &
        (meta._sample_type != 'Cell Line')
    ]

    # Expression matrices
    ## Xena
    meta_x = meta[meta._study.isin(['TCGA', 'GTEX'])]
    sids = meta_x.index.to_list()

    exprs_xena = read_exprs(cfg['xena_tpm_raw'], index_col='sample', samples=sids)
    exprs_xena.index = exprs_xena.index.str.rsplit('.').str[0]

    ## PACT
    meta_pact = meta[meta._study == 'PACT']
    sids = meta_pact.index.to_list()
    exprs_pact = read_exprs(cfg['pact_tpm_raw'], index_col='gene_id', samples=sids)

    ## Harmonize
    exprs_pact = exprs_pact.loc[:, exprs_pact.columns.isin(meta.index)]
    meta_pact = meta_pact.loc[exprs_pact.columns]
    exprs_pact = np.log2(exprs_pact + 0.001)

    # Merge expression matrices
    exprs_m = pd.concat([exprs_xena, exprs_pact], axis=1).fillna(np.log2(0.001))
    if not exprs_m.index.str.contains('ENSG').all():
        raise Exception('Expression matrix does not contain ENSG gene identifiers')


    if cfg['blacklist']:
        blacklist = [x.strip() for x in open(cfg['blacklist'])]
        dropme = exprs_m.filter(items = blacklist, axis=1).columns.to_list()
        if dropme:
            logging.info(f'Blacklisting {dropme} from {cfg["blacklist"]}')
            exprs_m.drop(dropme, axis=1, inplace=True)

    ## HGNC
    if args.gid == 'hgnc':
        hgnc_df = pd.read_table(cfg['hgnc'], index_col='ensembl_gene_id')
        hgnc_df.drop(columns='_merge')
        exprs_m = exprs_m.rpow(2)\
          .merge(hgnc_df, left_index=True, right_index=True)
        exprs_m['gene'] = np.where(exprs_m.gene.isna(), exprs_m.index, exprs_m.gene)
        exprs_m = np.log2(exprs_m.groupby('gene').sum())

    # Merge metadata
    meta_m = pd.concat([meta_x, meta_pact])
    sample_table = meta_m.groupby(['primary_site', '_study', 'tumor_normal']).primary_site.value_counts()
    logging.info(f'Merged sample counts:\n{sample_table}')
    logging.info(f'Merged genes : {exprs_m.shape[0]} samples: {exprs_m.shape[1]}')

    # Combat
    out = combat_split(exprs=exprs_m, meta=meta_m, split_name='primary_site', batch_name=args.batch_name, formula=args.formula)

    sites = ''.join(args.primary_site)
    logging.info(f'Writing to {args.outfile}')
    if 'parquet' in args.outfile:
        out.to_parquet(args.outfile, compression='gzip')
    else:
        out.to_csv(args.outfile, sep='\t')
    logging.info(f'Batch correction finished.')

if __name__ == "__main__":
    main()