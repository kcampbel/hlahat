import pandas as pd
import os
import logging
import s3fs
import re

def read_exprs(exprs_f:str, index_col:str, samples:list = None):
    if 'parquet' in exprs_f:
        if re.match('s3://', exprs_f):
            s3 = s3fs.S3FileSystem()
            infile = exprs_f.strip('s3://')
            with s3.open(infile, 'rb') as f:
                exprs = pd.read_parquet(f)
        else:
            exprs = pd.read_parquet(exprs_f)
    else:
        exprs = pd.read_table(exprs_f, index_col=index_col)
    if samples:
        exprs = exprs.filter(items=samples, axis=1)
    return(exprs)

def counts2mat(counts:list, metric:str, gene_name:str, eset, force:bool):
    if eset.index.name != gene_name:
        raise Exception(f'Eset gene column name {eset.index.name} does not match argument --gene_name {gene_name}. Aborting...')
    # Add new samples
    added, skipped = list(), list()
    for ii in counts:
        sampleId = os.path.basename(ii).split('_tumor')[0].strip('RNA_')
        if eset.columns.isin([sampleId]).any():
            if force:
                logging.info(f'--force enabled. Dropping {sampleId} from eset')
                eset = eset.drop(sampleId, axis=1)
            else:
                logging.warning(f'{sampleId} already exists in eset. Set --force to overwrite.')
                skipped.append(ii)
                continue
        dat = pd.read_csv(ii, sep='\t')

        # Sum rows by gene name and remove genes with zero expression
        tmp = dat.groupby(gene_name)[metric].sum()
        tmp.name = sampleId
        tmp = tmp[tmp > 0]

        # Merge with eset matrix
        logging.info(f'Adding {sampleId} to eset with {eset.shape[1]} samples')
        eset = eset.merge(tmp, how='outer', left_index=True, right_index=True).fillna(0)
        added.append(ii)

    logging.info(f'{metric} from {len(added)} samples added and {len(skipped)} skipped.')
    if added:
        return eset
    else:
        return None

def get_extension(filename):
    basename = os.path.basename(filename)  # os independent
    split = basename.split('.')
    prefix = split[0]
    ext = '.'.join(split[1:])
    return [prefix, '.' + ext] if ext else None
