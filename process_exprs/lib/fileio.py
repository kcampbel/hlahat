import os
import pandas as pd
import s3fs
import time
import re

def get_extension(filename:str):
    basename = os.path.basename(filename)  # os independent
    split = basename.split('.')
    prefix = split[0]
    ext = '.'.join(split[1:])
    return [prefix, '.' + ext] if ext else None

def file_time(filename:str, stamp:str='modified'):
    if stamp == 'modified':
        ftime = time.localtime(os.path.getmtime(filename))
    if stamp == 'created':
        ftime = time.localtime(os.path.getctime(filename))
    return time.strftime("%Y%m%dT%H%M%S", ftime)

def append_basename(filename:str, suffix:str):
    prefix, ext = get_extension(filename)
    new = f'{prefix}_{suffix}{ext}'
    orig = f'{prefix}{ext}'
    return(new)

def read_exprs(filename:str, index_col:str, samples:list = None):
    if 'parquet' in filename:
        if re.match('s3://', filename):
            s3 = s3fs.S3FileSystem()
            fn = filename.split('s3://')[-1]
            with s3.open(fn, 'rb') as f:
                exprs = pd.read_parquet(f)
        else:
            exprs = pd.read_parquet(filename)
    else:
        exprs = pd.read_table(filename, index_col=index_col)
    if samples:
        exprs = exprs.filter(items=samples, axis=1)
    return(exprs)

def write_exprs(df, filename:str, compression:str='gzip'):
    if 'parquet' in filename:
        if re.match('s3://', filename):
            s3 = s3fs.S3FileSystem()
            fn = filename.split('s3://')[-1]
            with s3.open(fn, 'wb') as f:
                exprs = df.to_parquet(f, compression=compression)
        else:
            exprs = df.to_parquet(filename, compression=compression)
    else:
        exprs = df.to_csv(filename, sep='\t')
    return(filename)




