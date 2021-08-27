import os
from importlib.resources import files
import re
import logging
import urllib.request
import pandas as pd
import time
import re
from commonLib.lib.search import locate

def check_paths_exist(paths):
    exists = dict()
    for ii in paths:
        if os.path.exists(ii):
            exists[ii] = True
        else:
            exists[ii] = False
    return exists

def package_file_path(package, filename):
    filepath = files(package).joinpath(filename)
    if filepath.exists():
        return str(filepath.absolute())
    else:
        return None

def get_file(src, dest, attempts=3, force=False):
    if not '//' in src:
        src = 'file://' + src
    fn = f'{dest}/{os.path.basename(src)}'
    logging.info(f'Copying {src} to {fn}')
    if not force and os.path.exists(fn):
        logging.info(f'{fn} found. Set force in {__name__} to overwrite.')
    else:
        #filename, msg = urllib.request.urlretrieve(src, fn)
        ii = 1
        while ii < attempts:
            try:
                filename, msg = urllib.request.urlretrieve(src, fn)
                break
            except:
                ii +=1
                logging.info(f'Retrying for the {ii}/{attempts} time')
        urllib.request.urlcleanup()
    return fn
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
    import s3fs
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
    import s3fs
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

def find_file(input_folder:str, filename:str, pattern:str = None, n:int = 1):
    """ Find a file
        Args:
            input_folder(str): folder to search
            filename(str): file name
            pattern(str): second pass pattern to search for within file name
            n(int): number of files with fn and pattern allowed
        Returns:
            list files matching filenames with pattern
    """
    hits = list(locate(filename, input_folder))
    if pattern:
        hits = [x for x in hits if re.search(pattern, x)]
    if len(hits) != n:
        raise ValueError(f'{n} file(s) not found searching {input_folder} for {filename}, ' + 
            f'pattern {pattern}\nhits: {hits}')
    else:
        return hits[0]

