import os
from importlib.resources import files
import re
import logging
import urllib.request

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

