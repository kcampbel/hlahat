""" TME Input stager
Writes a pipeline yaml from epic manifest and pipeline data
"""
import argparse
import os
import re
import yaml
import logging
from importlib.resources import files
from process_exprs import data
from commonLib.fileio import package_file_path
from commonLib.search import locate

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_folder', help='EPIC pipeline input folder')
    parser.add_argument('-o', '--outfile', default='process_exprs.yml', help='Pipeline yml to write')
    
    return parser.parse_args(args)

def pipeline_config(specimen_id:str, input_folder:str):
    def _find_file(input_folder, fn, pattern:str = None):
        hits = list(locate(fn, input_folder))
        if pattern:
            hits = [x for x in hits if re.search(pattern, x)]
        if len(hits) != 1:
            raise ValueError(f'{specimen_id} has != 1 input file:\n{fn} {hits}')
        else:
            return hits[0]

    # Pipeline datafiles and parameters
    config = {
    'tcga_gtex_map': package_file_path(data, 'tcga_gtex.tsv'),
    }
    params = yaml.safe_load(files(data).joinpath('pipeline.yml').read_text())
    config.update(params)
    
    # EPIC inputs
    # Counts for sample
    counts_f = f'{specimen_id}_genes.tsv'

    config.update(
        {
            'counts': counts_f,
        }
    )
    return config

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
    args = parse_args([
        './test_data/PACT056_T_196454/manifest.yml',
        './test_data',
		'-o', '/tmp/PACT056_T_196454_tme.yml',
    ])
    mf_d = yaml.safe_load(open(args.manifest))
    mf_d['manifest_f'] = os.path.abspath(args.manifest)
    specimen_id = mf_d['pipeline']['pact_id']

    config = pipeline_config(specimen_id, args.input_folder)
    config.update(mf_d)
    
    if not args.outfile:
        args.outfile = f'{specimen_id}_process_exprs.yaml'
    logging.info(f'Writing {args.outfile}')
    stream = open(args.outfile, 'w')
    yaml.dump(config, stream)
    stream.close()

    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
