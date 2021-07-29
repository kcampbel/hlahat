""" TME Input stager
Writes a pipeline yaml from epic manifest and pipeline data
"""
import argparse
import os
import yaml
import logging
from importlib.resources import files
from tme import R, data
from tme.search import locate

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('manifest', help='EPIC pipeline sample manifest')
    parser.add_argument('input_folder', help='EPIC pipeline input folder')
    parser.add_argument('-o', '--outfile', default='tme.yml', help='Pipeline yml to write')
    
    return parser.parse_args(args)

def package_file_path(package, filename):
    filepath = files(package).joinpath(filename)
    if filepath.exists():
        return str(filepath.absolute())
    else:
        return None

def tme_config(specimen_id, input_folder):
    def _find_file(input_folder, fn):
        hits = list(locate(fn, input_folder))
        if len(hits) != 1:
            raise Exception(f'{specimen_id} has >1 input files:\n{hits}')
        else:
            return hits[0]

    # Pipeline datafiles and parameters
    config = {
    'goi_f': package_file_path(data, 'pact_goi.tsv'),
    'geneset_f': package_file_path(data, 'genesets.txt'),
    'hgnc_f': package_file_path(data, 'hgnc_biomart.tsv.gz'),
    'xcell_celltypes_f': package_file_path(data, 'xcell_types.txt'),
    'hotspots_f': package_file_path(data, 'chang2007_hotspots.tsv'),
    }
    params = yaml.safe_load(files(data).joinpath('tme_report_params.yml').read_text())
    config.update(params)
    
    # EPIC inputs
    # Sequenza copy number segments file
    fn = f'{specimen_id}_segments.txt'
    cn_f = _find_file(input_folder, fn)

    # Epipope annotator
    fn = f'{specimen_id}_oncotator_results_Annotated.tsv'
    vars_f = _find_file(input_folder, fn)

    config.update(
        {
            'cn_f': cn_f,
            'vars_f': vars_f
        }
    )
    return config

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
#    args = parse_args([
#        './test_data/PACT056_T_196454/manifest.yml',
#        './test_data',
#		'-o', '/tmp/PACT056_T_196454_tme.yml',
#    ])
    mf_d = yaml.safe_load(open(args.manifest))
    mf_d['manifest_f'] = os.path.abspath(args.manifest)
    specimen_id = mf_d['pipeline']['pact_id']

    config = tme_config(specimen_id, args.input_folder)
    config.update(mf_d)
    
    if not args.outfile:
        args.outfile = f'{specimen_id}_tme.yaml'
    logging.info(f'Writing {args.outfile}')
    stream = open(args.outfile, 'w')
    yaml.dump(config, stream)
    stream.close()

    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == "__main__":
    main()
