#!/usr/bin/env python
""" Simple yaml to exported yaml parser

Takes a simple yaml without versions and compares it to an exported yaml with versions and returns
the matches
"""
import argparse
import re
import subprocess as sb
import yaml

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('simple', help='simple yaml file')
    parser.add_argument('--export', help='Exported yaml file')
    
    return parser.parse_args(args)

def conda_export(environment):
    cmd = ['conda', 'env', 'export', '-n', environment, '--no-builds']
    job = sb.run(cmd, check=True, capture_output=True, text=True)
    return(yaml.safe_load(job.stdout))

def main():
    args = parse_args()
#    args = parse_args(
#        [
#            'conda/fio.yml',
#        ]
#    )
    src = yaml.safe_load(open(args.simple))
    if args.export:
        export = yaml.safe_load(open(args.export))
    else:
        export = conda_export(src['name'])

    src_deps = src['dependencies']
    ex_deps = export['dependencies']

    hits = list()
    for ii in src_deps:
        if type(ii) == dict:
            continue
        spkg = ii.split('=')[0]
        hit = [ x for x in ex_deps if spkg == x.split('=')[0]]
        if not hit:
            print(f'{ii} not found in {args.export}')
        else:
            hits.extend(hit)
    [ print(f'- {x}') for x in hits ]

if __name__ == '__main__':
    main()
    


