"""PACT FIO pipeline
"""
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='fio',  # Required
    version='0.1',  # Required
    description='PACT FIO pipeline', 
    url='https://github.com/pact-pharma/bioinfo-fio',
    author='PACT Pharma', 
    author_email='no-reply@pactpharma.com', 
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Python :: 3.9',
    ],
    package_dir={'': 'src'},  # Optional
    packages=find_packages(where = 'src'),  # Required
    python_requires='>=3.9, <4',
    entry_points={  # Optional
        'console_scripts': [
            'process_expression.py=process_exprs.bin.process_expression:main',
            'batch_correct.py=process_exprs.bin.batch_correct:main',
            'counts2mat.py=process_exprs.bin.counts2mat:main',
            'timestamp_exprs.py=process_exprs.bin.timestamp_exprs:main',
            'update_exprs.py=process_exprs.bin.update_exprs:main',
            'tme_report.py=tme.bin.tme_report:main',
            'stage_input_tme.py=tme.bin.stage_input_tme:main',
            'prepare_ref_data.py=prepare_ref_data:main',
        ],
    },
)
