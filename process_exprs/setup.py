"""PACT Process RNAseq expression matrices
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='process_exprs',  # Required
    version='0.1',  # Required
    description='PACT RNAseq compilation', 
    url='https://github.com/pact-pharma/bioinfo-fio/process_exprs',
    author='PACT Pharma', 
    author_email='no-reply@pactpharma.com', 
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Python :: 3.7',
    ],
    package_dir={'': 'src'},  # Optional
    packages=find_packages(where = 'src'),  # Required
    python_requires='>=3.6, <4',
    entry_points={  # Optional
        'console_scripts': [
            'batch_correct.py=process_exprs.bin.batch_correct:main',
            'counts2mat.py=process_exprs.bin.counts2mat:main',
            'timestamp_exprs.py=process_exprs.bin.timestamp_exprs:main',
            'update_exprs.py=process_exprs.bin.update_exprs:main',
        ],
    },
)
