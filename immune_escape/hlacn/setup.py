"""HLA-CN 
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='hlacn',  # Required
    version='0.1',  # Required
    description='PACT HLA Copy number tool',
    url='https://github.com/pact-pharma/bioinfo-fio/immune_escape/hlacn',
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
            'compare_hla_typing.py=hlacn.bin.compare_hla_typing.py:main',
            'hisatgenotype_to_tsv.py=hlacn.bin.hisatgenotype_to_tsv.py:main',
            'imgt_reference_creator.py=hlacn.bin.imgt_reference_creator.py:main',
            'imgt_snp_finder.py=hlacn.bin.imgt_snp_finder.py:main',
            'snp2allele.py=hlacn.bin.snp2allele.py:main',
        ],
    },
)
