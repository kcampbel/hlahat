""" Tumor microenviroment reporting
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='tme',  # Required
    version='0.1',  # Required
    description='Tumor microenviroment reporting', 
    url='https://github.com/pact-pharma/bioinfo-fio/tme',
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
            'tme_report.py=tme.bin.tme_report:main',
            'stage_input_tme.py=tme.bin.stage_input_tme:main',
            'prepare_ref_data.py=prepare_ref_data:main',
        ],
    },
)
