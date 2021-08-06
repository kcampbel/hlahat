#!/usr/bin/env bash
set -x
#FIO
conda env create -f conda/fio.yml --force &&
#curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" &&
#unzip awscliv2.zip &&
#sudo ./aws/install

# PROCESS_EXPRS
conda env create -f conda/process_exprs.yml --force 

# TME
conda env create -f conda/tme.yml --force &&
set +eu && \
    . $(conda info --base)/etc/profile.d/conda.sh && \
    conda activate tme && \
    set -eu
Rscript -e 'devtools::install_github("Sage-Bionetworks/CMSclassifier", ref="ba02ec5c88194e34c2a2b6a0c81bc07c36c8e9fa", type="source")'

