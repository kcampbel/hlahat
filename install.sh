#!/usr/bin/env bash
set -x
#FIO
conda env create -f conda/fio.yml --force 

# PROCESS_EXPRS
conda env create -f conda/process_exprs.yml --force 

# TME
conda env create -f conda/tme.yml --force &&
set +eu && \
    . $(conda info --base)/etc/profile.d/conda.sh && \
    conda activate tme && \
    set -eu
Rscript -e 'install.packages("https://ftp.osuosl.org/pub/cran/src/contrib/randomForest_4.6-14.tar.gz", repos=NULL, type="source")'
Rscript -e 'devtools::install_github("Sage-Bionetworks/CMSclassifier", ref="ba02ec5c88194e34c2a2b6a0c81bc07c36c8e9fa", type="source", dependencies=FALSE)'

