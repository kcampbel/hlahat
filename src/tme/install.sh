conda env create -f environment.yml --force &&
set +eu && \
    . $(conda info --base)/etc/profile.d/conda.sh && \
    conda activate tme && \
    set -eu
Rscript -e 'devtools::install_github("Sage-Bionetworks/CMSclassifier", ref="ba02ec5c88194e34c2a2b6a0c81bc07c36c8e9fa", type="source")'

