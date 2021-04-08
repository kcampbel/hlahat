FROM r-base:4.0.0

RUN mkdir -p /code/

## copy files
COPY r_scripts /code/

## install R-packages
RUN R -e "install.packages(c('tidyverse','data.table','seqinr','BiocManager'), repos = 'https://cloud.r-project.org/')"
RUN R -e "BiocManager::install('msa')"
RUN R -e "BiocManager::install('Biostrings')"
