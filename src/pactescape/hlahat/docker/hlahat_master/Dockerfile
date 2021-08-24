FROM rocker/r-ubuntu:18.04

LABEL \
    description="Master image for HLA-HAT"

RUN apt-get update -y && apt-get install -y \
    build-essential \
    bzip2 \
    cmake \
    default-jdk \
    git \
    libnss-sss \
    libtbb2 \
    libtbb-dev \
    ncurses-dev \
    nodejs \
    python-dev \
    python-pip \
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev

ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2 && \
    cd /tmp/htslib-1.3.2 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar --bzip2 -xf samtools-1.3.1.tar.bz2 && \
    cd /tmp/samtools-1.3.1 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/samtools-1.3.1

ENV HISAT2_INSTALL_DIR=/opt/hisat2

RUN mkdir /opt/hisat2/ && \
    wget https://github.com/DaehwanKimLab/hisat2/archive/hisat2_v2.2.0_beta.zip && \
    unzip -d /opt/hisat2/ hisat2_v2.2.0_beta.zip && \
    cd /opt/hisat2/hisat2-hisat2_v2.2.0_beta/ && \
    make hisat2-align-s hisat2-build-s hisat2-inspect-s

ENV PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
ENV PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
ENV PATH=$PATH:/opt/samtools/bin

## Copy files
RUN mkdir -p /code/
COPY r_scripts /code/

RUN R -e "install.packages(c('dplyr','tidyr','data.table','seqinr','BiocManager'), repos = 'https://cloud.r-project.org/')"
RUN R -e "BiocManager::install('msa')"
RUN R -e "BiocManager::install('Biostrings')"

COPY scripts/hisatgenotype_typing_common_v_KC.py /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules/hisatgenotype_typing_common_v_KC.py
COPY scripts/hisatgenotype_typing_common_v_KC.pyc /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules/hisatgenotype_typing_common_v_KC.pyc
COPY scripts/hisatgenotype_extract_reads_v_KC.py /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py
COPY scripts/hisatgenotype_locus_v_KC.py /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_locus_v_KC.py
