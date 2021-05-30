Data processing
============


HISAT2 and HISAT-genotype
--------------

HLA read extraction
--------------

.. code-block:: bash
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin

    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
      -1 ${fq1} -2 ${fq2} ${true='--is-rna' false='' is_rna}\
      --database-list hla


HLA typing
--------------

docker: kcampbel/rnaseq_methods:v3

    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin

    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_locus_v_KC.py --genotype-genome ${hisat_prefix} \
      -U ${fq1} ${true='--is-rna' false='' is_rna}\
      --base hla --locus-list ${locus} --output-base ${name}.${locus} --keep-low-abundance-alleles

Construct HLA custom reference
--------------

docker: kcampbel/hlahat_r:v1

    grep "ranked" ${sep=" " hla_report_files} > ${name}.hla_types.txt
    Rscript /code/generate_reference_files.R ${name} ${hlatypes} ${sep="," gen_msf_list} ${sep="," nuc_msf_list}

Variant detection
--------------


Quantifying allelic imbalance
--------------


Paired tumor-normal data
--------------


Tumor-only datasets
--------------