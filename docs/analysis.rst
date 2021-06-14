Data processing using HISAT2 and HISAT-genotype
================================================

HLA read extraction
--------------------
Extract HLA-associated reads from DNA sequencing data::

		export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
		export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
		export PATH=$PATH:/opt/samtools/bin
		/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
			-1 ${fq1} -2 ${fq2} --database-list hla

To extract HLA-associated reads from RNA sequencing data, just add the ``--is-rna`` flag to the command::

    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
      -1 ${fq1} -2 ${fq2} --is-rna --database-list hla

Reads may also be extracted from single reads, using the ``-U`` option instead of paired ``-1`` and ``-2`` reads::
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
    -U ${fq1} --is-rna --database-list hla

HLA typing
--------------
HLA typing can be applied to 26 genes, spanning 7Mb in chromosome 6p21.3, including classical and non-classical Class I/II HLA genes, non-expressed Class I HLA pseudogenes, ATP binding cassette transporter genes, Class I chain-related and Class I-like genes: HLA-A, HLA-B, HLA-C, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DPB2, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-E, HLA-F, HLA-G, HLA-H, HFE, HLA-K, HLA-L, MICA, MICB, TAP1, TAP2, and HLA-V. These are specified as a comma-delimited list in the following HLA typing command by the `--locus-list` option. Note: Do not include "HLA-" when specifying HLA genes (e.g. use ``--locus-list A`` instead of ``--locus-list HLA-A``). Typing generally takes the longest on the HLA-A locus, and this command can be scattered across each loci individually, followed by report aggregation and summary (See suggested workflow).

Docker Commands (kcampbel/rnaseq_methods:v3)::

		export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
		export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
		export PATH=$PATH:/opt/samtools/bin
		/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_locus_v_KC.py --genotype-genome ${hisat_prefix} \
			-1 ${fq1} -- ${fq2} --base hla --locus-list ${locus} --output-base ${name}.${locus} --keep-low-abundance-alleles

Options for single-read sequencing (``U``) and RNA sequencing (``--is-rna``) can also be used for HLA typing.

Documentation for HISAT-genotype suggests using the alleles ranked 1 or 2, from hisatgenotype_locus.py for each gene, and the alleles may be reported up to the 4th field of resolution, which describes genomic differences in alleles outside of the coding regions. However, WES may not have sufficient sequencing coverage and RNAseq data would not be appropriate for detecting this level of information. HLA-HAT outputs the *${id}.all_types.tsv* file, indicating the ranked alleles, by abundance, to include the most comprehensive output from HISAT-genotype.

*${id}.all_types.tsv* file is a tab-delimited file derived from the report outputted by hisatgenotype_locus.py:
.. csv-table::
  :widths: auto
  :align: center
  :header: "Field", "Type", "Description"

  "ranks", "Integer", "Gene rank of allele, based upon percent abundance of reads assigned to corresponding HLA type"
  "alleles", "String", "Full resolution of ranked allele identified by HISAT-genotype"
  "gene", "String", "HLA gene"
  "perc_abundance", "Float", "Relative abundance of reads corresponding to allele"

By default, all alleles are reduced to their fullest resolution or up the third field of resolution (e.g. A*02:89 would remain A*02:89, while A*03:01:01:01 is reduced to A*03:01:01). Then, alleles up to the third field of resolution are summarized by the maximum percent abundance across those that are shared. Any alleles with less than 5% abundance are removed, and then the remaining one or top two alleles (at the third field of resolution) are chosen as the HLA types.

*Optional*: To summarize at the second field of resolution, the flag --field_of_resolution 2 can be used. If this parameter is set, the top two alleles at the second field of resolution are chosen.

*Example*: If the following Class I alleles are ranked in the report from HISAT-genotype:
.. csv-table::
   :widths: auto
	 :align: center
	 :header: "ranks", "alleles", "gene", "perc_abundance"

   "1", "A*02:01:01:01", "A", "40.85"
   "2", "A*33:01:01", "A", "31.63"
   "3", "A*33:03:23", "A", "13.97"
   "4", "A*34:01:01", "A", "4.52"
   "5", "A*34:05", "A", "4.52"
   "6", "A*34:14", "A", "4.52"
   "1", "B*14:02:01:01", "B", "50.79"
   "2", "B*15:01:01:01", "B", "37.33"
   "3", "B*15:01:01:03", "B", "11.87"
   "1", "C*08:02:01:01", "C", "51.18"
   "2", "C*03:03:01:01", "C", "48.82"

First, alleles are summarized to the third field of resolution:
.. csv-table::
   :widths: auto
	 :align: center
	 :header: "ranks", "alleles", "gene", "perc_abundance"

   "1", "A*02:01:01", "A", "40.85"
   "2", "A*33:01:01", "A", "31.63"
   "3", "A*33:03:23", "A", "13.97"
   "4", "A*34:01:01", "A", "4.52"
   "5", "A*34:05", "A", "4.52"
   "6", "A*34:14", "A", "4.52"
   "1", "B*14:02:01", "B", "50.79"
   "2", "B*15:01:01", "B", "37.33"
   "3", "B*15:01:01", "B", "11.87"
   "1", "C*08:02:01", "C", "51.18"
   "2", "C*03:03:01", "C", "48.82"

Alleles are summarized by the maximum percent abundance corresponding to each unique allele at the third field of resolution:
.. csv-table::
   :widths: auto
	 :align: center
	 :header: "ranks", "alleles", "gene", "perc_abundance"

   "1", "A*02:01:01", "A", "40.85"
   "2", "A*33:01:01", "A", "31.63"
   "3", "A*33:03:23", "A", "13.97"
   "4", "A*34:01:01", "A", "4.52"
   "5", "A*34:05", "A", "4.52"
   "6", "A*34:14", "A", "4.52"
   "1", "B*14:02:01", "B", "50.79"
   "2", "B*15:01:01", "B", "37.33"
   "1", "C*08:02:01", "C", "51.18"
   "2", "C*03:03:01", "C", "48.82"

Alleles with less than 5% abundance are removed:
.. csv-table::
   :widths: auto
	 :align: center
	 :header: "ranks", "alleles", "gene", "perc_abundance"

   "1", "A*02:01:01", "A", "40.85"
   "2", "A*33:01:01", "A", "31.63"
   "3", "A*33:03:23", "A", "13.97"
   "1", "B*14:02:01", "B", "50.79"
   "2", "B*15:01:01", "B", "37.33"
   "1", "C*08:02:01", "C", "51.18"
   "2", "C*03:03:01", "C", "48.82"

Finally, the top 1-2 ranked alleles are identified as the patient HLA type:
.. csv-table::
   :widths: auto
	 :align: center
	 :header: "ranks", "alleles", "gene", "perc_abundance"

   "1", "A*02:01:01", "A", "40.85"
   "2", "A*33:01:01", "A", "31.63"
   "1", "B*14:02:01", "B", "50.79"
   "2", "B*15:01:01", "B", "37.33"
   "1", "C*08:02:01", "C", "51.18"
   "2", "C*03:03:01", "C", "48.82"

The final list of HLA types is summarized by *${id}.top_hlatypes.tsv*, a tab-delimited file containing the filtered allele calls:
.. csv-table::
  :widths: auto
  :align: center
  :header: "Field", "Type", "Description"

  "gene", "String", "HLA gene"
  "allele", "String", "Filtered allele call"



Constructing a custom HLA reference
====================================

docker: kcampbel/hlahat_r:v1

		grep "ranked" ${sep=" " hla_report_files} > ${name}.hla_types.txt
		Rscript /code/generate_reference_files.R ${name} ${hlatypes} ${sep="," gen_msf_list} ${sep="," nuc_msf_list}

Variant detection
------------------


Quantifying allelic imbalance
------------------------------


Paired tumor-normal data
-------------------------


Tumor-only datasets
---------------------
