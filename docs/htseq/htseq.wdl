## Copyright Parker Institute for Cancer Immunotherapy, 2018
##
## # pre-processing_universal
##
## Preprocessing workflow for Bams and Fastqs
##
## This is meant to be able to read in a variety of different formats
## to prepare files for further analysis using the Broad best practices.
##
## ##Inputs
##  File bam_file
##  File reference_gtf
##
## ##Outputs
##  File analysis_ready_bam
##  File analysis_ready_bam_index
##  File analysis_ready_bam_md5
##
##
## Maintainer: Katie Campbell katiecampbell@mednet.ucla.edu
##
## Github: [https://github.com/ParkerICI/firecloud-methods/tree/master/methods/hisat](https://github.com/ParkerICI/firecloud-methods/tree/master/methods/htseq)
##
## Licensing :
## This script is released under the PICI Informatics License (GPL-3.0) Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script.

workflow HisatAlignmentWorkflow {
  # FILE INFO
  String base_file_name
  File aligned_bam
  String? strandness
  String? strandness_arg = if defined(strandness) then strandness else ""

  String? htseq_strandness_arg = if strandness_arg == "RF" || strandness_arg == "R" || strandness_arg == "reverse" then "reverse" else if strandness_arg == "FR" || strandness_arg == "F" || strandness_arg == "first" then "yes" else "no"


  # REFERENCE
  File reference_gtf

  call SortBam {
    input:
      base_file_name = base_file_name,
      bam_file = aligned_bam
  }

  call HtSeqCount {
    input:
      base_file_name = base_file_name,
      bam_file = SortBam.sorted_bam,
      gtf = reference_gtf,
      strandness = htseq_strandness_arg
  }

  output {
    File htseq_counts_file = HtSeqCount.counts_file
  }
}

task SortBam {
  String base_file_name
  File bam_file

  command <<<
    /usr/local/bin/samtools sort -n -o ${base_file_name}.namesort.bam ${bam_file}
  >>>

  output {
    File sorted_bam = "${base_file_name}.namesort.bam"
  }

  runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
      disks: "local-disk 500 SSD"
      memory: "8G"
      cpu: 2
  }
}

task HtSeqCount {
  String base_file_name
  File bam_file
  File gtf
  String strandness

  command <<<
    /usr/local/bin/htseq-count --format bam --order name --mode intersection-strict --stranded ${strandness} --minaqual 1 --type exon --idattr gene_id ${bam_file} ${gtf} > ${base_file_name}.HTSeqcounts.tsv
  >>>

  output {
    File counts_file = "${base_file_name}.HTSeqcounts.tsv"
  }

  runtime {
      docker: "jhart99/htseq"
      disks: "local-disk 100 SSD"
      memory: "16G"
      cpu: 2
  }
}
