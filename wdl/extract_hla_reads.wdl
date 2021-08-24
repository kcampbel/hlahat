workflow GetHlaReads {
  # For Agilent XT HS2 protocols
  Boolean? agilent_xths2
  Boolean agilent_xths2_default = select_first([agilent_xths2, false])
  File? trimmer_jar

  # Reference genome
  String hisat_prefix
  File hisat_index_file
  Array[String] hisat_index = read_lines(hisat_index_file)

  # Data input
  String name
  File? in_fq1 = "NA"
  File? in_fq2 = "NA"
  File? in_bam
  File? in_fqList
  Boolean? is_rna = false # Assumes DNA

  # If bam
  if (defined(in_bam) && in_bam != "NA" && (in_fq1 == "NA" || !defined(in_fq1) ) ) {
  	call BamtoFastq {
    	input:
        	inbam = in_bam,
          name = name
    }
  }

  if (agilent_xths2_default && defined(in_fq1)) {
      call AgentTrimFastqs as AgentTrimFastqs_fastqInput {
        input:
          fastq_1 = in_fq1,
          fastq_2 = in_fq2,
          trimmer_jar = trimmer_jar
      }
  }

  # Pick if bam or paired fastqs
  File? fq_f1 = select_first([BamtoFastq.fq1, AgentTrimFastqs_fastqInput.trimmed_fq1, in_fq1, "NA"])
  File? fq_f2 = select_first([BamtoFastq.fq2, AgentTrimFastqs_fastqInput.trimmed_fq2, in_fq2, "NA"])

  if ((defined(fq_f1) && fq_f1 != "NA")) {
  	if ((defined(fq_f2) && fq_f2 != "NA")) {
      call ExtractHlaReads_paired as ExtractHlaReads_Fq_paired {
        input:
          fq1 = fq_f1,
          fq2 = fq_f2,
          hisat_index = hisat_index,
          hisat_prefix = hisat_prefix,
          is_rna = is_rna
      }
    }
    if ((!defined(fq_f2) || fq_f2 == "NA")) {
      call ExtractHlaReads_single as ExtractHlaReads_Fq_single {
        input:
          fq1 = fq_f1,
          hisat_index = hisat_index,
          hisat_prefix = hisat_prefix,
          is_rna = is_rna
      }
    }
  }

  # If fastqList, scatter read extraction
  if (defined(in_fqList) && in_fqList != "NA") {
    Array[Array[String]] fqList = read_tsv(in_fqList)

    scatter (fq in fqList) {
      if (agilent_xths2_default) {
        call AgentTrimFastqs as AgentTrimFastqs_fastqListInput {
          input:
            fastq_1 = fq[0],
            fastq_2 = fq[1],
            trimmer_jar = trimmer_jar
        }
      }
      File? fqList_f1 = select_first([AgentTrimFastqs_fastqListInput.trimmed_fq1, fq[0], "NA"])
      File? fqList_f2 = select_first([AgentTrimFastqs_fastqListInput.trimmed_fq2, fq[1], "NA"])

      if ((defined(fqList_f2) && fqList_f2 != "NA")) {
        call ExtractHlaReads_paired as ExtractHlaReads_Fqlist_paired {
          input:
            fq1 = fqList_f1,
            fq2 = fqList_f2,
            hisat_index = hisat_index,
            hisat_prefix = hisat_prefix,
            is_rna = is_rna
        }
      }

      if ((!defined(fqList_f2) || fqList_f2 == "NA")) {
        call ExtractHlaReads_single as ExtractHlaReads_Fqlist_single {
          input:
            fq1 = fqList_f1,
            hisat_index = hisat_index,
            hisat_prefix = hisat_prefix,
            is_rna = is_rna
        }
      }
    }

    if (defined(ExtractHlaReads_Fqlist_paired.hla_fq1)) {
      call MergeFastq as MergeFastq_Fq1_paired {
        input:
          fqs = ExtractHlaReads_Fqlist_paired.hla_fq1
      }
    }
    if (defined(ExtractHlaReads_Fqlist_single.hla_fq1)) {
      call MergeFastq as MergeFastq_Fq1_single {
        input:
          fqs = ExtractHlaReads_Fqlist_single.hla_fq1
      }
    }
    if (defined(ExtractHlaReads_Fqlist_paired.hla_fq2)) {
      call MergeFastq as MergeFastq_Fq2 {
        input:
          fqs = ExtractHlaReads_Fqlist_paired.hla_fq2
      }
    }
  }

  File final_fq1 = select_first([MergeFastq_Fq1_paired.mergedfq, MergeFastq_Fq1_single.mergedfq, ExtractHlaReads_Fq_paired.hla_fq1, ExtractHlaReads_Fq_single.hla_fq1])
  File? final_fq2 = select_first([MergeFastq_Fq2.mergedfq, ExtractHlaReads_Fq_paired.hla_fq2, "NA"])

  output {
    File hla_fq1 = final_fq1
    File? hla_fq2 = final_fq2
  }
}

task AgentTrimFastqs {
  File fastq_1
  File fastq_2
  File trimmer_jar

  command <<<
    /usr/bin/java -jar ${trimmer_jar} -fq1 ${fastq_1} -fq2 ${fastq_2} -v2 -out_loc $PWD
  >>>

  output {
      File trimmed_fq1 = glob("*R1*")[0]
      File trimmed_fq2 = glob("*R2*")[0]
      File MBC_text = glob("*txt*")[0]
  }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    disks: "local-disk 200 SSD"
    memory: "16G"
    cpu: 1
  }
}

task BamtoFastq {
  File? inbam
  String name

  command <<<
    /usr/gitc/bamUtil/bin/bam bam2FastQ --in ${inbam} --gzip --outBase ./${name}
  >>>

  output {
    File fq1 = "${name}_1.fastq.gz"
    File fq2 = "${name}_2.fastq.gz"
  }

  runtime {
    memory: "48G"
    cpu: 16
    disks: "local-disk 100 SSD"
    docker: "pici/genomics"
  }
}

task MergeFastq {
  Array[File?]+ fqs

  command <<<
    cat ${sep=" " fqs} > merged.fq.gz
  >>>

  output {
    File mergedfq = "merged.fq.gz"
  }

  runtime {
    memory: "8G"
    cpu: 1
    disks: "local-disk 50 SSD"
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
  }
}

task ExtractHlaReads_single {
  File? fq1
  Array[File]+ hisat_index
  String hisat_prefix
  Boolean is_rna

  command <<<
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin

    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
      -U ${fq1} ${true='--is-rna' false='' is_rna}\
      --database-list hla
  >>>

  output {
    File hla_fq1 = glob("*hla.extracted*.fq.gz")[0]
  }

  runtime {
    memory: "16G"
    cpu: 1
    disks: "local-disk 100 SSD"
    docker: "kcampbel/rnaseq_methods:v3"
  }
}

task ExtractHlaReads_paired {
  File? fq1
  File? fq2
  Array[File]+ hisat_index
  String hisat_prefix
  Boolean is_rna

  command <<<
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    export PATH=$PATH:/opt/samtools/bin

    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix} \
      -1 ${fq1} -2 ${fq2} ${true='--is-rna' false='' is_rna}\
      --database-list hla
  >>>

  output {
    File hla_fq1 = glob("*hla.extracted*.fq.gz")[0]
    File? hla_fq2 = glob("*hla.extracted.2.fq.gz")[0]
  }

  runtime {
    memory: "16G"
    cpu: 1
    disks: "local-disk 100 SSD"
    docker: "kcampbel/rnaseq_methods:v3"
  }
}
