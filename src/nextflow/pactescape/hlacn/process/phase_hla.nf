/*
    HLA Allele SNP Phasing Tool
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process PHASE_HLA {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), val(data), file(hlatypes)
    file(genome_fasta)

    output:
    tuple val(meta), val(data), file('*alleles.tsv'), file('*snps.tsv'), emit: tsvs
    path 'prep'
    path '.command*'

    script:
    snps = meta.id + '_snps.tsv'
    """
    imgt_snp_finder.py ${hlatypes} -o ${snps}
    snp2allele.py\
      -i ${meta.specimen_id}\
      -n ${data.normal_bam}\
      -t ${data.tumor_bam}\
      -p ${snps}\
      -m ${data.sequenzaModelRData}\
      -g ${genome_fasta}\
      -s ${params.sequenza_solution}\
      -o ./
    """
}
