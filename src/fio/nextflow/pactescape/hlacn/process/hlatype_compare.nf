/*
    Extract HLA reads
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process HLATYPE_COMPARE {
 //   label 'process_medium'
    tag "${meta}"
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'merged_fastq', meta:meta, publish_by_meta:['id']) }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(epic_hlatypes), path(normal_dna_bam), path(tumor_dna_bam), path(snp_bed), path(snp_tsv), path(seqz_file), path(sequenzaModelRData)
    tuple val(hlahat_meta), path(hisatgt_hlatypes)
        
    output:
    path("*hlatypes_merged.txt"), emit: hla_types

    script:
    hlatypes_merged = meta + "_hlatypes_merged.txt"
    """
    hisatgenotype_to_tsv.py ${hisatgt_hlatypes} -o hlatypes.report.tsv
    compare_hla_typing.py ${epic_hlatypes} hlatypes.report.tsv -o ${hlatypes_merged}
    """
}
