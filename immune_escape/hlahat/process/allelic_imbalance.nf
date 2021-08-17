/*
    Calculate allelic imbalance
*/
include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process ALLELIC_IMBALANCE {
 //   label 'process_medium'
    tag "${meta.specimen_id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    //conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(flagstat_normal_dna), path(pileups_normal_dna)
    tuple val(meta), path(flagstat_tumor_dna), path(pileups_tumor_dna)
    tuple val(meta), path(flagstat_tumor_rna), path(pileups_tumor_rna)

    output:
    path "*.tsv"
    path ".command*"

    script:
    def software = getSoftwareName(task.process)

    prefix_dna = meta.specimen_id + '_tumor_dna'
    prefix_rna = meta.specimen_id + '_tumor_rna'
    """
    Rscript /code/generate_ld.R\
     --tumor_depth_cutoff=15\
     --tumor_pileups=${pileups_tumor_dna}\
     --tumor_flagstat=${flagstat_tumor_dna}\
     --normal_depth_cutoff=15\
     --normal_pileups=${pileups_normal_dna}\
     --normal_flagstat=${flagstat_normal_dna}\
     --out_prefix=${prefix_dna}

    Rscript /code/generate_ld.R\
     --tumor_depth_cutoff=15\
     --tumor_pileups=${pileups_tumor_rna}\
     --tumor_flagstat=${flagstat_tumor_rna}\
     --out_prefix=${prefix_rna}
    """
}