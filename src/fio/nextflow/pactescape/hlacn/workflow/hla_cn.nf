// Stage inputs

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input_hlacn
    // params.fasta, params.transcript_fasta, params.additional_fasta,
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input_hlacn) { ch_input_hlacn = file(params.input_hlacn) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//include { HISAT_GENOTYPE } from '../process/hisat_genotype'  addParams( params.modules["hisat_genotype"] )

/*
========================================================================================
    RUN WORKFLOW
========================================================================================
*/ 

include { HLATYPE_COMPARE } from '../process/hlatype_compare'  addParams( params.modules["hlatype_compare"] )

workflow HLA_CN { 
    Channel.from(ch_input_hlacn)
        .splitCsv(sep: ',', header: true)
        .map{ row-> tuple(
            row.specimen_id,
            file(row.epic_hlatypes),
            file(row.normal_dna_bam),
            file(row.tumor_dna_bam),
            file(row.snp_bed),
            file(row.snp_tsv),
            file(row.seqz_file),
            file(row.sequenzaModelRData),
        ) }
        //.view { it } 
        .groupTuple(by: [0])
        .set { ch_input_hlacn }

    take:
    hisatgt_hlatypes
    //ch_hisatgt_hlatypes = file(params.hisatgt_hlatypes)

    main:
    HLATYPE_COMPARE (
	ch_input_hlacn,
	hisatgt_hlatypes
	)
//   SNP_GENOTYPER (),
//   HLA_COPYNUMBER()
    
}
