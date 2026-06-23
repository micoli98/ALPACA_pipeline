#!/usr/bin/env nextflow

/*
Description of the pipeline

Input:
- sample_info: data frame with columns "patient" and "sample"
- pubDir: path to the output directory
*/

nextflow.enable.dsl=2

// Include modules
include { HARMONIZE } from "./modules/harmonize.nf"
include { GET_PATIENT_INFO } from "./modules/get-patient-info.nf"

workflow {
    log.info """\

        ALPACA PIPELINE     
        ===================================
        patient_list        : ${params.sample_info}
        pubdir              : ${params.pubDir}
        """
        .stripIndent()
    
    // Input channel
    sample_info = channel.fromPath(params.sample_info)
    hfun = channel.fromPath(params.harmonization_function).first()
    input_assets = channel.value([file(params.seg), 
                                file(params.pp)])

    // Prepare patient input channel
    sample_input_ch = sample_info
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.sample, row.patient)}
        .groupTuple(by:1)

    // Extract patient-specific segments
    patient_info_ch = GET_PATIENT_INFO(sample_input_ch, input_assets)

    // Rejoin with batch info
    // patient_input_conv_ch = converted_cf_ch
    //     .join(patient_input_ch)
    //     .map{ patient, file_cf, segs, pp, model, file_cf_old, file_tree, snp_files -> 
    //         tuple(patient, segs, pp, model, file_cf, file_tree, snp_files)
    //     }

    // rp_input_ch = HARMONIZE(hfun, patient_input_conv_ch)
    HARMONIZE(hfun, patient_info_ch)

}
