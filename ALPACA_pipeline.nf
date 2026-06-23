#!/usr/bin/env nextflow

/*
Description of the pipeline

Input:
- sample_info: data frame with columns "patient" and "sample"
- pubDir: path to the output directory
*/

nextflow.enable.dsl=2

// Include modules
include { NAME_CONVERSION} from "./modules/name-conversion.nf"
include { HARMONIZE} from "./modules/harmonize.nf"
include { CLONAL_INFO} from "./modules/clonal-info.nf"
include { REFPHASE} from "./modules/refphase.nf"
include { CALCULATE_CI} from "./modules/calculate-ci.nf"
include { ALPACA} from "./modules/alpaca.nf"
include { GET_STATS } from "./modules/get-stats.nf"

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
    conv_table = channel.fromPath(params.conversion_table).first()
    hfun = channel.fromPath(params.harmonization_function).first()
    tfun = channel.fromPath(params.tree_conversion_functions).first()
    cfun = channel.fromPath(params.ci_functions).first()
    afun = channel.fromPath(params.clones_comparison_functions).first()
    batch_info = channel.fromPath(params.batch_info)
    input_assets = channel.value([file(params.seg), 
                                file(params.pp)])

    // Prepare patient input channel
    sample_input_ch = sample_info
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.sample, row.patient)}

    batch_info_ch = batch_info
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.patient, row.batch, row.model, row.path_to_trees, row.file_cf, row.file_tree)}

    // Reformat segmentation, snps and purity ploidy input per patient to adhere RefPhase requirements
    patient_input_ch = sample_input_ch
        .map{ sample, patient -> 
            def snp_files = file("${params.snp_path}/${patient}/${sample}.amber.baf.tsv{,.gz}")
        
            // file() with glob returns a list, get the first (and should be only) match
            def snp_file = snp_files instanceof List ? snp_files[0] : snp_files
            
            // Check if file exists
            if (!snp_file.exists()) {
                error "SNP file not found for ${sample} (patient ${patient})"
            }
            
            tuple(patient, snp_file)

            // def snp_file = file("${params.snp_path}/${patient}/${sample}.amber.baf.tsv")
            // tuple(patient, snp_file)
        }
        .groupTuple(by:0)
        .join(batch_info_ch)
        .map{ patient, snp_files, batch, model, path_to_trees, file_cf, file_tree -> 
                tuple(patient, model, file_cf, file_tree, snp_files) 
            }

    cf_ch = patient_input_ch
        .map{ patient, model, file_cf, file_tree, snp_files -> 
            tuple(patient, file_cf)
        }

    // Name conversion step: evolution data have outdated sample names. To match with Purple output, convert them
    converted_cf_ch = NAME_CONVERSION(cf_ch, conv_table, input_assets)

    // Rejoin with batch info
    patient_input_conv_ch = converted_cf_ch
        .join(patient_input_ch)
        .map{ patient, file_cf, segs, pp, model, file_cf_old, file_tree, snp_files -> 
            tuple(patient, segs, pp, model, file_cf, file_tree, snp_files)
        }

    rp_input_ch = HARMONIZE(hfun, patient_input_conv_ch)
    
    // Reformat clonal information per patient to adhere to ALPACA requirements
    cl_info_ch = CLONAL_INFO(patient_input_conv_ch, tfun)
    
    // Run RefPhase to phase both copy number and snp data
    rp_ch = REFPHASE(rp_input_ch)

    // Confidence interval calculation for copy number segments
    ci_ch = CALCULATE_CI(rp_ch, cfun)

    // Run ALPACA
    alpaca_input_ch = ci_ch
        .join(cl_info_ch)
    
    clones_ch = ALPACA(alpaca_input_ch)

    // Generate statistics and plots
    analysis_ch = clones_ch
        .join(rp_ch)
        .join(cl_info_ch)
        .map{ pat, alpaca_out, ancestor_out, rp_segs, rp_snp, rp_pp, tree, cp ->
            tuple(pat, alpaca_out, rp_segs, cp)
        }

    GET_STATS(analysis_ch, afun)
}
