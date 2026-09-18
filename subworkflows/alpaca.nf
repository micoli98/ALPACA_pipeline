// ALPACA subworkflow: Purple CNV + RefPhase phasing + ALPACA optimizer,
// producing clonal copy-number profiles along a phylogenetic tree (per patient).

include { NAME_CONVERSION } from "../modules/NameConversion/main.nf"
include { HARMONIZE } from "../modules/Harmonize/main.nf"
include { CLONAL_INFO } from "../modules/ClonalInfo/main.nf"
include { REFPHASE } from "../modules/Refphase/main.nf"
include { CALCULATE_CI } from "../modules/CalculateCi/main.nf"
include { ALPACA } from "../modules/Alpaca/main.nf"
include { GET_STATS } from "../modules/GetStats/main.nf"
include { POST_PROCESS } from "../modules/PostProcess/main.nf"

workflow ALPACA_SUBWORKFLOW {
    take:
    sample_info_ch   // tuple(sample, patient)
    batch_info_ch    // tuple(patient, batch, model, path_to_trees, file_cf, file_tree)
    input_assets     // value([seg, pp])

    main:
    // Reformat segmentation, snps and purity ploidy input per patient to adhere RefPhase requirements
    patient_input_ch = sample_info_ch
        .map{ sample, patient ->
            def snp_files = file("${params.snp_path}/${patient}/${sample}.amber.baf.tsv{,.gz}")

            // file() with glob returns a list, get the first (and should be only) match
            def snp_file = snp_files instanceof List ? snp_files[0] : snp_files

            // Check if file exists
            if (!snp_file.exists()) {
                error "SNP file not found for ${sample} (patient ${patient})"
            }

            tuple(patient, snp_file)
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
    converted_cf_ch = NAME_CONVERSION(cf_ch, input_assets)

    // Rejoin with batch info
    patient_input_conv_ch = converted_cf_ch
        .join(patient_input_ch)
        .map{ patient, file_cf, segs, pp, model, file_cf_old, file_tree, snp_files ->
            tuple(patient, segs, pp, model, file_cf, file_tree, snp_files)
        }

    rp_input_ch = HARMONIZE(patient_input_conv_ch)

    // Reformat clonal information per patient to adhere to ALPACA requirements
    cl_info_ch = CLONAL_INFO(patient_input_conv_ch)

    // Run RefPhase to phase both copy number and snp data
    rp_ch = REFPHASE(rp_input_ch)

    // Confidence interval calculation for copy number segments
    ci_ch = CALCULATE_CI(rp_ch)

    // Run ALPACA
    alpaca_input_ch = ci_ch
        .join(cl_info_ch)

    alpaca_ch = ALPACA(alpaca_input_ch)

    joined_ch = alpaca_ch
        .join(rp_ch)
        .join(cl_info_ch)

    // Generate statistics and plots
    analysis_ch = joined_ch
        .map{ pat, alpaca_out, ancestor_out, ipynb, rp_segs, rp_snp, rp_pp, tree, cp ->
            tuple(pat, alpaca_out, rp_segs, cp)
        }

    GET_STATS(analysis_ch)

    // v0.3.1 post-processing
    POST_PROCESS(
        alpaca_input_ch
            .join(alpaca_ch)
            .map{ pat, ci_table, alpaca_input, tree_paths, cp_table, alpaca_out, ancestor_out, ipynb ->
                tuple(pat, alpaca_input, tree_paths, cp_table, alpaca_out)
            }
    )

    emit:
    joined_ch.map{ pat, alpaca_out, ancestor_out, ipynb, rp_segs, rp_snp, rp_pp, tree, cp ->
        tuple(pat, alpaca_out, ancestor_out, ipynb, rp_segs, cp)
    }
}
