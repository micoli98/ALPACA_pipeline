#!/usr/bin/env nextflow

/*
Description of the pipeline

Runs the ALPACA (clonal copy-number evolution, per patient) and SVclone
(structural variant CCF clustering, per sample) subworkflows over the same
cohort, and joins their outputs by patient into one terminal channel for
downstream consumption.

Input:
- sample_info: data frame with columns "patient" and "sample"
- outdir: path to the output directory
*/

nextflow.enable.dsl=2

include { ALPACA_SUBWORKFLOW } from "./subworkflows/alpaca.nf"
include { SVCLONE_SUBWORKFLOW } from "./subworkflows/svclone.nf"

workflow {
    log.info """\

        ALPACA + SVCLONE PIPELINE
        ===================================
        patient_list        : ${params.sample_info}
        outdir              : ${params.outdir}
        """
        .stripIndent()

    // Shared input: same cohort for both subworkflows
    sample_info_file = channel.fromPath(params.sample_info).first()

    // ----- ALPACA inputs -----
    batch_info_file = channel.fromPath(params.batch_info).first()
    input_assets = channel.value([file(params.seg),
                                file(params.pp)])

    sample_input_ch = sample_info_file
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.sample, row.patient)}

    batch_info_ch = batch_info_file
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.patient, row.batch, row.model, row.path_to_trees, row.file_cf, row.file_tree)}

    // ----- SVclone inputs -----
    coverage_info = channel.fromPath(params.coverage_info).first()

    // ----- Run subworkflows -----
    alpaca_ch = ALPACA_SUBWORKFLOW(
        sample_input_ch,
        batch_info_ch,
        input_assets
    )

    svclone_ch = SVCLONE_SUBWORKFLOW(
        sample_info_file,
        batch_info_file,
        coverage_info
    )

    // ----- Combine outputs by patient -----
    combined_ch = svclone_ch
        .map{ patient, sample, files -> tuple(patient, sample, files) }
        .groupTuple(by: 0)
        .join(alpaca_ch, by: 0)
}
