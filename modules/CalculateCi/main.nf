process CALCULATE_CI {
    tag "${pat}"
    publishDir "${params.outdir}/${pat}/alpaca", mode: "copy"
    conda params.alpaca_env

    input:
    tuple val(pat),
        path(rp_segs),
        path(rp_snps),
        path(rp_pp),
        path(cp_table)

    output:
    tuple val(pat),
        path("ci_table.csv"),
        path("ALPACA_input_table.csv")

    script:
    // ALPACA's refphase converter: bootstrap CIs from phased BAF + COBALT logR.
    // Threshold 0 keeps every segment (needed for comparison with Purple); segments with
    // <= 10 het SNPs or no usable logR get a fixed CI span of 0.5.
    """
    python ${params.convert_refphase} \\
        --tumour_id ${pat} \\
        --output_dir . \\
        --copy_number_tool refphase \\
        --refphase_segments ${rp_segs} \\
        --refphase_snps ${rp_snps} \\
        --refphase_purity_ploidy ${rp_pp} \\
        --conipher_cp_table ${cp_table} \\
        --heterozygous_SNPs_threshold 0
    """
}