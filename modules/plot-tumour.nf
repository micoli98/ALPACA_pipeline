process PLOT_TUMOUR {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}/plots", mode: "copy"

    input:
    tuple val(pat),
        path(ci_table),
        path(alpaca_input),
        path(tree_paths),
        path(cp_table),
        path(alpaca_output)

    output:
    tuple val(pat), path("plots_${pat}/*")

    script:
    """
    alpaca plot-tumour \\
        --input_directory . \\
        --output_directory plots_${pat} \\
        --alpaca_output_path ${alpaca_output} \\
        --plot_output_mode pdf \\
        --genome_build hg19
    """
}
