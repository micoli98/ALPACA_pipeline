process PLOT_TUMOUR {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}/plots", mode: "copy"
    conda "/mnt/storageBig8/work/micoli/miniconda3/envs/alpaca"

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
    mkdir -p input_${pat}
    cp ${alpaca_input} input_${pat}/ALPACA_input_table.csv
    cp ${ci_table} input_${pat}/ci_table.csv
    cp ${cp_table} input_${pat}/cp_table.csv
    cp ${tree_paths} input_${pat}/tree_paths.json

    alpaca plot-tumour \\
        --input_directory input_${pat} \\
        --output_directory plots_${pat} \\
        --alpaca_output_path ${alpaca_output} \\
        --plot_output_mode pdf \\
        --genome_build hg19
    """
}
