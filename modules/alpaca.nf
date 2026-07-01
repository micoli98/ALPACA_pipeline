process ALPACA {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"
    conda "/mnt/storageBig8/work/micoli/miniconda3/envs/alpaca"

    input:
    tuple val(pat),
        path(ci_table),
        path(alpaca_input),
        path(tree_paths),
        path(cp_table)

    output:
    tuple val(pat),
        path("ALPACA_output_${pat}.csv"),
        path("cn_change_to_ancestor.csv"),
        path("*_plots.ipynb")

    script:
    """
    alpaca run \\
        --input_tumour_directory . \\
        --output_directory \$PWD \\
        --solver gurobi \\
        --genome_build hg38 \\
        --extra_columns complexity CI_score D_score
    """
}