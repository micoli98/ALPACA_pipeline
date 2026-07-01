process CALCULATE_WGD {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"

    input:
    tuple val(pat), path(alpaca_output), path(tree_paths)

    output:
    tuple val(pat), path("wgd_ratio_scores.csv")

    script:
    """
    alpaca wgd \\
        --alpaca_output_path ${alpaca_output} \\
        --tree_path ${tree_paths} \\
        --output_directory .
    """
}
