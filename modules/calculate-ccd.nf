process CALCULATE_CCD {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"
    conda "/mnt/storageBig8/work/micoli/miniconda3/envs/alpaca"

    input:
    tuple val(pat), path(alpaca_output)

    output:
    tuple val(pat), path("clone_copy_number_diversity_scores.csv")

    script:
    """
    alpaca ccd \\
        --alpaca_output_path ${alpaca_output} \\
        --output_directory .
    """
}
