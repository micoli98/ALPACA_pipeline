process POST_PROCESS {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"

    input:
    tuple val(pat),
        path(alpaca_input),
        path(tree_paths),
        path(cp_table),
        path(alpaca_output)

    output:
    tuple val(pat), path("*")

    script:
    """
    # Recode chromosome X as 23 in the segment column so alpaca ccd/wgd
    # accept the format (they require purely numeric chromosome names)
    python3 -c "
    import pandas as pd, sys
    df = pd.read_csv('${alpaca_output}')
    df['segment'] = df['segment'].str.replace(r'^X_', '23_', regex=True)
    df.to_csv('alpaca_output_recoded.csv', index=False)
    "

    alpaca ccd \\
        --alpaca_output_path alpaca_output_recoded.csv \\
        --output_directory .

    alpaca wgd \\
        --alpaca_output_path alpaca_output_recoded.csv \\
        --tree_path ${tree_paths} \\
        --output_directory .

    alpaca plot-tumour \\
        --input_directory . \\
        --output_directory . \\
        --alpaca_output_path ${alpaca_output} \\
        --plot_output_mode notebook \\
        --genome_build hg38
    """
}
