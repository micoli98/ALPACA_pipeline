process ALPACA {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"

    input:
    tuple val(pat), 
        path(ci_table), 
        path(alpaca_input), 
        path(tree_paths), 
        path(cp_table)

    output:
    tuple val(pat), 
        path("ALPACA_output_${pat}.csv"),
        path("cn_change_to_ancestor.csv")

    script:
    """
    # Create input directory structure
    mkdir -p ${pat}
    cp ${alpaca_input} ${pat}/ALPACA_input_table.csv
    cp ${ci_table} ${pat}/ci_table.csv
    cp ${cp_table} ${pat}/cp_table.csv
    cp ${tree_paths} ${pat}/tree_paths.json
    
    \$CONDA_PREFIX/bin/python /mnt/storageBig8/work/micoli/miniconda3/envs/alpaca/bin/alpaca run \\
        --input_tumour_directory ${pat} \\
        --output_directory ${params.pubDir}/${pat}

    # Copy outputs back to work directory for Nextflow to track
    cp ${params.pubDir}/${pat}/ALPACA_output_${pat}.csv .
    cp ${params.pubDir}/${pat}/cn_change_to_ancestor.csv .
    """
}