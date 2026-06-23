process Prepare_RP_input {
    publishDir "$pubDir/${patient}", mode: "copy"

    input: 
    val input_rp
    val patient
    val pubDir
    val segmentation
    val snp_path
    val purity_ploidy
    val batch_info

    output:
    tuple val(patient), 
        path("${patient}-segments.tsv"), 
        path("${patient}-snps.tsv"), 
        path("${patient}-purity_ploidy.tsv"), emit: pat_rp_input

    script: 
    """
    Rscript $input_rp ${patient} $pubDir $segmentation $snp_path $purity_ploidy $batch_info
    """
}

process Get_clonal_info {
    publishDir "$pubDir/${patient}", mode: "copy"

    input: 
    val get_cf_tree
    val patient
    val pubDir
    val batch_info

    output:
    tuple val(patient), 
        path("tree_paths.json"), 
        path("cp_table.csv"), emit: pat_cf_tree

    script: 
    """
    Rscript $get_cf_tree ${patient} $pubDir $batch_info
    """
}

process Run_RP {
    cpus 4
    memory '8 GB'
    publishDir "$pubDir/${patient}", mode: "copy"

    input: 
    val run_rp
    tuple val(patient), path("${patient}-segments.tsv"), path("${patient}-snps.tsv"), path("${patient}-purity_ploidy.tsv")
    val pubDir
    val purity_ploidy

    output:
    tuple val(patient), 
        path("${patient}-refphase-segmentation.tsv"), 
        path("${patient}-refphase-phased-snps.tsv.gz"), 
        path("${patient}-refphase-sample-data-updated.tsv"), emit: pat_rp_output

    script: 
    """
    Rscript $run_rp ${patient} $pubDir "${patient}-purity_ploidy.tsv"
    """
}

process Get_CI {
    publishDir "$pubDir/${patient}", mode: "copy"

    input: 
    val get_ci
    tuple val(patient), path("${patient}-refphase-segmentation.tsv"), path("${patient}-refphase-phased-snps.tsv.gz"), path("${patient}-refphase-sample-data-updated.tsv")
    val pubDir

    output:
    tuple val(patient), path("ci_table.csv"), path("ALPACA_input_table.csv"), emit: pat_ci_output

    script: 
    """
    Rscript $get_ci ${patient} $pubDir
    """
}

process Run_ALPACA {
    cpus 8 
    publishDir "$pubDir/${patient}", mode: "copy"
    //conda 'alpaca'

    input:
    tuple val(patient), path("ci_table.csv"), path("ALPACA_input_table.csv"), path("tree_paths.json"), path("cp_table.csv")
    val pubDir

    output:
    tuple val(patient), path("*"), emit: pat_alpaca_output

    script:
    """
    alpaca run \
        --input_tumour_directory "$pubDir/${patient}" \
        --output_directory "$pubDir/${patient}"
    """

    stub:
    """
    touch $pubDir/${patient}/ALPACA_output_${patient}.csv
    touch $pubDir/${patient}/cn_change_to_ancestor.csv
    """
}

process Get_stats {
    publishDir "$pubDir/histograms", mode: "copy"

    input:
    val clones_eval
    val pubDir
    tuple val(patient), 
        path("*"),
        path("${patient}-refphase-segmentation.tsv"), 
        path("${patient}-refphase-phased-snps.tsv.gz"), 
        path("${patient}-refphase-sample-data-updated.tsv"),
        path("tree_paths.json"), 
        path("cp_table.csv")

    output:
    tuple val(patient), path("*"), emit: summary_stats

    script:
    """
    mkdir -p $pubDir/histograms
    Rscript $clones_eval ${patient} $pubDir "$pubDir/${patient}/ALPACA_output_${patient}.csv" ${patient}-refphase-segmentation.tsv cp_table.csv
    """
}
