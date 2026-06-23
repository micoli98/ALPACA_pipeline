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
        path("cn_change_to_ancestor.csv"),
        emit: results
    path "*_report.csv", optional: true, emit: reports
    path "run_gap_summary.csv", optional: true, emit: gap_summary

    script:
    """
    alpaca run \\
        --input_tumour_directory . \\
        --output_directory ./output_${pat} \\
        --solver gurobi \\
        --genome_build hg19 \\
        --extra_columns complexity CI_score D_score

    cp ./output_${pat}/ALPACA_output_${pat}.csv .
    cp ./output_${pat}/cn_change_to_ancestor.csv .
    for f in ci_modified_report.csv monoclonal_samples_report.csv run_gap_summary.csv infeasibility_report.csv; do
        [ -f ./output_${pat}/\$f ] && cp ./output_${pat}/\$f . || true
    done
    """
}