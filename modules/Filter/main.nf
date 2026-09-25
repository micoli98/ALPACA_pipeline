process FILTER {
    tag "${sample}"
    publishDir "${params.outdir}/${patient}/svclone/${sample}", mode: "copy"
    conda params.svclone_env
    // svclone rewrites read_params.txt in place; copy so the upstream work-dir file (and cache hash) stays untouched
    stageInMode 'copy'
    
    input: 
    tuple val(sample),
        val(patient),
        path(config),
        path(readinfo),
        path(svinfo),
        path(snv),
        path(cnv),
        path(pp)

    output:
    tuple val(sample),
        path("${sample}_filtered_svs.tsv"),
        path("${sample}_filtered_snvs.tsv"),
        path("purity_ploidy.txt")

    script:
    """
    svclone filter \
        -s "${sample}" \
        -i "${svinfo}" \
        -o "." \
        -cfg "${config}" \
        --params "${readinfo}" \
        -c "${cnv}" \
        -p "${pp}" \
        --snvs "${snv}" \
        --snv_format mutect
    """
}