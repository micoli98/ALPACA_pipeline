process COUNT {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: "copy"
    conda params.svclone_env
    
    input: 
    tuple val(sample),
        val(patient),
        path(bam),
        path(bai),
        path(svin),
        path(config),
        path(readinfo)

    output:
    tuple val(sample),
        path("${sample}_svinfo.txt")

    script:
    """
    svclone count \
        -i "${svin}" \
        -b "${bam}" \
        -s "${sample}" \
        -o "." \
        -cfg "${config}"
    """
}