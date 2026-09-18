process CLUSTERING {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: "copy"
    conda params.svclone_env
    
    input: 
    tuple val(sample),
        val(patient),
        path(config),
        path(readinfo),
        path(svs),
        path(snvs),
        path(pp)

    output:
    tuple val(sample), path("*")

    script:
    """
    svclone cluster \
        -s "${sample}" \
        -i "${svs}" \
        --snvs "${snvs}" \
        -p "${pp}" \
        --params "${readinfo}" \
        -o "." \
        -cfg "${config}" \
        --XX
    """
}