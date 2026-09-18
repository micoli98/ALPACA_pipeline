process ANNOTATION {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", pattern: "*.txt", mode: "copy"
    conda params.svclone_env
    errorStrategy 'ignore'
    
    input: 
    tuple val(sample),
        val(patient),
        path(bam),
        path(bai),
        path(vcf_sv),
        path(config)

    output:
    tuple val(sample),
        path("${sample}_svin.txt"),
        path("read_params.txt")

    script:
    """
    # Compress the VCF
    bgzip -c ${vcf_sv} > ${sample}.purple.breakpoints.vcf.gz
    tabix -p vcf ${sample}.purple.breakpoints.vcf.gz

    svclone annotate \
        -i "${sample}.purple.breakpoints.vcf.gz" \
        -b "${bam}" \
        -s "${sample}" \
        -o "." \
        --sv_format vcf \
        -cfg "${config}"
    """
}