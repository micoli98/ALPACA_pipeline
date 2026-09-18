process FILTER_SV {
    tag "${sample}"
    conda params.svclone_env
    publishDir "${params.outdir}/${sample}", mode: "copy"

    input: 
    tuple val(sample),
        val(patient),
        path(vcf_sv)

    output:
    tuple val(sample),
        val(patient),
        path("${sample}.purple.breakpoints.vcf")

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(StructuralVariantAnnotation)
    })

    sv_file <- readVcf("${vcf_sv}", genome="hg38")

    # Extract breakpoint representation as GRanges
    bp <- breakpointRanges(sv_file)

    # Retrieve the original VCF rows that produced valid breakpoints
    bp_ids <- unique(bp\$sourceId)
    vcf_bp <- sv_file[rownames(sv_file) %in% bp_ids, ]

    writeVcf(vcf_bp, "${sample}.purple.breakpoints.vcf")
    """
}