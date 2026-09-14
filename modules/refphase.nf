process REFPHASE {
    tag "${pat}"
    publishDir "${params.pubDir}/${pat}", mode: "copy"
    conda params.conda_env

    input: 
    tuple val(pat), 
        path(segs), 
        path(snps), 
        path(pp)

    output:
    tuple val(pat), 
        path("${pat}-refphase-segmentation.tsv"), 
        path("${pat}-refphase-phased-snps.tsv.gz"), 
        path("${pat}-refphase-sample-data-updated.tsv")

    script: 
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(refphase)

    ## Get sample list
    sample_list <- read.table("${pp}", sep="\\t", header=T) |>
        pull(sample)

    ## Run refphase
    rp_data <- refphase_load(data_format = "tsv",
                            samples = sample_list,
                            tsv_prefix = "${pat}-")
    rp_res <- refphase(rp_data)

    ## Write output
    write_segs(rp_res\$phased_segs, file = "${pat}-refphase-segmentation.tsv")
    write_snps(rp_res\$phased_snps, file = "${pat}-refphase-phased-snps.tsv.gz")
    write.table(rp_res\$sample_data, file = "${pat}-refphase-sample-data-updated.tsv", sep = "\\t", row.names = FALSE)

    """
}