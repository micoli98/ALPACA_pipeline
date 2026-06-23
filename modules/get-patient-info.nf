process GET_PATIENT_INFO {
    tag "${pat}"

    input:
    tuple val(samples), val(pat)
    tuple path(segs),
        path(pp)

    output:
    tuple val(pat), path("${pat}_segments_converted.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    # Filtering the segmentation and purity ploidy files to only include samples present in the converted cf file
    segments <- read.table("${segs}", sep = '\\t', header = TRUE) |>
        filter(patient == "${pat}")
    write.table(segments, file = "${pat}_segments_converted.tsv", sep = '\\t', row.names = FALSE)
    """
}