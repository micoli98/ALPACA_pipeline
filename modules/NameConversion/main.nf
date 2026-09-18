process NAME_CONVERSION {
    tag "${pat}"
    publishDir "${params.outdir}/${pat}", pattern: "*cf_converted.tsv", mode: "copy"

    input:
    tuple val(pat), path(cf_file)
    tuple path(segs),
        path(pp)

    output:
    tuple val(pat), 
        path("${pat}_cf_converted.tsv"), 
        path("${pat}_segments_converted.tsv"),
        path("${pat}_purity_ploidy_converted.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    conv_table <- read.table("${params.conversion_table}", sep="\\t", header=T) |>
        select(bamName, platform, id, current)

    # Name conversion
    cf_orig <- read.table("${cf_file}", sep = '\\t', header = TRUE)
    cf_names <- cf_orig |>
        select(sample.id) |>
        unique() |>
        rename(old_name = sample.id)
    for(n in 1:nrow(cf_names)) {
        old_name <- cf_names[n, "old_name"]
        
        initial_row <- as.integer(rownames(conv_table[conv_table\$bamName == old_name & conv_table\$current, ]))
        if (length(initial_row) == 0) {
            initial_row <- as.integer(rownames(conv_table[conv_table\$bamName == old_name, ]))
        }
        if (length(initial_row) > 1) {
            initial_row <- initial_row[-1]
        }
        while(!conv_table[initial_row, "current"] ) {
            initial_row <- initial_row - 1
        }
        cf_names[n, "new_name"] <- conv_table[initial_row, "bamName"]
    }
    
    cf_new <- cf_orig |>
        left_join(cf_names, by = c("sample.id" = "old_name")) |>
        mutate(sample.id = new_name) |>
        select(-new_name)

    # Filtering the segmentation and purity ploidy files to only include samples present in the converted cf file
    segments <- read.table("${segs}", sep = '\\t', header = TRUE) |>
        filter(sample %in% cf_new\$sample.id)
    write.table(segments, file = "${pat}_segments_converted.tsv", sep = '\\t', row.names = FALSE)

    pp <- read.table("${pp}", sep = '\\t', header = TRUE) |>
        filter(sample %in% cf_new\$sample.id)
    write.table(pp, file = "${pat}_purity_ploidy_converted.tsv", sep = '\\t', row.names = FALSE)

    # Flter again cellular fractions to exclude cell lines samples
    cf_new <- cf_new |>
        filter(sample.id %in% unique(segments\$sample))
    write.table(cf_new, file = "${pat}_cf_converted.tsv", sep = '\\t', row.names = FALSE)
    """
}


//     # Read input file
//     conv_table <- read.table("${conv_table}", sep="\\t", header=T) |>
//         group_by(sample) %>%
//         mutate(bamName = bamName[current == TRUE][1]) %>%  # take the TRUE one
//         ungroup() %>%
//         filter(!is.na(prevBamName)) |>
//         select(sample, bamName, prevBamName) |>
//         filter(!prevBamName %in% bamName)

//     # Name conversion
//     cf_new <- read.table("${cf_file}", sep = '\\t', header = TRUE) |>
//         # match sample.id to old_name
//         left_join(conv_table, by = c("sample.id" = "prevBamName")) %>%
//         # if there's a new_name, use it; otherwise keep sample.id
//         mutate(sample.id = coalesce(bamName, sample.id)) |>
//         select(model.num, sample.id, cloneID, cell.freq) |>
//         distinct()

//     # Write output file
//     write.table(cf_new, file = "${pat}_cf_converted.tsv", sep = '\\t', row.names = FALSE)
// ######################################################################################################