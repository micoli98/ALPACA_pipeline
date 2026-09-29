process HARMONIZE {
    tag "${pat}"
    publishDir "${params.outdir}/${pat}/alpaca", mode: "copy" 
    conda params.alpaca_env
    
    input:
    tuple val(pat),
        path(segs),
        path(pp),
        val(model),
        path(cf_file),
        path(tree_file),
        path(snp_files),
        path(rd_files)

    output:
    tuple val(pat), 
        path("${pat}-segments.tsv"),
        path("${pat}-snps.tsv"),
        path("${pat}-purity_ploidy.tsv")

    script: 
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(tidyverse)
    })
    source("${params.harmonization_function}")

    ### Pt1: Harmonization of CNVs between samples ###
    # From Purple output, derive a similar one to GATK, so all samples get the same number of segments

    #valid_samples <- read.table("${cf_file}", sep = '\\t', header = TRUE) |>
    #    pull(sample.id) |>
    #    unique()
    
    sp <- read.table("${segs}", sep = '\\t', header = TRUE) #|>
    #    filter(sample %in% valid_samples)

    # Create GRanges object
    gr_all <- GRanges(
        seqnames = sp\$chromosome,
        ranges = IRanges(start = sp\$start, end = sp\$end),
        major = sp\$majorAlleleCopyNumber,
        minor = sp\$minorAlleleCopyNumber,
        sample = sp\$sample
        )
    
    # Get unified segmentation
    # Get all chromosomes
    all_chroms <- unique(seqnames(gr_all))

    # Compute unified segments per chromosome
    unified_list <- lapply(all_chroms, function(chr) {
        gr_chr <- gr_all[seqnames(gr_all) == chr]
        breaks <- unique(sort(c(start(gr_chr), end(gr_chr))))
        # Avoid intervals of length 0
        valid <- which(diff(breaks) > 1)
        GRanges(
            seqnames = chr,
            ranges = IRanges(
            start = breaks[valid],
            end = breaks[valid + 1]  # make end inclusive
            )
        )
    })

    # Combine into one GRanges object
    unified_ranges <- do.call(c, unified_list)

    # 2. Assign unified segments to each sample
    # Apply to all samples
    valid_samples <- unique(sp\$sample)
    harmonized <- map_dfr(valid_samples, get_harmonized_segments) |>
        select(sample, everything()) |>
        setNames(c("sample_id", "chrom", "start", "end", "cn_major", "cn_minor")) |>
        mutate(chrom = gsub("chr", "", chrom))
    write.table(harmonized, "${pat}-segments.tsv", sep="\\t", col.names = T, row.names = F)
    
    ### Pt2: Refine SNPs and purity/ploidy format ###
    # SNPs input
    # Get list of SNP files
    # Get list of SNP files from the working directory
    snp_files <- unlist(strsplit("${snp_files}", split = " "))
    print(snp_files)

    # Extract sample names from filenames
    sample_names <- basename(snp_files) |>
        str_remove(".amber.baf.tsv(.gz)*\$")

    # Keep only samples present in harmonized segments
    keep <- sample_names %in% harmonized\$sample_id
    sample_names <- sample_names[keep]
    snp_files    <- snp_files[keep] 

    # COBALT read-depth files, named by sample
    rd_files <- unlist(strsplit("${rd_files}", split = " "))
    names(rd_files) <- str_remove(basename(rd_files), fixed(".cobalt.ratio.tsv.gz"))

    # GC-corrected logR of the COBALT 1 kb window containing each SNP.
    # Windows with ratio -1 (masked) or 0 (no coverage) are skipped, so their SNPs get NA.
    cobalt_logr <- function(snps, rd_file) {
        rd <- read_tsv(rd_file, show_col_types = F,
                       col_select = c(chromosome, position, tumorGCRatio, referenceGCDiploidRatio)) |>
            filter(tumorGCRatio > 0, referenceGCDiploidRatio > 0)
        windows <- GRanges(gsub("chr", "", rd\$chromosome), IRanges(rd\$position, width = 1000))
        hits <- findOverlaps(GRanges(snps\$chrom, IRanges(snps\$pos, width = 1)), windows, select = "first")
        log2(rd\$tumorGCRatio / rd\$referenceGCDiploidRatio)[hits]
    }

    # Read and process all files at once
    snps_all <- map2_dfr(snp_files, sample_names, function(file, sample) {
        snps <- read_tsv(file, show_col_types = F) |>
            transmute(
                chrom = gsub("chr", "", chromosome),
                pos = position,
                baf = tumorBAF,
                germline_baf = normalModifiedBAF,
                sample_id = sample
            )
        snps\$logR <- cobalt_logr(snps, rd_files[[sample]])
        select(snps, chrom, pos, baf, germline_baf, logR, sample_id)
    })

    # RefPhase's validate_snps() rejects NA logR: drop SNPs in masked COBALT windows
    snps_all |>
        group_by(sample_id) |>
        summarise(pct_dropped = round(100 * mean(is.na(logR)), 2)) |>
        print()
    snps_all <- filter(snps_all, !is.na(logR))

    # Save output
    write.table(snps_all, "${pat}-snps.tsv", sep="\\t", quote=FALSE, row.names=FALSE)

    # Purity ploidy input
    pp <- read.table("${pp}", sep="\\t", header=T) |>
        filter(sample %in% valid_samples) |>
        select(sample, purity, ploidy) 
    write.table(pp, "${pat}-purity_ploidy.tsv", sep="\\t", col.names = T, row.names = F)

    """
}



