################################ Preparation of Refphase input ################################ 

library(GenomicRanges)
library(tidyverse)

# home <- "~/mnt/work/micoli"
# data_path <- file.path(home, "SCNA_Purple/results/250701")
# output_path <- file.path(home, "tumor_evolution/ALPACA-model/CI_intervals")

#' Input required:
#' - patient
#' - output path
#' - segmentation
#' - SNPs path
#' - purity-ploidy

arg <- commandArgs(trailingOnly = TRUE)

pat <- arg[1]
output_path <- arg[2]
seg <- read.table(arg[3], sep="\t", header=T)
snp_path <- arg[4]
pp <- read.table(arg[5], sep="\t", header=T)
cl_info <- read.table(arg[6], sep="\t", header=T)

### Pt1: Harmonization of CNVs between samples ###
# From Purple output, derive a similar one to GATK, so all samples get the same number of segments

# Function to map and average copy numbers
get_harmonized_segments <- function(sample_name) {
  gr_sample <- gr_all[gr_all$sample == sample_name]
  
  overlaps <- findOverlaps(unified_ranges, gr_sample)
  
  df <- tibble(
    chrom = as.character(seqnames(unified_ranges[queryHits(overlaps)])),
    start = start(unified_ranges[queryHits(overlaps)]),
    end = end(unified_ranges[queryHits(overlaps)]),
    orig_major = mcols(gr_sample[subjectHits(overlaps)])$major,
    orig_minor = mcols(gr_sample[subjectHits(overlaps)])$minor,
  ) %>%
    group_by(chrom, start, end) %>%
    summarize(
      majorAlleleCopyNumber = mean(orig_major),
      minorAlleleCopyNumber = mean(orig_minor),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_name)
  
  return(df)
}

# Patient harmonization
valid_samples <- read.table(cl_info[cl_info$patient == pat, "file_cf"], sep = '\t', header = TRUE) |>
  pull(sample.id) |>
  unique()

sp <- seg[seg$patient == pat & seg$sample %in% valid_samples, ]

# Create GRanges object
gr_all <- GRanges(
  seqnames = sp$chromosome,
  ranges = IRanges(start = sp$start, end = sp$end),
  major = sp$majorAlleleCopyNumber,
  minor = sp$minorAlleleCopyNumber,
  sample = sp$sample
)

# 1. Get unified segmentation
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
harmonized <- map_dfr(valid_samples, get_harmonized_segments) |>
  select(sample, everything()) |>
  setNames(c("sample_id", "chrom", "start", "end", "cn_major", "cn_minor")) |>
  mutate(chrom = gsub("chr", "", chrom))
write.table(harmonized, paste0(pat, "-segments.tsv"), sep="\t", col.names = T, row.names = F)

### Pt2: Refine SNPs and purity/ploidy format ###
# SNPs input
snps_all <- data.frame()
for(s in valid_samples) {
  snps <- read.table(file.path(snp_path, pat, paste0(s, ".amber.baf.tsv.gz")), sep="\t", header=T)
  snps_mod <- snps |>
    mutate(
      chrom = gsub("chr", "", chromosome),
      pos = position,
      baf = tumorModifiedBAF,
      germline_baf = normalModifiedBAF,
      logR = log2((tumorDepth + 1) / (normalDepth + 1)),  # with pseudocount
      sample_id = s
    ) %>%
    select(chrom, pos, baf, germline_baf, logR, sample_id)
  snps_all <- bind_rows(snps_all, snps_mod)
}
write.table(snps_all,  paste0(pat, "-snps.tsv"), sep="\t", col.names = T, row.names = F)

# Purity ploidy input
pp <- pp |>
  filter(sample %in% valid_samples) |>
  select(sample, purity, ploidy) 
write.table(pp, paste0(pat, "-purity_ploidy.tsv"), sep="\t", col.names = T, row.names = F)



