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
