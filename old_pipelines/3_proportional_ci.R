################################ Proportional width confidence intervals ################################
library(tidyverse)
source("/mnt/storageBig8/work/micoli/tumor_evolution/ALPACA_run/src/ci_functions.R")

# home <- "/mnt/storageBig8/work/micoli/"
# output <- file.path(home, "tumor_evolution/ALPACA-model/CI_intervals")

#' Input required:
#' - patient
#' - pubDir

arg <- commandArgs(trailingOnly = TRUE)

pat <- arg[1]
output_path <- arg[2]

# Load data and correct negative copy number values
segs <- read.table(file.path(output_path, pat, paste0(pat, "-refphase-segmentation.tsv")), sep="\t", header=T) |>
  mutate(cn_a = ifelse(cn_a < 0, 0, cn_a))

### Calcualte confidence intervals
# If you want to avoid extremes, add caps (e.g., never narrower than 0.05)
cis <- get_proportional_ci(segs, avg_width = 0.25, min_width = 0.05)|>
  rename(sample = sample_id) |>
  mutate(ci_value = 0.5,
         tumour_id = pat) |>
  select(segment, sample, lower_CI_A, upper_CI_A, lower_CI_B, upper_CI_B, tumour_id, ci_value)
write.table(cis, "ci_table.csv", sep=",", col.names = T, row.names = F)

### Fix the segmentation file
s <- segs |>
  rename(cpnA = cn_a,
         cpnB = cn_b,
         sample = sample_id) |>
  mutate(segment = paste(chrom, start, end, sep="_"),
         tumour_id = pat) |>
  select(tumour_id, sample, segment, cpnA, cpnB)
write.table(s, "ALPACA_input_table.csv", sep=",", col.names = T, row.names = F)
