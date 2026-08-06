################################ Refphase ################################ 
library(tidyverse)
library(refphase)

# home <- "/mnt/storageBig8/work/micoli/"
# output <- file.path(home, "tumor_evolution/ALPACA-model/CI_intervals")

#' Input required:
#' - patient
#' - path to RP files

arg <- commandArgs(trailingOnly = TRUE)

pat <- arg[1]
output_path <- arg[2]
pp <- read.table(arg[3], sep="\t", header=T)

## Get sample list
sample_list <- pp |>
  pull(sample)

## Run refphase
rp_data <- refphase_load(data_format = "tsv",
                         samples = sample_list,
                         tsv_prefix = file.path(output_path, pat, paste0(pat, "-")))
rp_res <- refphase(rp_data)

## Write output
write_segs(rp_res$phased_segs, file = paste0(pat, "-refphase-segmentation.tsv"))
write_snps(rp_res$phased_snps, file = paste0(pat, "-refphase-phased-snps.tsv.gz"))
write.table(rp_res$sample_data, file = paste0(pat, "-refphase-sample-data-updated.tsv"), sep = "\t", row.names = FALSE)
