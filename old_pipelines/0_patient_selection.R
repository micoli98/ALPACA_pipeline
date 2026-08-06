## Determine which patients to analyze with ALPACA
library(tidyverse)

setwd("~/mnt/work/micoli/tumor_evolution")

# Load data
batch_info <- read.table("clones_analysis/batch_info.tsv", sep="\t", header=T)
av_samples <- read.table("../SCNA_Purple/results/250701/purity_ploidy_estimates.tsv", sep="\t", header=T)

# Select samples with relapses and cellular frequencies avalable
r_patients <- av_samples |>
  mutate(phase = sub( '^(\\w+\\d+)_([piro])(.*)', '\\2', sample)) |>
  group_by(patient) |>
  filter(any(phase == "r")) |>
  ungroup() |> 
  left_join(batch_info[, c("patient", "batch")]) |>
  filter(!is.na(batch))

# Patient list
p <- r_patients |> select(sample, patient) |> unique()
write.table(p, "ALPACA_run/results/patients_251112.txt", sep="\t", col.names = T, row.names = F, quote = F)
