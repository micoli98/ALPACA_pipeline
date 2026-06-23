################################ Reconstruction check ##################################
library(tidyverse)
library(rlang)

source("/mnt/storageBig8/work/micoli/tumor_evolution/ALPACA_run/src/clones_comparison_functions.R")

arg <- commandArgs(trailingOnly = TRUE)

pat <- arg[1]
output_path <- arg[2]
segsA <- read.table(arg[3], sep=",", header=T)
segsH <- read.table(arg[4], sep="\t", header=T)
cf <- read.table(arg[5], sep=",", header=T)


# get refphase segments
segsH <- segsH |>
  rename(sample = sample_id)|>
  mutate(copyNumber = cn_a+cn_b,
         sample = paste0(sample, "_og"),
         segment = paste(chrom, start, end, sep="_"))

# Separation of the copy number values in 3 data frames
segsAR <- segsA |>
  rename(cn_a = pred_CN_A,
         cn_b = pred_CN_B) |>
  mutate(tmp_a = if_else(cn_b > cn_a, cn_b, cn_a), # invert cn_a and cn_b if cn_b > cn_a
         tmp_b = if_else(cn_b > cn_a, cn_a, cn_b),
         cn_a = tmp_a,
         cn_b = tmp_b) |>
  select(-tmp_a, -tmp_b) |>
  mutate(copyNumber = cn_a + cn_b)

cols <- c("copyNumber", "cn_a", "cn_b")
cnv_list <- lapply(cols, function(val_col) {
  # inner loop: per sample
  res_per_clone <- segsAR |>
    select(segment, clone, !!val_col) |>
    spread(clone, !!val_col)
  res_per_clone
})
names(cnv_list) <- cols

# Reconstruction of copy number for the samples
top_levels <- names(cnv_list)

reconstruction <- lapply(top_levels, function(lev) {
  process_folder(cnv_list[[lev]], cf, segsH)
})
names(reconstruction) <- top_levels

# average and histograms
hist_info <- plot_diff_hists(reconstruction, pat, top_levels = c("copyNumber", "cn_a","cn_b"))

