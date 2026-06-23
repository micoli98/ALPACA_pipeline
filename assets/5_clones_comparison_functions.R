# Function to reconstruct copynumber
# ----- Example usage -----
# segments_tbl <- your first table (segment + clone1..clone9)
# freqs_tbl    <- your second table (clone + sample columns)
# out <- reconstruct_copy_number(segments_tbl, freqs_tbl)
reconstruct_copy_number <- function(segments_tbl, freqs_tbl,
                                    segment_col = segment, clone_col = clone) {
  seg_col_nm <- as_name(enquo(segment_col))
  clone_col_nm <- as_name(enquo(clone_col))
  
  # clone names in each table
  clones_seg  <- setdiff(names(segments_tbl), seg_col_nm)
  clones_freq <- freqs_tbl[[clone_col_nm]]
  
  # overlap set and sanity check
  common <- intersect(clones_seg, clones_freq)
  if (length(common) == 0) stop("No overlapping clone names between the two tables.")
  
  # build segment-by-clone matrix (rows = segments, cols = clones)
  M <- segments_tbl %>%
    select(all_of(seg_col_nm), all_of(common)) %>%
    mutate(across(all_of(common), ~ replace_na(.x, 0))) %>%
    { as.matrix(select(., all_of(common))) }
  rownames(M) <- segments_tbl[[seg_col_nm]]
  
  # build clone-by-sample frequency matrix (rows = clones, cols = samples)
  F <- freqs_tbl %>%
    filter(.data[[clone_col_nm]] %in% common) %>%
    mutate(across(-all_of(clone_col_nm), ~ replace_na(.x, 0))) %>%
    arrange(factor(.data[[clone_col_nm]], levels = common)) %>%
    select(-all_of(clone_col_nm)) %>%
    { as.matrix(.) }
  rownames(F) <- common
  
  # matrix product: (segments x clones) %*% (clones x samples) -> (segments x samples)
  R <- M %*% F
  
  # back to tibble
  tibble(!!seg_col_nm := rownames(R)) %>%
    bind_cols(as_tibble(R, .name_repair = "minimal"))
}

# Process one folder table (e.g., clonal_info[["copyNumber"]][["005"]])
process_folder <- function(seg_by_clone_tbl, cf, rp_segs, obs_suffix = "_og") {
  # 1) reconstructed (segment × samples)
  ex <- reconstruct_copy_number(seg_by_clone_tbl, cf)
  
  # 2) observed RP: long -> wide (segment × samples)
  obs_wide <- rp_segs %>%
    select(segment, sample, copyNumber) %>%
    pivot_wider(names_from = sample, values_from = copyNumber)
  
  # 3) ensure observed columns are suffixed with _og to pair with reconstructed
  if (!any(grepl(paste0(obs_suffix, "$"), names(obs_wide)[-1]))) {
    names(obs_wide)[-1] <- paste0(names(obs_wide)[-1], obs_suffix)
  }
  
  # 4) join and compute absolute diffs: |obs_og - reconstructed|
  out <- ex %>%
    inner_join(obs_wide, by = "segment") %>%
    mutate(
      across(
        ends_with(obs_suffix),
        ~ abs(. - get(sub(paste0(obs_suffix, "$"), "", cur_column()))),
        .names = "diff_{.col}"
      )
    ) %>%
    rename_with(~ sub(paste0(obs_suffix, "$"), "", .x), starts_with("diff_"))
  
  out
}

# Function to get histograms of the difference columns
# ---- Run it ----
# hist_info <- plot_diff_hists(result, top_levels = c("cn_a","cn_b"))
# This will save PNGs under ./diff_histograms and return a nested list with plot objects & file paths.
plot_diff_hists <- function(result,
                            patient,
                            top_levels = c("cn_a", "cn_b"),
                            bins = 30) {

  out <- result[top_levels] |>
    imap(function(tbl, top_name) {   # <-- accept BOTH tbl and name
      # grab all diff_* columns
      diff_cols <- grep("^diff_", names(tbl), value = TRUE)
      if (length(diff_cols) == 0) {
        message(sprintf("[%s] No diff_* columns found, skipping.", top_name))
        return(NULL)
      }
      
      # long format: sample, diff
      diffs_long <- tbl |>
        select(all_of(diff_cols)) |>
        pivot_longer(everything(),
                     names_to = "sample",
                     values_to = "diff") |>
        mutate(sample = sub("^diff_", "", sample))
      
      # per-sample means for vline
      means <- diffs_long |>
        group_by(sample) |>
        summarise(mu = mean(diff, na.rm = TRUE), .groups = "drop")
      
      p <- ggplot(diffs_long, aes(x = diff)) +
        geom_histogram(bins = bins, na.rm = TRUE) +
        geom_vline(data = means, aes(xintercept = mu), linetype = "dashed") +
        geom_text(data = means,
                  aes(x = mu, y = 0, label = sprintf("mean=%.2f", mu)),
                  angle = 90, vjust = -0.5, hjust = 0, size = 3, inherit.aes = FALSE) +
        facet_wrap(~ sample, scales = "free_y") +
        labs(
          title = sprintf("%s: absolute difference (observed vs reconstructed)", top_name),
          x = "Absolute difference",
          y = "Count"
        ) +
        theme_minimal(base_size = 12)
      
      file_png <- file.path(paste0(patient, "_", top_name, "_diff_hist.png"))
      ggsave(filename = file_png, plot = p, width = 12, height = 8, dpi = 150, bg = "white")
      
      list(plot = p, file = file_png,
           n_segments = nrow(tbl),
           n_samples  = length(diff_cols))
    })
  
  invisible(out)
}

### Function to add segment length in the dataframes
# Add a length column to one data frame that has a `segment` column "chrom_start_end"
add_seg_length <- function(df, seg_col = segment, new_col = "seg_len") {
  seg_col <- rlang::ensym(seg_col)
  
  df %>%
    # split "chrom_start_end" into temporary columns; keep original `segment`
    separate({{ seg_col }},
             into = c(".chrom", ".start", ".end"),
             sep = "_", remove = FALSE, convert = TRUE) %>%
    mutate(!!new_col := .end - .start + 1) %>%    # interval length
    select(-.chrom, -.start, -.end)               # drop temps
}
