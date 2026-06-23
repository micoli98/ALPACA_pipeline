# Function to get segment confidence intervals proportional to the length
get_proportional_ci <- function(df,
                                avg_width = 0.25,
                                min_width = NA_real_,
                                max_width = NA_real_) {
  stopifnot(is.data.frame(df))
  req_cols <- c("cn_a", "cn_b")
  missing_req <- setdiff(req_cols, names(df))
  if (length(missing_req)) stop("Missing columns: ", paste(missing_req, collapse = ", "))
  
  # Obtain segment length
  if ("width" %in% names(df)) {
    seg_len <- df$width
  } else if (all(c("start", "end") %in% names(df))) {
    seg_len <- df$end - df$start + 1
  } else {
    stop("Provide either a 'length' column or both 'start' and 'end'.")
  }
  
  # Clean/guard
  if (anyNA(seg_len)) stop("Segment lengths contain NA.")
  if (any(seg_len <= 0)) stop("All segment lengths must be > 0.")
  if (!is.numeric(df$cn_a) || !is.numeric(df$cn_b)) {
    stop("'cn_a' and 'cn_b' must be numeric.")
  }
  
  # Scale factor so mean(width) = avg_width
  scale <- avg_width / mean(seg_len)
  
  width <- seg_len * scale
  
  # Optional caps
  if (!is.na(min_width)) width <- pmax(width, min_width)
  if (!is.na(max_width)) width <- pmin(width, max_width)
  
  # Build ranges centered at the CN values
  halfw <- width / 2
  major_lo <- ifelse(df$cn_a - halfw <0 , 0, df$cn_a - halfw)
  major_hi <- df$cn_a + halfw
  minor_lo <-  ifelse(df$cn_b - halfw <0 , 0, df$cn_b - halfw)
  minor_hi <- df$cn_b + halfw
  
  out <- df |>
    mutate(segment = paste(chrom, start, end, sep ="_")) |>
    select(sample_id, segment, cn_a, cn_b, width)
  out$width    <- width
  out$lower_CI_A <- major_lo
  out$upper_CI_A <- major_hi
  out$lower_CI_B <- minor_lo
  out$upper_CI_B <- minor_hi
  out$width_A <- major_hi - major_lo
  out$width_B <- minor_hi - minor_lo
  
  out
}
