# Function to trace from leaf to root
trace_path <- function(label) {
  path <- character()
  current <- label
  while (!is.na(current) && label_to_parent[[current]] != "clone-1") {
    path <- c(current, path)
    current <- label_to_parent[[current]]
  }
  path <- c("clone1", path)
  return(path)
}

# Transform to a string
format_list_as_string <- function(lineages) {
  items <- lapply(lineages, function(path) {
    paste0('["', paste(path, collapse = '", "'), '"]')
  })
  result <- paste0("[", paste(items, collapse = ", "), "]")
  return(result)
}
