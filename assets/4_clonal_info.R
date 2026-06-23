################################ Clonal information reformatting ################################
library(tidyverse)
source("/mnt/storageBig8/work/micoli/tumor_evolution/ALPACA_run/src/trees_functions.R")
#home <- "~/mnt/work/micoli/tumor_evolution/ALPACA-model/CI_intervals"

#' Input required:
#' - patient
#' - pubDir
#' - clonal_info

arg <- commandArgs(trailingOnly = TRUE)

pat <- arg[1]
output_path <- arg[2]
cl_info <- read.table(arg[3], sep="\t", header=T)

####### Processing ####### 
file_path <- cl_info[cl_info$patient == pat, "file_cf"]

### Clone proportions table ###
# Read cellular frequency file
cf <- read.table(file_path, sep = '\t', header = TRUE)

# Get the correct model and corresponding frequencies
model_num <- cl_info[cl_info$patient == pat, "model"]
cf <- cf[cf$model.num == model_num, ]

cf_table <- cf |> 
  select(-model.num) |>
  mutate(cell.freq = cell.freq/100) |>
  spread(sample.id, cell.freq, fill = 0) |>
  rename(clone = cloneID) %>%
  mutate(clone = paste0("clone", clone),
         across(-1, ~ . / sum(.))) # ensure sum of columns is 1

### Tree structure ###
fileTree_path <- cl_info[cl_info$patient == pat, "file_tree"]

tree <- readRDS(fileTree_path)
tree <- tree$matched$merged.trees[[model_num]]

tree$lab <- paste0("clone", tree$lab)
tree <- tree |> 
  select(lab, parent) |>
  mutate(parent=paste0("clone", parent))

# Create lookup: lab -> parent
label_to_parent <- setNames(tree$parent, tree$lab)

# Find leaves: nodes that are not parents
parent_labels <- unique(tree$parent)
leaf_labels <- setdiff(tree$lab, parent_labels)

# Trace all paths
lineages <- map(leaf_labels, trace_path)
output_string <- format_list_as_string(lineages)

# Check all clones in the tree are also in cf_table
cf_table <- tibble(clone = unique(c(parent_labels[-1], leaf_labels))) %>%
  left_join(cf_table, by = "clone") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

####### Printing ####### 
write_file(output_string, "tree_paths.json")
write.table(cf_table, "cp_table.csv", sep=",", col.names = T, row.names = F)










