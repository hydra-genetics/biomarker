suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) < 3) {
  stop("Usage: Rscript calculate_depth.R <exon_sizes_file> <input_file> <output_file>")
}

# Parameterized resources
exon_sizes_file <- input_args[1]
input_file <- input_args[2]
output_file <- input_args[3]

if (!file.exists(exon_sizes_file)) {
    stop(str_c("Exon sizes file not found: ", exon_sizes_file))
}

# need to get list of all unique gene_exon values with their exon sizes
all_ids <- read_tsv(exon_sizes_file,
                    col_types = cols()) %>%
  mutate(id = str_c(gene, exon, sep = "_")) %>%
  dplyr::select(id, size)

# function to take depth data and normalize it by exon length and reads per sample
normalize_depth_data <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt", "")
  
  depth_data <- read_tsv(input_file,
                        col_names = c("gene", "exon", "count"),
                        col_types = cols())
  
  num_reads <- sum(depth_data$count)
  scale_factor <- if (num_reads == 0) NA_real_ else num_reads / 1e6
  
  output <- depth_data %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    left_join(all_ids, by = "id") %>% 
    mutate(norm_depth = if_else(is.na(scale_factor), NA_real_, (count / size) / scale_factor)) %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, norm_depth)
  
  return(output)
}

# normalize sample
sample_normalized <- normalize_depth_data(input_file)

# set output file and write to file
write_tsv(sample_normalized, output_file)
