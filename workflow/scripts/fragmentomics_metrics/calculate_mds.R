suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

left_file <- snakemake@input[["left"]]
right_file <- snakemake@input[["right"]]
output_file <- snakemake@output[[1]]

calculate_mds <- function(left_data, right_data) {
  
  sample_name <- str_replace(basename(left_data), "_left_4mers.txt.gz", "")
  
  left_motifs <- read_table(left_data,
                          col_names = c("count", "gene", "exon", "motif"),
                          col_types = cols())
  right_motifs <- read_table(right_data,
                           col_names = c("count", "gene", "exon", "motif"),
                           col_types = cols())
  
  if (nrow(left_motifs) == 0 && nrow(right_motifs) == 0) {
      stop(paste("FATAL: Both left and right motif inputs are empty for sample:", sample_name))
  }

  combined_data <- bind_rows(left_motifs, right_motifs) %>% 
    group_by(gene, exon, motif) %>% 
    summarise(
      count = sum(count),
      .groups = "drop"
    ) %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    dplyr::select(id, motif, count)
  
  output <- combined_data %>% 
    group_by(id) %>% 
    summarise(
      mds = entropy(count)
    ) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, id, mds)
  
  if (nrow(output) == 0) {
      stop(paste("FATAL: MDS calculation produced zero results for sample:", sample_name))
  }

  return(output)
}

sample_mds <- calculate_mds(left_data = left_file, right_data = right_file)

# write to file
write_tsv(sample_mds, output_file)
