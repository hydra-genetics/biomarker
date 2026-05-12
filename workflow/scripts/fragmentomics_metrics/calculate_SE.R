suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# function to read in file and calculate Shannon entropy
calculate_SE <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt.gz", "")
  
  read_data <- read_tsv(input_file,
                         col_names = c("gene", "exon", "fragsize", "count"),
                         col_types = cols())
  
  if (nrow(read_data) == 0) {
      message(paste("WARNING: Input data is empty for sample:", sample_name))
      return(tibble(sample = character(), id = character(), se = numeric()))
  }

  output <- read_data %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    group_by(id) %>% 
    summarise(se = entropy(count), .groups = "drop") %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, se)
  
  if (nrow(output) == 0) {
      message(paste("WARNING: Shannon entropy calculation produced zero results for sample:", sample_name))
      return(tibble(sample = character(), id = character(), se = numeric()))
  }
    
  return(output)
}

# normalize sample
sample_SE <- calculate_SE(input_file)

# set output file and write to file
write_tsv(sample_SE, output_file)
