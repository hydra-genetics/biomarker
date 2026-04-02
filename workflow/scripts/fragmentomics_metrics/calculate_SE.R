suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) > 0) {
  input_file <- input_args[1]
  output_file <- input_args[2]
}

# function to read in file and calculate Shannon entropy
calculate_SE <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt.gz", "")
  
  read_data <- read_tsv(input_file,
                         col_names = c("gene", "exon", "fragsize", "count"),
                         col_types = cols())
  
  output <- read_data %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    group_by(id) %>% 
    summarise(se = entropy(count)) %>% 
    ungroup() %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, se)
    
  return(output)
}

# normalize sample
sample_SE <- calculate_SE(input_file)

# set output file and write to file
write_tsv(sample_SE, output_file)
