suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) > 0) {
  input_file <- input_args[1]
  output_file <- input_args[2]
}

# function to read in file and calculate Shannon entropy
calculate_small_frag_prop <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt.gz", "")
  
  read_data <- read_tsv(input_file,
                         col_names = c("gene", "exon", "fragsize", "count"),
                         col_types = cols())
  
  output <- read_data %>% 
    mutate(bin = cut(fragsize, breaks = c(0, 150, 500), labels = c("small", "large"))) %>% 
    group_by(gene, exon, bin) %>% 
    summarise(bin_count = sum(count)) %>% 
    ungroup() %>% 
    filter(!is.na(bin)) %>%
    group_by(gene, exon) %>% 
    mutate(bin_prop = bin_count / sum(bin_count)) %>% 
    ungroup() %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    filter(bin == "small") %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, bin_prop)
    
  return(output)
}

# normalize sample
sample_output <- calculate_small_frag_prop(input_file)

# set output file and write to file
write_tsv(sample_output, output_file)
