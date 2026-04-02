suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) > 0) {
  sample_file <- input_args[1]
  output_file <- input_args[2]
}

calculate_atac_entropy <- function(input_file) {
  
  sample_name <- str_replace(basename(input_file), "_ATAC_frag_count.txt.gz", "")
  
  output <- read_table(input_file,
                             col_names = c("count", "cancer", "fragsize"), 
                             col_types = cols()) %>% 
    group_by(cancer) %>% 
    summarise(ATAC_entropy = entropy(count)) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, cancer, ATAC_entropy)
  
  return(output)
}

# calculate entropy by TF
sample_ATAC_entropy <- calculate_atac_entropy(sample_file)

# write to file
write_tsv(sample_ATAC_entropy, output_file)
