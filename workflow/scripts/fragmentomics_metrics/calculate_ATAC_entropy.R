suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

sample_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

calculate_atac_entropy <- function(input_file) {
  
  sample_name <- str_replace(basename(input_file), "_ATAC_frag_count.txt.gz", "")
  
  data <- read_table(input_file,
                     col_names = c("count", "cancer", "fragsize"), 
                     col_types = cols())
  
  if (nrow(data) == 0) {
      stop(paste("FATAL: Input data is empty for sample:", sample_name))
  }

  output <- data %>% 
    group_by(cancer) %>% 
    summarise(ATAC_entropy = entropy(count)) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, cancer, ATAC_entropy)
  
  if (nrow(output) == 0) {
      stop(paste("FATAL: Entropy calculation results are empty for sample:", sample_name))
  }

  return(output)
}

# calculate entropy by TF
sample_ATAC_entropy <- calculate_atac_entropy(sample_file)

# write to file
write_tsv(sample_ATAC_entropy, output_file)
