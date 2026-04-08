suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

sample_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

calculate_tfbs_entropy <- function(input_file) {
  
  sample_name <- str_replace(basename(input_file), "_TFBS_frag_count.txt.gz", "")
  
  output <- read_table(input_file,
                             col_names = c("count", "frag_size", "tfbs"),
                             col_types = cols()) %>% 
    group_by(tfbs) %>% 
    summarise(TF_entropy = entropy(count)) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, tfbs, TF_entropy)
  
  return(output)
}

# calculate entropy by TF
sample_TFBS_entropy <- calculate_tfbs_entropy(sample_file)

# write to file
write_tsv(sample_TFBS_entropy, output_file)
