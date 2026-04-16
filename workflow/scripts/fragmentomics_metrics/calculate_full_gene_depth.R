suppressPackageStartupMessages(library(tidyverse))

# Parameterized resources
exon_sizes_file <- snakemake@params[["exon_sizes"]]
cds_bed_file <- snakemake@params[["cds_bed"]]
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

if (is.null(exon_sizes_file) || !file.exists(exon_sizes_file)) {
    stop(paste("FATAL: Exon sizes file not found or NULL:", exon_sizes_file))
}

# Load Targeted Exons from BED
targeted_exons <- read_tsv(cds_bed_file, 
                          col_names = c("chr", "start", "stop", "transcript", "refseq", "gene", "exon", "strand"),
                          col_types = cols()) %>%
  mutate(id = str_c(gene, exon, sep = "_")) %>%
  pull(id) %>%
  unique()

# need to get list of all unique gene_exon values with their exon sizes
all_ids <- read_tsv(exon_sizes_file,
                    col_types = cols()) %>% 
  filter(str_detect(gene, "_", negate = TRUE)) %>% 
  mutate(id = str_c(gene, exon, sep = "_")) %>% 
  filter(id %in% targeted_exons) %>% # Filter to targeted exons only
  dplyr::select(id, size)

gene_sizes <- all_ids %>% 
  separate(id, into = c("gene", "exon"), sep = "_") %>% 
  group_by(gene) %>% 
  summarise(gene_size = sum(size))

# function to take depth data and normalize it by exon length and reads per sample
normalize_full_depth_data <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt", "")
  
  depth_data <- read_tsv(input_file,
                        col_names = c("gene", "exon", "count"),
                        col_types = cols())
  
  if (nrow(depth_data) == 0) {
      stop(paste("FATAL: Input data is empty for sample:", sample_name))
  }

  num_reads <- sum(depth_data$count)
  scale_factor <- if (num_reads == 0) NA_real_ else num_reads / 1e6
  
  output <- depth_data %>% 
    group_by(gene) %>% 
    summarise(gene_count = sum(count)) %>% 
    left_join(gene_sizes, by = "gene") %>% 
    filter(!is.na(gene_size)) %>% 
    mutate(norm_depth = (gene_count / gene_size) / scale_factor) %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, gene, norm_depth)
  
  if (nrow(output) == 0) {
      stop(paste("FATAL: Normalization results are empty for sample:", sample_name))
  }

  return(output)
}

# normalize sample
sample_normalized <- normalize_full_depth_data(input_file)

# set output file and write to file
write_tsv(sample_normalized, output_file)
