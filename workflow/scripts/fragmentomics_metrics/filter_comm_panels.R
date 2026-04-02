suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) < 1) {
  stop("Usage: Rscript filter_comm_panels.R <input_file>")
}

input_file <- input_args[1]
output_file <- str_replace(basename(input_file), ".gz", "")

# Load reads
reads <- read_tsv(input_file, show_col_types = F,
         col_names = c("chr", "start", "stop", "ens_id", "refseq", "gene", "exon", "strand"))

# Output directly to data directory
# In Hydra, this might need to be adjusted depending on the rule's output path
output_dir <- "data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(reads, file.path(output_dir, output_file), col_names = FALSE)
