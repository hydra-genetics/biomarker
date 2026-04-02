suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) < 2) {
  stop("Usage: Rscript filter_comm_panels.R <input_file> <output_file>")
}

input_file <- input_args[1]
output_file <- input_args[2]

# Load reads
reads <- read_tsv(input_file, show_col_types = F,
         col_names = c("chr", "start", "stop", "ens_id", "refseq", "gene", "exon", "strand"))

# Output to target path
output_dir <- dirname(output_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(reads, output_file, col_names = FALSE)
