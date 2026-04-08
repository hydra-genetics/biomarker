suppressPackageStartupMessages(library(tidyverse))

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Load reads
reads <- read_tsv(input_file, show_col_types = F, delim = "\t",
         col_names = c("chr", "start", "stop", "ens_id", "refseq", "gene", "exon", "strand"))

if (ncol(reads) < 8) {
    stop(paste("Input file", input_file, "has fewer than 8 columns. Data may be malformed."))
}

reads <- reads %>% select(chr, start, stop, ens_id, refseq, gene, exon, strand)

# Output to target path
output_dir <- dirname(output_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(reads, output_file, col_names = FALSE)
