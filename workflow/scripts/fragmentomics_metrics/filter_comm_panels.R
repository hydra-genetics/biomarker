suppressPackageStartupMessages(library(tidyverse))

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Load reads (8-column format as per specification:
# chr, start, stop, ens_id, refseq, gene, exon, strand)
reads <- read_tsv(input_file, show_col_types = FALSE, quote = "",
         col_names = c("chr", "start", "stop", "ens_id", "refseq", "gene", "exon", "strand"),
         col_types = cols(
             chr = col_character(),
             start = col_double(),
             stop = col_double(),
             ens_id = col_character(),
             refseq = col_character(),
             gene = col_character(),
             exon = col_character(),
             strand = col_character()
         ))

reads <- reads %>% filter(!is.na(start) & !is.na(stop) & start < stop)

if (ncol(reads) < 8) {
    stop(paste("Input file", input_file, "has fewer than 8 columns. Data may be malformed."))
}

reads <- reads %>% select(chr, start, stop, ens_id, refseq, gene, exon, strand)

# Output to target path
output_dir <- dirname(output_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(reads, output_file, col_names = FALSE, quote = "none")
