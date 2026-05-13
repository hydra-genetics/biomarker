suppressPackageStartupMessages(library(tidyverse))

# Get Strand Info from snakemake params
strand_file <- snakemake@params[["strand_mapping"]]
if (is.null(strand_file) || !file.exists(strand_file)) {
    stop(paste("Strand mapping file not found or NULL:", strand_file))
}

gene_strands <- read_tsv(strand_file, show_col_types = F) %>% 
  dplyr::select(gene = gene_symbol, strand) %>% 
  unique()


feature_tables <- c("depth","mds", "se", "small_frags")

# Extract First Exon
for (ft in feature_tables) {
    print(str_c("Filtering E1 for", ft, sep = " "))
    
    rds_inpath <- snakemake@input[[ft]]
    outpath <- snakemake@output[[paste0(ft, "_E1")]]
    
    if (is.null(rds_inpath) || !file.exists(rds_inpath)) {
      warning(str_c("Input RDS not found: ", rds_inpath))
      next
    }
    
    full_exon_data <- read_rds(rds_inpath)
    
    first_exon_data <- full_exon_data %>% 
      pivot_longer(cols = -sample, names_to = "feature", values_to = "value") %>% 
      separate(feature, into = c("gene", "exon"), sep = "_", remove = FALSE) %>% 
      left_join(gene_strands, by = "gene") %>% 
      mutate(exon = as.integer(exon)) %>% 
      group_by(gene) %>% 
      mutate(
        first_exon = ifelse(strand == "+" & exon == min(exon), TRUE, FALSE),
        first_exon = ifelse(strand == "-" & exon == max(exon), TRUE, first_exon)
      ) %>% 
      filter(first_exon) %>% 
      ungroup() %>% 
      dplyr::select(sample, feature, value) %>% 
      pivot_wider(names_from = feature, values_from = value)
    
    saveRDS(first_exon_data, outpath)
}
