suppressPackageStartupMessages(library(tidyverse))

feature_tables <- c("ATAC_entropy", "depth", "frag_bins", "full_gene_depth", "mds", "se", "small_frags", "TFBS_entropy")

for (ft in feature_tables) {
    
    print(str_c("Building feature table for", ft, sep = " "))
    
    file_list <- snakemake@input[[ft]]
    out_path <- snakemake@output[[ft]]
    
    if (length(file_list) == 0) {
      warning(str_c("No files provided for feature category: ", ft))
      next
    }
    
    all_data <- purrr::map_dfr(file_list, ~ read_tsv(.x, show_col_types = FALSE))
    
    all_data <- all_data %>% 
      dplyr::rename("sample" = 1, "feature" = 2, "value" = 3) %>% 
      filter(str_detect(feature, "CRLF2", negate = TRUE)) %>%   # this gene ambiguously maps to both X and Y chr so is removed
      filter(str_detect(feature, "P2RY8", negate = TRUE)) %>%       # this gene ambiguously maps to both X and Y chr so is removed
      group_by(sample, feature) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    
    output_data <- all_data %>% 
      pivot_wider(names_from = feature, values_from = value)
    
    # remove or replace NAs with 0
    output_data <- output_data %>%
      mutate(across(everything(), ~replace_na(., 0)))
    
    saveRDS(output_data, out_path)
}
