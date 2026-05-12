suppressPackageStartupMessages(library(tidyverse))

feature_tables <- c("ATAC_entropy", "depth", "frag_bins", "full_gene_depth", "mds", "se", "small_frags", "TFBS_entropy")

for (ft in feature_tables) {
    
    print(str_c("Building feature table for", ft, sep = " "))
    
    file_list <- snakemake@input[[ft]]
    out_path <- snakemake@output[[ft]]
    
    if (length(file_list) == 0) {
      stop(str_c("FATAL: No input files provided for feature category: ", ft))
    }
    
    all_data <- purrr::map_dfr(file_list, ~ read_tsv(.x, show_col_types = FALSE))
    
    if (nrow(all_data) == 0) {
      message(str_c("WARNING: No data found in input files for feature category: ", ft))
      # Create an empty RDS with expected structure if possible, or just skip
      output_data <- tibble(sample = character())
    } else {
        all_data <- all_data %>% 
          dplyr::rename("sample" = 1, "feature" = 2, "value" = 3) %>% 
          filter(str_detect(feature, "CRLF2", negate = TRUE)) %>%   # this gene ambiguously maps to both X and Y chr so is removed
          filter(str_detect(feature, "P2RY8", negate = TRUE)) %>%       # this gene ambiguously maps to both X and Y chr so is removed
          group_by(sample, feature) %>%
          summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
        
        output_data <- all_data %>% 
          pivot_wider(names_from = feature, values_from = value)
        
        # Check if we still have data after filtering/pivoting
        if (nrow(output_data) == 0) {
          message(str_c("WARNING: Resulting feature table is empty after filtering for category: ", ft))
          output_data <- tibble(sample = unique(all_data$sample))
        }

        # replace NAs with 0
        output_data <- output_data %>%
          mutate(across(everything(), ~replace_na(., 0)))
    }
    
    saveRDS(output_data, out_path)
}
