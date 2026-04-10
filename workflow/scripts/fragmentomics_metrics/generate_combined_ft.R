suppressPackageStartupMessages(library(tidyverse))

get_top_15k_by_variance <- function(df) {
  output <- sapply(df %>% dplyr::select(-c(sample)), var) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "feature") %>% 
    dplyr::rename(variance = 2) %>% 
    arrange(desc(variance)) %>% 
    head(n = 15000) %>% 
    pull(feature)
  return(output)
}

append_suffix_to_colnames <- function(df, suffix) {
  output <- df %>% 
    pivot_longer(cols = -c(sample), names_to = "feature", values_to = "value") %>% 
    mutate(feature = str_c(feature, suffix, sep = "__"))
  return(output)
}

print("Generating Combined Feature Table")

outpath <- snakemake@output[["combined_rds"]]

dfs <- list()

# Only iterate over named inputs to avoid technical crashes with empty Snakemake indices
input_names <- names(snakemake@input)[names(snakemake@input) != ""]

if (length(input_names) == 0) {
    stop("FATAL: No feature tables provided to combine. Check upstream rules.")
}

for (suffix in input_names) {
  f <- snakemake@input[[suffix]]
  
  # We removed the file.exists check to ensure read_rds fails loudly if a required file is missing
  data <- read_rds(f)
  
  if (nrow(data) == 0) {
      stop(paste("FATAL: Feature table is empty for category:", suffix))
  }
  
  dfs[[suffix]] <- data %>% append_suffix_to_colnames(suffix = suffix)
}

combined_df <- bind_rows(dfs) %>% 
  pivot_wider(names_from = feature, values_from = value)

if (nrow(combined_df) == 0) {
    stop("FATAL: Combined feature table ended up with 0 rows.")
}

top_features <- get_top_15k_by_variance(combined_df)

combined_df_filtered <- combined_df %>% dplyr::select(sample, all_of(top_features))

saveRDS(combined_df_filtered, outpath)
