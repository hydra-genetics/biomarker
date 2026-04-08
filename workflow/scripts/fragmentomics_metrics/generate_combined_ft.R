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

for (suffix in names(snakemake@input)) {
  f <- snakemake@input[[suffix]]
  if (file.exists(f)) {
    dfs[[suffix]] <- read_rds(f) %>% append_suffix_to_colnames(suffix = suffix)
  } else {
    warning(str_c("Feature file not found, skipping: ", f))
  }
}

if (length(dfs) == 0) {
    stop("No feature tables found to combine.")
}

combined_df <- bind_rows(dfs) %>% 
  pivot_wider(names_from = feature, values_from = value)

top_features <- get_top_15k_by_variance(combined_df)

combined_df_filtered <- combined_df %>% dplyr::select(sample, all_of(top_features))

saveRDS(combined_df_filtered, outpath)
