library(dplyr)
library(ape)

## Data preparation
format_sample_ids <- function(sample_vector){
  formatted_vector <- sample_vector %>% gsub("\\.","_",.)%>% gsub("-|\\.","_",.) %>% 
    gsub("tasmul_|_[0-9].*","",.) %>% gsub("^X", "", .)
  return(formatted_vector)
}

remove_rare_features <- function(feat_table, 
                                 prevalence_cutoff = 0.1, 
                                 avg_abundance_cutoff = 0.005,
                                 rescale = TRUE) {
  
  if (prevalence_cutoff < 0 | prevalence_cutoff > 1) 
    error('Provided an invalid prevalence_cutoff value')
  if (avg_abundance_cutoff < 0 | avg_abundance_cutoff > 1) 
    error('Provided an invalid avg_abundance_cutoff value')
  
  # Initialize table with all current features
  new_feat_table <- feat_table
  n_samples <- nrow(feat_table)
  n_taxa_before_filter <- ncol(feat_table)
  
  # Prevalence calculations (number of non-zero values per feature)
  prevalence_min_samples <- n_samples * prevalence_cutoff
  
  prevalence_mask <- colSums(new_feat_table > 0) >= prevalence_min_samples
  prevalent_features <- names(prevalence_mask)[prevalence_mask]
  new_feat_table <- new_feat_table[,prevalent_features]
  print(length(prevalent_features))
  
  # Average abundance calculations
  abundance_mask <- (colSums(new_feat_table) / n_samples) >= avg_abundance_cutoff
  abundant_features <- names(abundance_mask)[abundance_mask]
  print(length(abundant_features))
  
  # Keep only genera with sufficient prevalence and abundance
  target_features <- intersect(prevalent_features, abundant_features)
  n_taxa_after_filter <- length(target_features)
  message(n_taxa_before_filter-n_taxa_after_filter, ' of ',n_taxa_before_filter, ' features were removed.')
  new_feat_table <- new_feat_table[, target_features]
  
  # Rescale the table so each sample (row) sums to 1
  if (rescale) {
    new_feat_table <- new_feat_table / rowSums(new_feat_table)
  }
  return(new_feat_table)
}