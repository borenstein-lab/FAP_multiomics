### General visualization tools
plot_box_and_jitter <- function(palette = NULL,font_size = 12, legend_position = "none",
                                point_size = 2) {
  geom =   list(
    geom_boxplot(alpha = 0.3, outlier.shape = NA), 
    geom_point(position = position_jitterdodge(), size = point_size, alpha = 0.9),
    theme_bw(base_size = font_size),
    theme(legend.position = legend_position) 
  )
  if (!is.null(palette)){
    geom <- c(geom, list(scale_fill_manual(values = palette), scale_color_manual(values = palette)))
  }
  return(geom)
  
}

### Relative abundance functions
plot_abundance_barplot <- function(ggdata, palette, x_var = "patient_id", y_var = "ra", fill_var = "taxon" ){
  ggp <- ggplot(ggdata, aes_string(x = x_var, y = y_var, fill = fill_var)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = palette) +
    labs(x = "Sample", y = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    guides(fill = guide_legend(ncol = 2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.direction = "vertical") +
    theme(legend.text = element_text(size = 9)) +
    theme(legend.key.size = unit(9, "points")) 
  return(ggp)
}

### PCA functions
get_pca <- function(df, n_dims, center = F, scale = F){
  pca_result <- df %>% prcomp(.,center = center, scale = scale)
  pca_df <- pca_result$x %>% as.data.frame() %>% select(all_of(paste0("PC", 1:n_dims))) %>% tibble::rownames_to_column("sample_id")
  eig_vals <- (summary(pca_result)$importance[2,][1:n_dims] %>% unname() * 100) %>% round(.,1)
  loadings <- pca_result$rotation[, paste0("PC", 1:n_dims)] %>% as.data.frame() %>%  tibble::rownames_to_column("pathway")
  return(list(eig_vals, pca_df,loadings))
}

plot_pca<- function(ggdata, x_var, y_var,  eig_vals,color_var = "sample_id"){
  ggp <- ggplot(ggdata,aes_string(paste0("PC", x_var), paste0("PC", y_var), color = color_var)) +
    geom_point(size = 2, alpha = 0.7) + 
    ggtitle(paste0('PCA (Axes ',x_var,' vs. ',y_var,')')) +
    xlab(paste0('PC ',x_var,' (', eig_vals[x_var], '% variance explained)')) +
    ylab(paste0('PC ',y_var,' (', eig_vals[y_var], '% variance explained)')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  return(ggp)
}

### PCoA functions
get_pcoa <- function(dist_df, n_dims){
  PCOA <- dist_df %>% pcoa()
  eig_vals <- round(100* PCOA$values$Relative_eig[1:3],n_dims)
  pcoa_df <- PCOA$vectors[,1:n_dims] %>% as.data.frame() %>% tibble::rownames_to_column('sample_id')
  colnames(pcoa_df) <- gsub("Axis.", "axis_", colnames(pcoa_df))
  return(list(eig_vals, pcoa_df))
}

plot_pcoa<- function(ggdata, x_var, y_var,  eig_vals,metric, color_Var = "sample_id"){
  ggp <- ggplot(ggdata,aes_string(paste0("axis_", x_var), paste0("axis_", y_var), color = color_Var)) +
    geom_point(size = 2, alpha = 0.7) + 
    ggtitle(paste0('PCoA - ',metric,' (Axes ',x_var,' vs. ',y_var,')')) +
    xlab(paste0('Axis ',x_var,' (', eig_vals[x_var], '% variance explained)')) +
    ylab(paste0('Axis ',y_var,' (', eig_vals[y_var], '% variance explained)')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  return(ggp)
  
}
plot_pcoa_from_dist_mat <- function(dist_df, metric_name, axis_x = 1, axis_y = 2, color_Var = "sample_id", n_dims = 3, mdata = NULL, sample_col = "sample_id") {
  pcoa_df <- get_pcoa(dist_df, n_dims)
  eig_vals <- pcoa_df[[1]]
  pcoa_df <- pcoa_df[[2]] 
  if (!is.null(mdata)){
    pcoa_df <- pcoa_df %>% left_join(mdata, by = c('sample_id' = sample_col))
  }
  ggp <- plot_pcoa(pcoa_df, axis_x, axis_y, eig_vals, metric_name, color_Var)
  
}

plot_pcoa_2d <- function(dist_df, core_mdata, palette,  x_var = 1, y_var =2,
                         shape_var = "batch", dist_metric = ""){
  pcoa_df <- get_pcoa(dist_df, 2)
  eig_vals <- pcoa_df[[1]]
  pcoa_df <- pcoa_df[[2]] %>% left_join(core_mdata, by = "sample_id")
  
  multi_sample_patient = pcoa_df %>% count(patient_id) %>% filter(n>1) %>% pull(patient_id)
  pcoa_df <- pcoa_df %>% mutate(patient_cat = ifelse(patient_id %in% multi_sample_patient, patient_id, "single_sample"))
  # Plot PCs 1 and 2 by phenotype
  ggp <- ggplot(pcoa_df,aes_string(paste0("axis_", x_var), paste0("axis_", y_var), color = "phenotype", shape = shape_var)) +
    geom_point(size = 2, alpha = 0.7) + 
    ggtitle(paste0('PCoA ',dist_metric,'(Axes ',x_var,' vs. ',y_var,')')) +
    xlab(paste0('Axis ',x_var,' (', eig_vals[x_var], '% variance explained)')) +
    ylab(paste0('Axis ',y_var,' (', eig_vals[y_var], '% variance explained)')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values = palette) 
  plot(ggp)
  
  
  # Plot PCs 1 and 2 by patient and batch
  ggp <- ggplot(pcoa_df,aes_string(paste0("axis_", x_var), paste0("axis_", y_var), color = "patient_cat", group = "patient_cat", shape = shape_var)) +
    geom_point(data = pcoa_df %>% filter(patient_cat == "single_sample"), size = 3, alpha = 0.5, color = "grey") +
    geom_point(data = pcoa_df %>% filter(patient_cat != "single_sample"), size = 3, alpha = 0.7) + 
    geom_line(data = pcoa_df %>% filter(patient_cat != "single_sample")) + 
    ggtitle(paste0('PCoA ',dist_metric,'(Axes ',x_var,' vs. ',y_var,')')) +
    xlab(paste0('Axis ',x_var,' (', eig_vals[x_var], '% variance explained)')) +
    ylab(paste0('Axis ',y_var,' (', eig_vals[y_var], '% variance explained)')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  plot(ggp)
  
  return(pcoa_df)
}

## MaAslin2 functions
plot_volcano <- function(df, x_var = "coef", y_var = "-log10(FDR)",
                         color_var = "color", label_var = "label",
                         fdr_min = 0.1, coef_min = 1,n_label_overlaps = 5 ){
  ggp <- ggplot(df, aes_string(x = x_var, y = y_var,  color = color_var)) +
    geom_point(size = 2) +
    labs(x = "Effect size [log2(FC)]", y = "-log10(FDR)") + 
    guides(alpha = "none",   color = guide_legend(
      title = "Enrichment group",
      override.aes = aes(label = ""))) + 
    geom_hline(yintercept = -log10(fdr_min), linetype = "dashed", color = "black") +
    geom_vline(xintercept = coef_min, linetype = "dashed", color = "black") + 
    geom_vline(xintercept = -coef_min, linetype = "dashed", color = "black") + 
    ggrepel::geom_text_repel(aes_string(label = label_var),
                             force = 10, max.overlaps = n_label_overlaps,
                             min.segment.length = 0) +
    theme_bw(base_size = 12)# + theme(legend.position = "top")
  return(ggp)
}


prepare_maaslin2_for_volcano <- function(df,metadata_vars,  fdr_min = 0.1, coef_min = 1){
  ggdata <- df %>% filter(metadata %in% metadata_vars) %>% 
    mutate(is_significant = FDR < fdr_min, 
           is_coef = abs(coef) > coef_min,
           is_highlight = is_significant & is_coef,
           label = ifelse(is_highlight, feature, ""),
           color = ifelse(is_highlight & coef > coef_min, pair_a, pair_b),
           color = ifelse(!is_highlight,"Non-Significant", color))
  return(ggdata)
}


# to be deleted
run_analysis_diff_abundance <- function(df, mdata, metadata_cols_to_analyze ,palette = NULL,
                                        metadata_cols_fix_effects = c("gender", "age", "bmi", "is_surgery"),
                                        metadata_cols_rand_effects = c(), 
                                        boxplot_x_order = NULL) {
  # Run Maaslin2 differential abundance analysis
  diff_ra <- get_differential_abundance(df, mdata,
                                        metadata_cols_to_analyze = metadata_cols_to_analyze,
                                        metadata_cols_fix_effects = metadata_cols_fix_effects,
                                        metadata_cols_rand_effects = metadata_cols_rand_effects)
  # Plot volcano
  target_col <- metadata_cols_to_analyze[1]
  ggdata <- prepare_maaslin2_for_volcano(diff_ra$all_results, target_col)
  signif_labs <- ggdata %>% filter(is_highlight) %>% pull(label) %>% unique()
  
  
  ggp_volcano <- plot_volcano(ggdata)
  if (! is.null(palette))  {
    ggp_volcano <- ggp_volcano + scale_color_manual(values = palette)
  }
  
  # plot boxplots of significant differences
  boxplots <- diff_ra$plots
  names(boxplots) <- gsub(paste0(target_col," "), "", names(boxplots)) # remove target column name from plot names
  boxplots <- boxplots[names(boxplots)[names(boxplots) %in% signif_labs]] # only keep significant plots
  plot_list = list()
  #for (ggp in boxplots) {
  for (name in names(boxplots)) {
    ggp <- boxplots[[name]]
    ggp <- ggp + labs(x = "") +  theme_bw(base_size = 10)+theme(legend.position = "none")  
    if (! is.null(palette)) {
      ggp <- ggp + scale_fill_manual(values = palette)
    }
    if (! is.null(boxplot_x_order)) {
      ggp <- ggp + scale_x_discrete(limits = boxplot_x_order)
    }
    plot_list <- c(plot_list, list(ggp))
  }
  names(plot_list) <- names(boxplots)
  return(list("diff_ra" = diff_ra$all_results, "ggp_volcano" = ggp_volcano, "ggp_boxplots" = plot_list))
  
}