library(glmnet)
library(rsample)
library(pROC)
library(dplyr)
library(permuco)
library(emmeans)
library(tibble)
library(vegan)

perm_ancova_pairs <- function(df, X, group, covariates, np = 2000,
                              emmeans_adjust = "tukey",  # multiplicity control for parametric pairwise p
                              method = "freedman_lane") {
  
  d <- df
  d[[group]] <- gsub("-", " ",d[[group]]) %>% as.factor() 
  for (cn in covariates) if (is.character(d[[cn]])) d[[cn]] <- as.factor(d[[cn]])
  
  # Model formula: X ~ group + covariates
  form <- as.formula(paste(X, "~", paste(c(group, covariates), collapse = " + ")))
  fit_perm <- aovperm(
    formula = form,
    data    = d,
    np      = np,
    method  = method
  )
  s <- summary(fit_perm)
  p_perm_omnibus <- s[rownames(s) == group, "resampled P(>F)"]
  
  # Parametric LM for EMMs (estimated marginal means)
  fit_lm <- lm(form, data = d)
  emm <- emmeans(fit_lm, as.formula(paste0("~", group)))
  
  pair_emm <- contrast(emm, method = "pairwise", adjust = emmeans_adjust)
  pair_emm_df <- as.data.frame(pair_emm) %>%
    rename(estimate = estimate, se = SE, df = df, t = t.ratio, p_emm = p.value) %>%
    mutate(feature = X) %>%
    select(feature, contrast, estimate, se, df, t, p_emm)
  
  # Calculate pairwise effect size
  split_contrast <- function(lbl) {
    parts <- strsplit(lbl, " - ", fixed = TRUE)[[1]]
    list(a = parts[1], b = parts[2])
  }
  
  pair_perm_list <- lapply(pair_emm_df$contrast, function(lbl) {
    parts <- split_contrast(lbl)
    subd <- d %>% filter(.data[[group]] %in% c(parts$a, parts$b)) %>%
      droplevels()
    
    # two-level ANCOVA with same covariates
    fit_perm_pair <- aovperm(
      formula = form,
      data    = subd,
      np      = np,
      method  = method
    )
    s_pair <- summary(fit_perm_pair)
    p_pair <- s_pair[rownames(s_pair) == group, "resampled P(>F)"]
    
    # Effect size standardization using residual SD from parametric LM on the same subset
    fit_lm_pair <- lm(form, data = subd)
    # get the adjusted mean difference that matches the direction in emmeans
    emm_pair <- emmeans(fit_lm_pair, as.formula(paste0("~", group)))
    cont_pair <- as.data.frame(contrast(emm_pair, list(this = c(ifelse(levels(subd[[group]])[1] == parts$a, 1, -1),
                                                                ifelse(levels(subd[[group]])[2] == parts$b, -1, 1))), adjust = "none"))
    est <- cont_pair$estimate[1]
    sigma_hat <- sigma(fit_lm_pair)
    df_res    <- df.residual(fit_lm_pair)
    J         <- 1 - 3 / (4*df_res - 1)
    cohen_d   <- est / sigma_hat
    hedges_g  <- J * cohen_d
    
    data.frame(contrast = lbl,
               p_value_perm_pair = p_pair,
               effect_size = est,
               cohen_d = cohen_d,
               hedges_g = hedges_g,
               stringsAsFactors = FALSE)
  })
  pair_perm_df <- bind_rows(pair_perm_list)
  pairwise_results <- pair_emm_df %>%
    left_join(pair_perm_df, by = "contrast")
  return(pairwise_results)
}


## Wilcox test between pairs of values in a character column in a df 
add_significance <- function(df, pval_col = "p_value"){
  df$fdr_p_value <- p.adjust(df[[pval_col]], method = "fdr")
  df <- add_sigstars(df,"fdr_p_value")
  return(df)
}

add_sigstars <- function(df, fdr_p_value = "fdr_p_value") {
  df$significance <- case_when(
    df[fdr_p_value] < 0.001 ~ "***",
    df[fdr_p_value] < 0.01  ~ "**",
    df[fdr_p_value] < 0.05  ~ "*",
    TRUE ~ ""
  )
  return(df)
}

wilcox_test_by_pair <- function(df, group_var, formula, adjust_p_value = TRUE){
  group_pairs <- combn(unique(df[[group_var]]), 2, simplify = FALSE)
  results_df <- data.frame(group1 = character(),
                           group2 = character(),
                           p_value = numeric())
  
  # Perform Mann-Whitney test in a for loop
  for (pair in group_pairs) {
    test_data <- df[df[[group_var]] %in% pair,]
    test_result <- wilcox.test(as.formula(formula), data = test_data)
    
    # Store results in the dataframe
    results_df <- rbind(results_df, data.frame(
      group1 = pair[1],
      group2 = pair[2],
      p_value = test_result$p.value
    ))
  }
  if (adjust_p_value) results_df <- results_df %>% add_significance()
  
  return (results_df)
}

## Helper functions for regressing confounders

remove_factors <- function(df){
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) {
      if (all(!is.na(as.numeric(as.character(x))))) {
        as.numeric(as.character(x))
      } else if (all(levels(x) %in% c("TRUE", "FALSE"))) {
        as.logical(as.character(x))
      } else {
        as.character(x)
      }
    } else {
      x
    }
  })
  return(df)
}


get_residuals <- function(df, regressors){
  residuals_df <- df
  target_cols <- names(which(sapply(residuals_df, is.numeric))) %>% setdiff(.,regressors)
  for (x_var in target_cols){
    formula_string <- paste(x_var, "~", paste(regressors, collapse = " + "))
    model <- lm(formula = as.formula(formula_string), data = df)
    residuals_df[x_var] <- residuals(model)
  }
  return(residuals_df)
}

test_residuals <- function(df, regressors, x_var, y_var, show_plots = F, method = "wilcoxon"){
  # set strings as factors
  df[] <- lapply(df, function(x) if (is.character(x)) as.factor(x) else x)
  
  # Calculate the residuals between the regressors and the feature of interest
  formula_string <- paste(x_var, "~", paste(regressors, collapse = " + "))
  model <- lm(formula = as.formula(formula_string), data = df)
  df$residuals <- residuals(model)
  
  ## Helper plots to interpret the model
  
  if (show_plots) {
    added_vars <- car::avPlots(model)
    imp <- relaimpo::calc.relimp(model, type = "lmg", rela = TRUE)
    barplot(imp$lmg, main = "Relative Importance (LMG)", las = 2)
  }
  # Test and visualize the differences in residuals between the groups defined by y_var
  test_result <- if (method == "wilcoxon") {
    wilcox.test(df[["residuals"]] ~ df[[y_var]], exact = FALSE)
  } else {
    t.test(df[["residuals"]] ~ df[[y_var]])
  }
  
  #wilcox <- wilcox.test(df[["residuals"]] ~ df[[y_var]], exact = FALSE)
  ggp <- ggplot(df, aes_string(y_var, "residuals", color = y_var, fill = y_var)) +
    plot_box_and_jitter() + labs(y = "Residuals")
  return(list(test_result$p.value, ggp, df))
}

# Logistic regression and ROC curves 

run_glm_prediction <- function(feat_table, labels,
                               n_folds = 5, n_repeats = 5, n_shuffles = 10,
                               seed = 111, alpha = 1, scale_features = F) {
  set.seed(seed)
  # Setup output data frames
  auc_df <- data.frame()
  roc_df <- data.frame()
  feature_importance_df <- data.frame()
  
  # Setup CV folds and data
  folds <- vfold_cv(data.frame(id = 1:nrow(feat_table)), v = n_folds, repeats = n_repeats)
  X <- as.matrix(feat_table)
  y <- as.factor(labels)
  
  # Evaluate shuffled labels
  shuffled_res <- list()
  for (i in 0:n_shuffles) {
    y_for_pred <- if (i == 0) y else sample(y)
    eval_result <- evaluate_roc(X, y_for_pred, folds, alpha, scale_features)
    feature_importance_df <- rbind(feature_importance_df, eval_result$feature_importance_df %>% mutate(shuffle = i))
    roc_result <- eval_result$roc_df
    auc_df <- rbind(auc_df, list(i, mean(eval_result$auc_vec), sd(eval_result$auc_vec)))
    roc_df <- rbind(roc_df, roc_result %>% mutate(shuffle = i)) 
  }
  colnames(auc_df) <- c("shuffle", "mean_auc", "sd_auc")
  return(list(auc_df = auc_df, roc_df = roc_df, feature_importance_df = feature_importance_df))
}


evaluate_roc <- function(X, y, folds,alpha = 1, scale_features = F) {
  library(glmnet)
  library(pROC)
  
  predictions <- data.frame()
  feature_importances <- list()
  auc_vec <- c()
  
  for (i in seq_along(folds$splits)) {
    split <- folds$splits[[i]]
    train_idx <- split$in_id
    test_idx <- setdiff(1:nrow(X), train_idx)
    X_train <- X[train_idx, ] ; y_train <- y[train_idx]
    X_test <- X[test_idx, ] ; y_test <- y[test_idx]
    
    if (scale_features){
      X_scaled <- min_max_scale_features(X_train, X_test)
      X_train <- X_scaled[["train"]]
      X_test <- X_scaled[["test"]]
    }
    model <- cv.glmnet(x = X_train, y = y_train, family = "binomial", alpha = alpha)
    # Get non-zero coefficients at lambda.min
    coef_i <- as.matrix(coef(model, s = "lambda.min"))
    coef_i <- coef_i[-1, , drop = FALSE]  # remove intercept
    feature_importances[[i]] <- as.numeric(coef_i)
    
    pred <- predict(model, X_test, s = "lambda.min", type = "response")[, 1]
    predictions_local <- data.frame(pred = pred, actual = y_test)
    predictions <- rbind(predictions, predictions_local)
    
    # Add AUC of current split to the auc vector
    auc_local <- as.numeric(auc(roc(predictions_local$actual,predictions_local$pred, quiet = TRUE)))
    auc_vec <- c(auc_vec, auc_local)
  }
  
  # Compute average feature importance
  feature_importances_mat <- do.call(rbind, feature_importances)
  feature_importance_df <- data.frame(
    feature = colnames(X),
    selection_rate =  colMeans(feature_importances_mat >0),
    mean_coef = colMeans(feature_importances_mat),
    sd_coef = apply(feature_importances_mat, 2, sd))
  
  # Compute ROC with CI
  roc_obj <- roc(predictions$actual, predictions$pred, quiet = TRUE, ci = TRUE, ci.alpha = 0.95)
  
  ci_obj <- ci.se(roc_obj, specificities = seq(0, 1, length.out = 25), boot.n = 100)
  roc_df <- data.frame(
    specificity = as.numeric(rownames(ci_obj)),
    lower = ci_obj[, 1],
    median = ci_obj[, 2],
    upper = ci_obj[, 3]
  )

 return(list(
   roc_df = roc_df,
   auc_vec = auc_vec, 
   feature_importance_df = feature_importance_df))
}


# Scaling helper: learns min/max from training, applies to both train/test, clamps to [0,1]
min_max_scale_features <- function(X_train, X_test) {
  mins <- apply(X_train, 2, min, na.rm = TRUE)
  maxs <- apply(X_train, 2, max, na.rm = TRUE)
  ranges <- maxs - mins
  to_scale <- names(which(ranges > 0 & (mins < 0 | maxs > 1)))
  
  if (length(to_scale) > 0) {
    # Apply scaling
    scale_cols <- function(mat) {
      mat[, to_scale] <- sweep(
        sweep(mat[, to_scale, drop = FALSE], 2, mins[to_scale], "-"),
        2, ranges[to_scale], "/"
      )
      mat[, to_scale] <- pmin(pmax(mat[, to_scale], 0), 1)
      mat
    }
    X_train <- scale_cols(X_train)
    X_test  <- scale_cols(X_test)
  }
  
  list(train = X_train, test = X_test)
}

## ORA
# Basic Over-Representation Analysis (ORA) function
ora <- function(target, pathway, universe = NULL) {
  # If universe not provided, use union of all KOs
  if (is.null(universe)) {
    universe <- unique(c(target, pathway))
  }
  
  a <- length(intersect(target, pathway))          # in target & pathway
  b <- length(setdiff(target, pathway))            # in target only
  c <- length(setdiff(pathway, target))            # in pathway only
  d <- length(universe) - (a + b + c)                  # neither in target nor pathway
  mat <- matrix(c(a, c, b, d), nrow = 2)
  
  test <- fisher.test(mat, alternative = "greater")
  result <- list(
    p_value = test$p.value,
    odds_ratio = unname(test$estimate),
    pathway_size = length(pathway),
    n_intersect = a,
    n_target_only = b,
    n_in_pathway_only = c,
    intersect_vec = intersect(target, pathway))
  return(result)
}

run_ora_on_pathway_list <- function(target, pathway_list, universe){
  ora_results <- lapply(names(pathway_list), function(pw) {
    pw_comp <- pathway_list[[pw]]
    res <- ora(target, pw_comp, universe)
    intersect_vec <- ifelse(length(res$intersect_vec)>0, list(res$intersect_vec), NA)
    
    data.frame(
      pathway    = pw,
      p_value    = res$p_value,
      odds_ratio = as.numeric(res$odds_ratio),
      pathway_size = res$pathway_size,
      n_intersect = res$n_intersect,
      n_target_only = res$n_target_only,
      n_in_pathway_only = res$n_in_pathway_only
    ) %>% mutate(intersect_vec = intersect_vec)
  })
  ora_results_df <- do.call(rbind, ora_results)
  return(ora_results_df)
}

## Permanova
get_permanova <- function(disct_metadata, 
                          dist_matrix, 
                          dist_metric_name,
                          sample_id_column = 'sample_id', 
                          metadata_cols_to_analyze = NULL, 
                          metadata_cols_fix_effects = NULL,
                          fdr_threshold = 0.1,
                          permutations = 9999,
                          plot_labels = NULL,
                          pcoa_fill = NULL,
                          ord = NULL,
                          add_centroid = TRUE,
                          seed = 111,
                          level = 0.9,
                          pairs_mode = FALSE
){
  
  # Required libraries loading
  set.seed(seed)
  
  # Define metadata_cols_to_analyze as all columns if not provided, and check if exist if provided.
  if (!is.null(metadata_cols_to_analyze)) {
    if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  } else {
    metadata_cols_to_analyze <- colnames(disct_metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.\n")
  }
  
  # More data verifications
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in metadata")
  if (!all(metadata_cols_to_analyze %in% colnames(disct_metadata)))
    stop("One of metadata_cols_to_analyze not found in metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the metadata. Only samples with both features and metadata will be used for the following analysis.\n")
  
  # Make sure all random/fixed effects are indeed in data
  if (!is.null(metadata_cols_fix_effects)) {
    if (!all(metadata_cols_fix_effects %in% colnames(disct_metadata)))
      stop("One of metadata_cols_fix_effects not found in metadata")
  }
  
  # Merge the dataframes
  df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE)  %>%
    remove_rownames %>% 
    column_to_rownames(var=sample_id_column)
  
  if (nrow(df) == 0) stop('No common samples found.')
  
  # Add fixed effects to a string that later will be added to the formula and check that they are not one of the metadata_cols_to_analyze.
  mixed_string = ""
  if (!is.null(metadata_cols_fix_effects)) {
    for (fixed_col in metadata_cols_fix_effects) {
      # Validate that all metadata_cols_fix_effects not in metadata_cols_to_analyze.
      if (fixed_col %in% metadata_cols_to_analyze)
        stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                    fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
      
      # Remove rows with missing fixed effect values
      df <- df %>% filter(!is.na(get({{fixed_col}})))
      mixed_string =  paste0(mixed_string, " + tmp_df[,'", fixed_col, "'] ")
    }
  }

  # Run adonis2
  pvals <- c()
  for (col in metadata_cols_to_analyze) {
    message('Working on variable: ', col)
    tmp_df <- df[!is.na(df[,col]), ]
    
    # Check that `col` is not a constant
    if (n_distinct(tmp_df[[col]]) == 1) stop(paste('Column', col, 'has only one value, cannot run PERMANOVA'))
    
    formula_string <- paste("tmp_df[, rownames(tmp_df)] ~ tmp_df[[col]]", mixed_string)
    permanova <- adonis2(as.formula(formula_string),
                         permutations = permutations)
    pvals[[col]] <- permanova$`Pr(>F)`[1]
  }
  
  # Define pcoa_values varible as NULL, in order to return it when there is no significant results.
  pcoa_values <- NULL
  
  # Call plot_colored_pcoa_discrete for significant features.
  all_fdrs <- p.adjust(pvals, 'fdr')
  plots <- list()
  if(!all(all_fdrs > fdr_threshold)) {
    for (metadata_feature in names(all_fdrs[all_fdrs <= fdr_threshold])){
      
      # Generate PCoA plots of all the data
      if (!pairs_mode){
        p12 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, plot_labels, pcoa_fill, ord, add_centroid, level)
        plots[[paste0(metadata_feature,"_pc1&pc2")]] <- p12$p
        pcoa_values <- p12$pcoa_values
        p13 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, plot_labels, pcoa_fill, ord, add_centroid, level)$p
        plots[[paste0(metadata_feature,"_pc1&pc3")]] <- p13
        p23 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, plot_labels, pcoa_fill, ord, add_centroid, level)$p
        plots[[paste0(metadata_feature,"_pc2&pc3")]] <- p23
        
        # Generate PCoA plots for all pairs of comparisons (although the permanova is for all the data),
      } else {
        disct_metadata[, "metadata_feature"] <- disct_metadata[, metadata_feature]
        pairs <- combn(unique(disct_metadata$metadata_feature), 2)
        for (i in 1:ncol(pairs)) {
          pair <- pairs[, i]
          p12 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, plot_labels, pcoa_fill, ord, add_centroid, level, pair)
          plots[[paste0(metadata_feature, "_pc1&pc2_", pair[1], "_", pair[2])]] <- p12$p
          pcoa_values <- p12$pcoa_values
          p13 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, plot_labels, pcoa_fill, ord, add_centroid, level, pair)$p
          plots[[paste0(metadata_feature, "_pc1&pc3_", pair[1], "_", pair[2])]] <- p13
          p23 <- plot_colored_pcoa_discrete(dist_matrix, disct_metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, plot_labels, pcoa_fill, ord, add_centroid, level, pair)$p
          plots[[paste0(metadata_feature, "_pc2&pc3_", pair[1], "_", pair[2])]] <- p23
          
        }
      }
    }
  }
  
  # Return a list of a data frame of the p-values and the FDRs, and the PCoA plots.
  results = data.frame(unlist(pvals), unlist(all_fdrs)) %>%
    rename(pvals = 1, FDR = 2) %>%
    tibble::rownames_to_column("feature")
  
  return(list(results = results, plots = plots, pcoa_values = pcoa_values))
}

## PCoA colored by a discrete metadata column (called by get_permanova)
# Adapted from DAP_microbiome/analysis/categorical_metadata_analysis.R: added `level` (ellipse confidence level) and
# `pair` (restrict to two groups, used when pairs_mode = TRUE), and returns list(p, pcoa_values) as get_permanova expects.
plot_colored_pcoa_discrete <- function(dist_matrix,
                                       disct_metadata,
                                       sample_id_column,
                                       var_to_color_by,
                                       dist_metric_name,
                                       axis_x = 1,
                                       axis_y = 2,
                                       plot_labels = NULL,
                                       pcoa_fill = NULL,
                                       ord = NULL,
                                       add_centroid = TRUE,
                                       level = 0.95,
                                       pair = NULL){

  # Validation
  if(! var_to_color_by %in% names(disct_metadata))
    stop(var_to_color_by, " column is missing from metadata")
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(disct_metadata))
    stop("sample_id_column not found in disct_metadata")
  if (nrow(disct_metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the disct_metadata. Only samples with features and disct_metadata will be used for the following analysis.\n")

  # Keep only the two groups of the requested pair
  if (!is.null(pair))
    disct_metadata <- disct_metadata[disct_metadata[[var_to_color_by]] %in% pair, ]

  df <- merge(x = dist_matrix, y = disct_metadata, by = sample_id_column, all.x=FALSE)
  df <- df %>% tidyr::drop_na(all_of(var_to_color_by))

  # Compute PCoA using the ape package
  dist_matrix <- as.dist(df[df[,sample_id_column]])
  PCOA <- ape::pcoa(dist_matrix)
  pcoa_values <- PCOA$vectors[,1:10] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(sample_id_column)

  # Format titles for the axes
  ax_x_title <- paste0('Axis ', axis_x, ' (', round(PCOA$values$Relative_eig[axis_x]*100,1) , '% var. explained)')
  ax_y_title <- paste0('Axis ', axis_y, ' (', round(PCOA$values$Relative_eig[axis_y]*100,1) , '% var. explained)')

  # Add metadata to pcoa values, for plotting
  pcoa_values <- pcoa_values %>%
    left_join(disct_metadata, by = sample_id_column)

  # Convert discrete values to given plot labels by plot_labels.
  if (! is.null(plot_labels)){
    for (val in names(plot_labels))
      pcoa_values[, var_to_color_by][pcoa_values[, var_to_color_by] == strtoi(val)] <- plot_labels[[val]]
  }

  if (!is.null(ord))
    pcoa_values[, var_to_color_by] <- factor(pcoa_values[, var_to_color_by] , levels=ord)

  pcoa_values[, "var_to_color_by"] <- pcoa_values[, var_to_color_by]

  # Plot
  p <- ggplot(pcoa_values,
              aes_string(x = paste0('Axis.', axis_x),
                         y = paste0('Axis.', axis_y),
                         fill = var_to_color_by)) +
    geom_point(color = 'black', alpha = 0.8, shape = 21, size = 2) +
    ggtitle(paste0(dist_metric_name, ' ~ ', var_to_color_by)) +
    xlab(ax_x_title) +
    ylab(ax_y_title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(legend.title = element_text(size = 8))

  if (add_centroid){
    p <- p + stat_ellipse(aes_string(x = paste0('Axis.', axis_x),
                                     y = paste0('Axis.', axis_y),
                                     fill = var_to_color_by),
                          geom = "polygon",
                          level = level,
                          alpha = 0.08)

    means <- pcoa_values %>%
      group_by(var_to_color_by) %>%
      summarise(mean_x = mean(!!sym(paste0('Axis.', axis_x))), mean_y = mean(!!sym(paste0('Axis.', axis_y))))
    means[, paste0('Axis.', axis_x)] <- means[, "mean_x"]
    means[, paste0('Axis.', axis_y)] <- means[, "mean_y"]
    means[, var_to_color_by] <- means[, "var_to_color_by"]

    p <- p + geom_point(data = means,
                        aes_string(x = paste0('Axis.', axis_x),
                                   y = paste0('Axis.', axis_y),
                                   fill = var_to_color_by),
                        size = 5, shape =  21) +
      geom_point(data = means,
                 aes_string(x = paste0('Axis.', axis_x),
                            y = paste0('Axis.', axis_y),
                            fill = var_to_color_by),
                 size = 4, shape =  4)
  }

  if (!is.null(pcoa_fill))
    p <- p + scale_fill_manual(values = pcoa_fill)
  return(list(p = p, pcoa_values = pcoa_values))
}


get_differential_abundance <- function(feat_table, 
                                       disct_metadata, 
                                       sample_id_column = 'sample_id', 
                                       metadata_cols_to_analyze = NULL,
                                       metadata_cols_fix_effects = NULL,
                                       metadata_cols_rand_effects = NULL,
                                       fdr_threshold = 0.1,
                                       output_dir = "MaAsLin2_output",
                                       return_maaslin_plots = FALSE,
                                       keep_maaslin_plots = FALSE,
                                       normalization_method = "NONE",
                                       transformation_method = "LOG",
                                       analysis_method = "LM",
                                       category_name_map = NULL,
                                       skip_fdr_correction = FALSE,
                                       logscale_y_for_plot = FALSE,
                                       quiet = FALSE,
                                       ord=NULL) {
  # Required libraries loading
  require(Maaslin2)
  require(ggplot2)
  require(ggsignif)
  require(grid)
  require(png)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Data verifications
  # In case maaslin plots required, change keep_maaslin_plots to TRUE.
  if (return_maaslin_plots & !(keep_maaslin_plots)){
    message("Can not return maaslin plots if keep_maaslin_plots = FALSE. keep_maaslin_plots was changed to TRUE.\n")
    keep_maaslin_plots = TRUE
  }
  
  # If no metadata_cols_to_analyze provided, run the MaAsLin2 on all the metadata columns.
  if (is.null(metadata_cols_to_analyze))
    metadata_cols_to_analyze = colnames(disct_metadata)
  
  # Check if metadata columns are discrete (1, 0 and NA).
  for (col in metadata_cols_to_analyze){ 
    if (length(unique(disct_metadata[,col])) > 3)
      warning(paste0("Too many values in disct_metadata in the ", col, " column"))
  }
  
  # Check if normalization method fits the data (if only one column: normalization method cannot be 'TSS')
  if ((ncol(feat_table) ==  2) & (normalization_method == "TSS")){
    stop("TSS normlaization is not allowed in case there is only one feature in feat_table. 
         Try again with normalization_method = 'NONE'.")
  }
  
  # Create output directory (will be deleted if not required).
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  disct_metadata <- disct_metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(disct_metadata), row.names(feat_table))
  disct_metadata <- disct_metadata[samples_to_include,]
  feat_table <- feat_table[samples_to_include,,drop=FALSE]
  
  # TODO: warn if some samples were not in both tables
  
  # Check that all the metadata_cols_rand_effects are not in metadata_cols_to_analyze.
  for (rand_col in metadata_cols_rand_effects){
    if (rand_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_rand_effects and metadata_cols_to_analyze: ", 
                  rand_col, "\nTry again with different metadata_cols_rand_effects!\n"))
  }
  
  # Go over fixed effects and perform one-hot encoding if needed
  fixed_cols_final_names <- c()
  for (fixed_col in metadata_cols_fix_effects) {
    
    # Validate that metadata_cols_fix_effects are not in metadata_cols_to_analyze.
    if (fixed_col %in% metadata_cols_to_analyze)
      stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                  fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
    
    # Check if fixed column is numeric with less than 5 unique values or categorical. Perform one-hot encoding in these cases.
    if (((is.numeric(disct_metadata[, fixed_col])) & 
         (n_distinct(disct_metadata[, fixed_col]) < 5)) | 
        (!is.numeric(disct_metadata[, fixed_col]))) {
      
      # Print a message if there are more than 5 categories in the metadata_cols_fix_effects because it might affect the FDR correction.
      if (n_distinct(disct_metadata[[fixed_col]]) > 5)
        message("Too many unique values in ", fixed_col, " column. It might affect the FDR correction and omit some significant results!")
      
      # Perform one-hot encoding for columns if they have more then 2 categories.
      if (n_distinct(disct_metadata[[fixed_col]]) > 2){
        message("Performing one-hot encoding on ", fixed_col)
        fixed_cols_final_names <- c(fixed_cols_final_names, paste(fixed_col, unique(disct_metadata[[fixed_col]]), sep = '_'))
        tmp_metadata_rownames <- rownames(disct_metadata)
        disct_metadata <- disct_metadata %>% mutate(value = 1)  %>% spread(fixed_col, value,  fill = 0, sep = '_')
        rownames(disct_metadata) <- tmp_metadata_rownames
        
      } else { # Categorical with 2 values or less. One-hot is not needed 
        fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
      }
      
    } else { # Numeric case, one-hot is not needed  
      fixed_cols_final_names <- c(fixed_cols_final_names, fixed_col)
    }
  }
  
  # Rum MaAsLin2
  maaslin2_signif_results <- data.frame() # Initialize
  maaslin2_all_results <- data.frame() # Initialize
  
  for (col in metadata_cols_to_analyze) {
    if (!quiet) message('Working on column: ', col)
    
    # For multi-category variables, we run differential abundance separately for each pair of categories
    for (pair in combn(na.omit(unique(disct_metadata[[col]])), 2, simplify = FALSE)) {
      
      # Keep only the two selected values (also drop missing values)
      tmp_metadata <- disct_metadata[(!is.na(disct_metadata[[col]])) & (disct_metadata[[col]] %in% pair), ]
      
      # Get the final list of samples for which DA will be calculated
      samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
      tmp_metadata <- tmp_metadata[samples_to_include,,drop=F]
      tmp_feat_table <- feat_table[samples_to_include,,drop=F]
      
      # Remove constant columns if exist
      tmp_metadata <- tmp_metadata[, apply(tmp_metadata, 2, n_distinct) > 1]
      
      # Fix column names if needed
      colnames(tmp_metadata) <- make.names(colnames(tmp_metadata))
      if (!is.null(metadata_cols_rand_effects))
        metadata_cols_rand_effects <- make.names(metadata_cols_rand_effects)
      
      invisible(capture.output(tmp_results <- Maaslin2(
        input_data = tmp_feat_table, 
        input_metadata = tmp_metadata, 
        output = paste(output_dir, "/", col, sep=""),
        fixed_effects = make.names(c(col, fixed_cols_final_names)),
        random_effects = metadata_cols_rand_effects,
        normalization = normalization_method,
        transform = transformation_method,
        analysis_method = analysis_method,
        plot_scatter = keep_maaslin_plots,
        max_significance = fdr_threshold,
        plot_heatmap = FALSE
      )))
      
      tmp_results$results$pair_a <- as.character(pair[[1]])
      tmp_results$results$pair_b <- as.character(pair[[2]])
      
      # Save significant results only
      maaslin2_signif_results <- bind_rows(
        maaslin2_signif_results,
        tmp_results$results %>% filter(qval <= fdr_threshold)
      )
      
      # Save all results
      maaslin2_all_results <- bind_rows(
        maaslin2_all_results,
        tmp_results$results
      )
    }
  }
  
  # New FDR 
  if (skip_fdr_correction){
    maaslin2_all_results$FDR <- maaslin2_all_results$pval
  } else {
    maaslin2_all_results$FDR <- p.adjust(maaslin2_all_results$pval, 'fdr')
  }
  maaslin2_signif_results <- maaslin2_all_results[maaslin2_all_results$FDR < fdr_threshold,]
  
  # Return significant results only for column in metadata_cols_to_analyze (and not in metadata_cols_fix_effects).
  maaslin2_signif_results <- maaslin2_signif_results[maaslin2_signif_results$metadata %in% metadata_cols_to_analyze, ]
  
  # Generate plots to plots list
  plots <- list()
  
  if (nrow(maaslin2_signif_results) > 0) {
    if (return_maaslin_plots) {
      png_list <- fs::dir_ls(path = output_dir, recurse = TRUE, type = "file", glob = "*.png")
      for (png_dir in png_list){
        print(png_dir)
        plots[[png_dir]] <- rasterGrob(readPNG(png_dir))
      }
    } else {
      
      # Fix feature names to match how maaslin outputs them
      feat_table2 <- feat_table
      names(feat_table2) <- make.names(names(feat_table2))
      
      
      for (col in unique(maaslin2_signif_results$metadata)){
        for (feat in unique(maaslin2_signif_results$feature)){
          
          unique_signif_results <- maaslin2_signif_results[(maaslin2_signif_results$metadata == col) & (maaslin2_signif_results$feature == feat), ]
          if (dim(unique_signif_results)[1] == 0)
            next
          
          
          tmp_metadata <- disct_metadata %>% tidyr::drop_na(col)
          samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table2))
          tmp_metadata <- tmp_metadata[samples_to_include,,drop=FALSE]
          feat_table2 <- feat_table[samples_to_include,,drop=FALSE]
          
          tmp <- data.frame(
            col_vals = tmp_metadata[[col]],
            feat_vals = feat_table2[[feat]]
          )
          
          if (!is.null(category_name_map)) 
            tmp$col_vals <- factor(category_name_map[as.character(tmp$col_vals)],
                                   levels = unname(category_name_map))
          
          y_m <- max(tmp$feat_vals) 
          tmp$col_vals<- as.character(tmp$col_vals)
          
          # Add number of samples in each group to the x-axis ticks.
          xlabs <- list()
          for (val in as.character(unique(tmp$col_vals)))
            xlabs[[val]] = paste0(val,"\n(", dim(tmp[as.character(tmp$col_vals) == val,])[1], ")")
          
          # Sort the tmp dataset based on the specified order
          if (!is_null(ord)){     
            tmp$col_vals <- factor(tmp$col_vals, levels = ord)
            tmp <- tmp[order(tmp$col_vals), ]
          }
          
          p <- ggplot(tmp, aes(x = col_vals, group = col_vals, y = feat_vals)) +
            geom_boxplot(aes(fill = col_vals), color = 'black', alpha = 0.7, outlier.shape = NA) +
            geom_jitter(height = 0, width = 0.1, alpha = 0.5, color = 'black', size = 1) +
            theme_classic() + labs(x = col, y = "Abundance", subtitle = feat) +
            scale_x_discrete(labels=xlabs)
          
          
          unique_signif_results$group1 <- unique_signif_results$pair_a
          unique_signif_results$group2 <- unique_signif_results$pair_b
          
          makeStars <- function(x){
            stars <- c("***", "**", "*")
            vec <- c(0, fdr_threshold/100, fdr_threshold/10 ,fdr_threshold)
            i <- findInterval(x, vec)
            stars[i]
          }
          
          unique_signif_results$FDR <- makeStars(unique_signif_results$FDR)
          
          p <- p + stat_pvalue_manual(unique_signif_results,
                                      y.position = max(tmp$feat_vals),
                                      step.increase = 0.1 ,
                                      label = "FDR")
          
          # Log scale y axis if needed
          if (logscale_y_for_plot) {
            p <- p + yscale("log10", .format = TRUE)
          }
          plots[[paste(col, feat, collapse = ';')]] <- p
        }
      }
    } 
  }
  
  if (! keep_maaslin_plots){
    unlink(output_dir, recursive=TRUE)
  }
  
  
  return(list(sig_results = maaslin2_signif_results, all_results = maaslin2_all_results, plots = plots))
}
