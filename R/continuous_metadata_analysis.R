###############################################################################
#' The following script contains various functions for the analysis of 
#' microbiome-related feature tables against continuous metadata features.
#'  
#' Last updated: 01/01/2023
#' 
#' Implemented functions:  
#' - get_correlations_to_metadata
#' - get_predictability_regression
#' - plot_colored_pcoa 
#' - plot_multiple_colored_pcoa 
#' - get_mantel
#' - get_mixed_effect
###############################################################################

#' Calculate simple correlations between each feature in the feature table, and
#'    each metadata feature specified. Produce scatter plots with regression
#'    lines for significantly correlated features. 
#' 
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column.
#' @param metadata_cols_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. If NULL, all available 
#'    columns will be analyzed.
#' @param corr_method a character string indicating which correlation coefficient
#'    is to be used for the test. One of "pearson", "kendall", or "spearman", can
#'    be abbreviated.
#' @param fdr_threshold FDR threshold used to define significant findings. Plots
#'    will only be generated for significant findings.
#' @param r_threshold r threshold used to define what plots to return (only
#'    plots with |r| > r_threshold).
#' @param remove_zeros remove samples with zeros from the correaltion of specific
#'    feature (default: FALSE).
#' @param logscale_y log scale of the Y axis in the plots (default: FALSE).
#' @param xlab Manually change the x label title (by default it will be be the 
#'    relevant col metadata_cols_to_analyze.
#' @param color change the regression line color (default: "blue").
#' @param stratify_by_col column in the metadata to stratify the correlations b,
#'    need to be string.
#' @param strat_colors in case 'stratify_by_col', add list with the names of the
#'    unique values in metadata$stratify_by_col and colors as values to change
#'    the regression line color by group. For example:
#'    list("group.A" = "red", ...)
#' 
#' 
#' @return A list with two elements: 
#'    `results` contains statistical test results. Data frame with rows that 
#'     represent features and columns with r, p and adjusted p.
#'    `plots` is a list of relevant plots: correlation plots with r and adjusted 
#'     p as title. Plots generetaed for adjusted p-value < fdr_threshold.
#'
#' @export
get_correlations_to_metadata <- function(feat_table, 
                                         metadata, 
                                         sample_id_column = 'sample_id', 
                                         metadata_cols_to_analyze = NULL, 
                                         corr_method = 'spearman',
                                         fdr_threshold = 0.1,
                                         r_threshold = 0, 
                                         remove_zeros = FALSE,
                                         logscale_y = FALSE,
                                         xlab = NULL,
                                         color = "blue",
                                         stratify_by_col = NULL,
                                         strat_colors=NULL) {
  # Required libraries loading
  require(ggplot2)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Data verifications 
  if (!sample_id_column %in% colnames(feat_table))
    stop("sample_id_column not found in feat_tale")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  
  if (nrow(metadata) != nrow(feat_table))
    warning("Note that nrow(metadata) != nrow(feat_table). Only samples with features and metadata will be used for the following analysis.")
  if (! corr_method %in% c('pearson', 'spearman'))
    stop("Invalid 'corr_method' argument. Only 'pearson' and 'spearman' supported.")

  if (is.null(metadata_cols_to_analyze)){
    metadata_cols_to_analyze <- colnames(metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.")
  } else {
    if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  }
  
  if(!is.null(stratify_by_col)){
    if (!stratify_by_col %in% colnames(metadata))
      stop("stratify_by_col not found in metadata")
    if (n_distinct(metadata[,stratify_by_col]) > 10)
      message(paste0("Note that stratify_by_col: ", stratify_by_col, " have more then 10 distinct values, it might affect your FDR."))
  } else {
    # Add dummy stratification column.
    metadata[,"stratify_by_col"] = "no_strat"
    stratify_by_col = "stratify_by_col"
  }
  
  metadata[, "stratify_by_col"] <- metadata[, stratify_by_col]  
  metadata <- metadata %>% drop_na(stratify_by_col)

  # Merge the metadata and the feature_table by sample_id_column.
  df <- merge(x = feat_table, y = metadata, by = sample_id_column, all.x=FALSE)

  # Generating outputs lists 
  plots <- list()
  pvals <- list()
  rhos <- list()

  safe_corr <- safely(cor.test)
  
  # Main calculations
  for (strat_val in unique(metadata$stratify_by_col)){
    strat_df <- df[df[, stratify_by_col] == strat_val, ]
    for (m_col in metadata_cols_to_analyze){
      for (f_col in setdiff(colnames(feat_table), sample_id_column)){
        tmp_df <- strat_df
        if (remove_zeros)
          tmp_df <- strat_df[strat_df[, f_col] != 0, ]
        res <- safe_corr(as.numeric(unlist(tmp_df[, f_col])), as.numeric(unlist(tmp_df[, m_col])), method=c(corr_method), exact=FALSE)
        if (is.null(res$result))
          message(paste0(res$error, strat_val))
        
        pvals[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)] <- res$result$p.value
        rhos[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)] <- res$result$estimate[[1]]
      }
    }
  }
  all_fdrs <- p.adjust(pvals, 'fdr')

  results = data.frame(unlist(rhos), unlist(pvals), unlist(all_fdrs))
  names(results) = c("rhos", "p-values", "FDR")
  
  # Generate plots (only for significant results after FDR)
  for (strat_val in unique(metadata$stratify_by_col)){
    strat_df <- df[df[, stratify_by_col] == strat_val, ]
    for (m_col in metadata_cols_to_analyze){
      for (f_col in colnames(feat_table[, -which(names(feat_table) == sample_id_column)])){
        if (is.na(all_fdrs[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)]))
          next
        if (all_fdrs[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)] > fdr_threshold)
          next
        
        r <- rhos[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)]
        if (abs(as.numeric(r)) < r_threshold)
          next
        
        tmp_df <- strat_df
        if (remove_zeros)
          tmp_df <- strat_df[strat_df[, f_col] != 0, ]
        
        f_print <- f_col
        if (grepl(";g__", f_col, fixed=TRUE))
          f_print <- paste("f__", strsplit(f_col, ";f__")[[1]][2], sep="")
        
        m_print <- m_col
        if (!is.null(xlab))
          m_print <- xlab
        
        subtitle = ""
        if(strat_val != "no_strat"){
          subtitle =  paste0(strat_val, " (", dim(tmp_df)[1], ")")
          if (!is.null(strat_colors))
          color =  strat_colors[[strat_val]]
        }
        
        title <- paste( m_print, "vs.\n", f_print, 
                       "\nFDR:", round(all_fdrs[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)], 3), " r:", round(rhos[paste0(f_col, "_vs._", m_col, "_strat_", strat_val)][[1]], 3))
        p <- ggscatter(tmp_df, x = m_col, y = f_col, title=title, subtitle = subtitle,
                       color = "black", shape = 19, size = 1,
                       add = "reg.line",
                       add.params = list(color = color, fill = "lightgray"),
                       conf.int = TRUE,
                       cor.coeff.args = list(method = corr_method,  label.sep = "\n"),
                       xlab = m_print , ylab = f_print
        ) + theme(plot.title = element_text(hjust = 0.5)) + font("xlab", size = 10) + font("ylab", size = 10) +
          font("xy.text", size = 10, color = "gray", face = "bold") 
          if (logscale_y)
            p <- p + yscale("log10", .format = TRUE)
          p <- p + theme(plot.subtitle=element_text(size=10, hjust=0.5, color=color))
        
        plots[[paste0(f_col, "_vs._", m_col)]] <- p
      }
    }
  }
  
  results <- as_tibble(results, rownames = "feat.vs.metadata")
  return(list(results = results, plots = plots))
}

#' Train machine learning regression models based on all features in the 
#'    feature table, to assess predictability of each continuous metadata 
#'    feature. Model evaluation is carried out using a Spearman correlation 
#'    between predicted and actual values. 
#'    
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. 
#' @param ml_method One of "glm" (generalized linear model), or "rf" (random 
#'    forest regressor)
#' @param fdr_threshold FDR threshold used to define significant findings 
#' @param corr_threshold Correlation threshold used to define significant findings 
#' @param corr_method Correlation method for testing predicted values vs. true 
#'    values. Either 'pearson' or 'spearman'
#' @param n_ml_repeats Number of repeats for repeated cross-validation
#' @param n_ml_folds Number of folds for repeated cross-validation
#' @param feat_importance For RF models, 'altmann' method can be computed. Use 
#'   'none' to avoid feature importance (and decrease run times). 
#'   Default: 'none'.
#' @param rand_seed Random seed for reproducible results
#' 
#' @return A list with one element: 
#'    `results` contains summary results: correlations between predicted and 
#'    actual values, summarized over cv-folds.
#' @export
get_predictability_regression <- function(feat_table, 
                                          metadata, 
                                          sample_id_column, 
                                          vars_to_analyze, 
                                          ml_method = 'glm',
                                          fdr_threshold = 0.1,
                                          corr_threshold = 0.3,
                                          corr_method = "pearson",
                                          n_ml_repeats = 10,
                                          n_ml_folds = 5,
                                          feat_importance = 'none',
                                          rand_seed = 1111,
                                          quiet = FALSE) {
  # Required libraries loading
  require(rsample)
  require(dplyr)
  require(glmnet)
  require(ranger)
  # require(metap) 
  
  # Argument validations 
  if(! all(vars_to_analyze %in% names(metadata)))
    stop("Invalid 'vars_to_analyze' argument. Some of the variables provided are missing from the metadata table")
  if(is.null(vars_to_analyze))
    stop("Invalid 'vars_to_analyze' argument. Cannot be NULL")
  if(! ml_method %in% c('glm', 'rf'))
    stop("Invalid 'ml_method' argument. Only 'glm' and 'rf' supported.")
  if (! corr_method %in% c('pearson', 'spearman'))
    stop("Invalid 'corr_method' argument. Only 'pearson' and 'spearman' supported.")
  
  # Main calculations
  n_samples <- nrow(feat_table)
  
  # Generate folds for cross-validation
  set.seed(rand_seed)
  folds <- vfold_cv(
    data.frame(sample_num = 1:n_samples),
    v = n_ml_folds,
    repeats = n_ml_repeats
  )
  
  # Collect evaluation stats here
  cv_correlations <- data.frame()
  cv_feat_imp <- data.frame()
  
  # Fix names of features (causes errors in some models)
  names(feat_table) <- make.names(names(feat_table))
  
  # Iterate over variables
  for (measure in vars_to_analyze) {
    if (!quiet) message("Training model for variable: ", measure)
    
    # Prepare temporary table for training (features + label only)
    tmp <- feat_table %>% 
      left_join(
        metadata[,c(sample_id_column, measure)] %>% 
          rename(Label = 2), 
        by = sample_id_column) %>%
      select(-any_of(c(sample_id_column)))
    
    # Cross validation loop
    for (i in 1:nrow(folds)) {
      fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
      
      # Get train/test samples of this fold
      train_samples <- folds$splits[[i]]$in_id
      train_data <- tmp[train_samples,] %>% filter(!is.na(Label))
      test_data <- tmp[-train_samples,] %>% filter(!is.na(Label))
      
      # Train logistic regression / RF
      if (ml_method == 'glm') {
        lasso_model <- cv.glmnet(
          x = train_data %>% select(-Label) %>% as.matrix(), 
          y = train_data$Label, 
          grouped = (n_samples > 50)) # grouped = FALSE relevant for very small sample sets only
      } else if (ml_method == 'rf') {
        rf_model <- ranger(Label ~ ., data = train_data, importance = ifelse(feat_importance=='altmann','permutation','none'), seed = rand_seed)
        if (feat_importance != 'none') {
          feat_imp <- importance_pvalues(rf_model, method="altmann", formula = Label ~ ., data = train_data, seed = rand_seed)
        }
      }
        
      # Make prediction on held out samples
      if (ml_method == 'glm') {
        preds <- predict(lasso_model, 
                         test_data %>% select(-Label) %>% as.matrix(), 
                         s = "lambda.min") %>% unname()
      } else if (ml_method == 'rf') {
        preds <- predict(rf_model, test_data %>% select(-Label))$predictions
      }
      
      # Calculate model performance for this fold
      tmp_eval <- data.frame(preds = preds, true_label = test_data$Label)   
      cv_correlations <- bind_rows(
        cv_correlations,
        data.frame(
          measure = measure,
          fold_id = fold_id,
          test_cor = cor(tmp_eval$preds, tmp_eval$true_label, method = corr_method),
          test_cor_p = cor.test(tmp_eval$preds, 
                                tmp_eval$true_label, 
                                alternative = 'greater', # Check positive association
                                method = corr_method,
                                exact = FALSE)$p.value
        )
      )
      
      if (feat_importance != 'none')
        cv_feat_imp <- bind_rows(
          cv_feat_imp,
          as.data.frame(feat_imp) %>% 
            tibble::rownames_to_column('feature') %>%
            mutate(measure = measure, fold_id = fold_id)
        )
    }
  }
  
  # sumlog2 <- function(ps) {ps <- ps[!is.na(ps)]; return(sumlog(ps)$p)}
  sumlog2 <- function(pvals) {
    X2 <- -2 * sum(log(pvals))
    df <- 2 * length(pvals)
    p_combined <- pchisq(X2, df, lower.tail = FALSE)
    return(p_combined)
  }
  
  # Summarise over folds
  summary_correlations <- cv_correlations %>%
    group_by(measure) %>%
    summarise(n_folds = n(),
              mean_cor = mean(test_cor, na.rm = TRUE),
              sd_cor = sd(test_cor, na.rm = TRUE),
              p = sumlog2(test_cor_p),
              .groups = "drop") %>%
    mutate(fdr = p.adjust(p)) %>%
    mutate(Significant = fdr < fdr_threshold & mean_cor > corr_threshold)
  
  # Summarize feature importance (if requested by user)
  summary_feat_importance <- data.frame()
  if (feat_importance != 'none')
    summary_feat_importance <- cv_feat_imp %>%
      group_by(measure, feature) %>%
      summarise(n_folds = n(),
                mean_feat_importance = mean(importance, na.rm = TRUE),
                sd_feat_importance = sd(importance, na.rm = TRUE),
                p = sumlog2(pvalue),
                .groups = "drop") %>%
      mutate(fdr = p.adjust(p))
  
  if (! any(summary_correlations$Significant)) {
    message('No significantly well-predicted metadata measures')
  } else {
    message('The following measures were well-predicted by the microbiome using a ', toupper(ml_method), ' model:')
    message(paste(summary_correlations$measure[summary_correlations$Significant], collapse = ', '))
  }
  
  return(list(results = summary_correlations, feat_importance = summary_feat_importance))
}

# TODO Tal: update the documentation.
#' Train machine learning regression models based on all features in the 
#'    feature table, to assess predictability of each continuous metadata 
#'    feature. Model evaluation is carried out using a Spearman correlation 
#'    between predicted and actual values. TODO: EFRAT
#'    
#' @param feat_table A feature table (taxa abundances, pathways, alpha-diversity 
#'    metrics, etc.) with rows as samples and columns as features. A 'sample id' 
#'    column is expected in order to merge the table with the metadata.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. 
#' @param ml_method One of "glm" (generalized linear model), or "rf" (random 
#'    forest regressor)
#' @param fdr_threshold FDR threshold used to define significant findings 
#' @param corr_threshold Correlation threshold used to define significant findings 
#' @param corr_method Correlation method for testing predicted values vs. true 
#'    values. Either 'pearson' or 'spearman'
#' @param n_ml_repeats Number of repeats for repeated cross-validation
#' @param n_ml_folds Number of folds for repeated cross-validation
#' @param rand_seed Random seed for reproducible results
#' 
#' @return A list with one element: 
#'    `results` contains summary results: correlations between predicted and 
#'    actual values, summarized over cv-folds.
#' @export
get_predictability_regression_stratify <- function(feat_table, 
                                          metadata, 
                                          sample_id_column, 
                                          vars_to_analyze, 
                                          ml_method = 'glm',
                                          fdr_threshold = 0.1,
                                          corr_threshold = 0.3,
                                          corr_method = "pearson",
                                          n_ml_repeats = 10,
                                          n_ml_folds = 5,
                                          rand_seed = 1111,
                                          quiet = FALSE,
                                          stratify_by_col = NULL) {
  # Required libraries loading
  require(rsample)
  require(dplyr)
  require(glmnet)
  require(ranger)
  require(metap) 
  require(xgboost)
  
  # Argument validations 
  if(! all(vars_to_analyze %in% names(metadata)))
    stop("Invalid 'vars_to_analyze' argument. Some of the variables provided are missing from the metadata table")
  if(is.null(vars_to_analyze))
    stop("Invalid 'vars_to_analyze' argument. Cannot be NULL")
  if(! ml_method %in% c('glm', 'rf', 'lgbm', 'xgboost'))
    stop("Invalid 'ml_method' argument. Only 'glm', 'rf', 'lgbm' or 'xgboost' supported.")
  if (! corr_method %in% c('pearson', 'spearman'))
    stop("Invalid 'corr_method' argument. Only 'pearson' and 'spearman' supported.")
  
  if(!is.null(stratify_by_col)){
    if (!stratify_by_col %in% colnames(metadata))
      stop("stratify_by_col not found in metadata")
    if (n_distinct(metadata[,stratify_by_col]) > 10)
      message(paste0("Note that stratify_by_col: ", stratify_by_col, " have more then 10 distinct values, it might affect your FDR."))
  } else {
    # Add dummy stratification column.
    metadata[,"stratify_by_col"] = "no_strat"
    stratify_by_col = "stratify_by_col"
  }
  
  metadata[, "stratify_by_col"] <- metadata[, stratify_by_col]  
  metadata <- metadata %>% drop_na(stratify_by_col)
  
  # Fix names of features (causes errors in some models)
  names(feat_table) <- make.names(names(feat_table))
  
  # Collect evaluation stats here
  cv_correlations <- data.frame()
  
  preds_by_strat <- list()
  
  for (strat_val in unique(metadata$stratify_by_col)){
  
    strat_metadata <- metadata %>% filter(stratify_by_col == strat_val)
    
    # Generate strat_df only for checking number of rows 
    strat_df <- feat_table %>% inner_join(strat_metadata, by = sample_id_column)
    n_samples <- nrow(strat_df)
    
    # Main calculations
    
    # Generate folds for cross-validation
    set.seed(rand_seed)
    folds <- vfold_cv(
      data.frame(sample_num = 1:n_samples),
      v = n_ml_folds,
      repeats = n_ml_repeats
    )
    preds_by_measure <- list()
    
    # Iterate over variables
    for (measure in vars_to_analyze) {
      if (!quiet) message("Training model for variable: ", measure)
      
      # Prepare temporary table for training (features + label + id only) the id will be further remove
      tmp <- feat_table %>% inner_join( 
                   strat_metadata[,c(sample_id_column, measure)] %>% 
                     rename(Label = 2),
                   by = sample_id_column) 

      # Cross validation loop
      pred_dataframe_list <- list()
      for (i in 1:nrow(folds)) {
        fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
        
        # Get train/test samples of this fold
        train_samples <- folds$splits[[i]]$in_id
        train_data <- tmp[train_samples,] %>% filter(!is.na(Label))
        test_data <- tmp[-train_samples,] %>% filter(!is.na(Label))

        # Export sample ids from the train and the test data, to new dataframes.
        train_samples_ids <- train_data %>% select(any_of(c(sample_id_column)))
        test_samples_ids <- test_data %>% select(any_of(c(sample_id_column)))
        train_data <- train_data %>% select(-any_of(c(sample_id_column))) # Prepare table for training (features + label only) 
        test_data <- test_data %>% select(-any_of(c(sample_id_column))) # Prepare table for testing (features + label only) 
        
        # Train logistic regression / RF
        if (ml_method == 'glm') {
          lasso_model <- cv.glmnet(
            x = train_data %>% select(-Label) %>% as.matrix(), 
            y = train_data$Label, 
            grouped = (n_samples > 50)) # grouped = FALSE relevant for very small sample sets only
        } else if (ml_method == 'rf') {
          rf_model <- ranger(Label ~ ., data = train_data)
        } else if (ml_method == 'lgbm') {
          dtrain <- lgb.Dataset(train_data %>% select(-Label) %>% as.matrix(),
                                label = train_data$Label)
          dtest <- lgb.Dataset.create.valid(dtrain, 
                                            test_data %>% select(-Label) %>% as.matrix(),
                                            label = test_data$Label)
          params <- list(
            objective = "regression",
            metric = "l2",
            min_data = 1L,
            learning_rate = .3)
          invisible(capture.output(lgbm_model <- lgb.train(
            params = params,
            data = dtrain,
            nrounds = 100L,
            #valids = list(test = dtest),
            verbose = 0L)))
        } else if (ml_method == 'xgboost') {
          watchlist <- list(train = xgb.DMatrix(data = train_data %>% select(-Label) %>% as.matrix(), 
                                                label = train_data$Label),
                            valid = xgb.DMatrix(data = test_data %>% select(-Label) %>% as.matrix(), 
                                                label = test_data$Label))
          numberOfClasses <- length(unique(train_data$Label))
          invisible(capture.output(bst_model <- xgb.train(data = xgb.DMatrix(data = train_data %>% select(-Label) %>% as.matrix(), 
                                                                             label = train_data$Label),
                                                          objective = "reg:squarederror",
                                                          #watchlist = watchlist, 
                                                          nrounds = 50)))
        }
        
        # Make prediction on held out samples
        if (ml_method == 'glm') {
          preds <- predict(lasso_model, 
                           test_data %>% select(-Label) %>% as.matrix(), 
                           s = "lambda.min") %>% unname()
        } else if (ml_method == 'rf') {
          preds <- predict(rf_model, test_data %>% select(-Label))$predictions
        } else if (ml_method == 'lgbm') {
          preds <- predict(lgbm_model, test_data %>% select(-Label) %>% as.matrix()) 
        } else if (ml_method == 'xgboost') {
          test_matrix <- xgb.DMatrix(data = test_data %>% select(-Label) %>% as.matrix(), label = test_data$Label)
          preds <- predict(bst_model, newdata = test_matrix)
        }
        
        test_samples_ids$preds <- preds
        test_samples_ids$strat <- strat_val
        pred_dataframe_list[[i]] <- test_samples_ids
        
        # Calculate model performance for this fold
        tmp_eval <- data.frame(preds = preds, true_label = test_data$Label)  

        cv_correlations <- bind_rows(
          cv_correlations,
          data.frame(
            group = strat_val,
            measure = measure,
            fold_id = fold_id,
            test_cor = cor(tmp_eval$preds, tmp_eval$true_label, method = corr_method),
            test_cor_p = cor.test(tmp_eval$preds, 
                                  tmp_eval$true_label, 
                                  alternative = 'greater', # Check positive association
                                  method = corr_method)$p.value
          )
        )
      }
      final_preds <- do.call(rbind, pred_dataframe_list) %>%
        group_by(sample_id, strat) %>%
        summarize(across(.cols = everything(), .fns = mean))
      final_preds$measure <- measure
      preds_by_measure[[measure]] <- final_preds
    }
    preds_by_strat[[strat_val]] <- do.call(rbind, preds_by_measure) 
  }
  preds_by_strat <- do.call(rbind, preds_by_strat) 
  
  sumlog2 <- function(ps) {ps <- ps[!is.na(ps)]; return(sumlog(ps)$p)}
  
  # Summarise over folds
  summary_correlations <- cv_correlations %>%
    group_by(group, measure) %>%
    summarise(n_folds = n(),
              mean_cor = mean(test_cor, na.rm = TRUE),
              sd_cor = sd(test_cor, na.rm = TRUE),
              p = sumlog2(test_cor_p),
              .groups = "drop") %>%
    mutate(fdr = p.adjust(p)) %>%
    #select(-p) %>%
    mutate(Significant = fdr < fdr_threshold & mean_cor > corr_threshold)
  
   if(stratify_by_col == "stratify_by_col")
     summary_correlations <- summary_correlations %>% select(-group)

  if (! any(summary_correlations$Significant)) {
    message('No significantly well-predicted metadata measures')
  } else {
    message('The following measures were well-predicted by the microbiome using a ', toupper(ml_method), ' model:')
    if(stratify_by_col == "stratify_by_col"){
      message(paste(summary_correlations$measure[summary_correlations$Significant], collapse = ', '))
    } else {
      message(paste0(summary_correlations$measure[summary_correlations$Significant], " (",
                    summary_correlations$group[summary_correlations$Significant], ")", collapse = ', '))
      
    }
  }
  
  return(list(results = summary_correlations, preds = preds_by_strat))
}



#' Creates a PCoA plot, where samples are colored by some continuous variable.
#'    The function also fits a linear model on the PC's with the variable of 
#'    interest as the dependent variable, to assess whether overall microbiome 
#'    composition associates with the variable. A p value for the linear model 
#'    is returned as well as the plot.
#'    
#' @param dist_matrix A microbiome distance matrix (column names correspond 
#'    to sample id's in the first column and are in the same order).
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param var_to_color_by The continuous feature by which to color the PCOA 
#'    points (name of column in metadata).
#' @param samples_to_include List of samples to include. Null if all samples 
#'    should be used.
#' @param dist_metric_name Name of distance metric (for plot title). 
#' @param axis_x Which PC should be used for the x axis? Default: 1. 
#' @param axis_y Which PC should be used for the y axis? Default: 2.
#' 
#' @return A list with 3 elements: 
#'    `pcoa_plot` containing the PCoA plot, 
#'    `lm_r_squared` containing the R^2 value for the linear model
#'    `lm_p_value` containing the p value for the linear model
#' @export
plot_colored_pcoa <- function(dist_matrix,
                              metadata,
                              sample_id_column,
                              var_to_color_by,
                              samples_to_include = NULL,
                              dist_metric_name,
                              axis_x = 1,
                              axis_y = 2,
                              num_pcs_for_lm = 10) {
  require(dplyr)
  require(ggplot2)
  require(ape)
  
  # Validation
  if(! var_to_color_by %in% names(metadata))
    stop(var_to_color_by, " column is missing from metadata")
  
  # Utility function
  get_lm_p_value <- function(lm_model) {
    if (class(lm_model) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(lm_model)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  dist_matrix <- dist_matrix %>%
    tibble::column_to_rownames(var = sample_id_column)
  
  if (!is.null(samples_to_include)) {
    if(all(samples_to_include %in% colnames(dist_matrix))) {
      dist_matrix <- dist_matrix[samples_to_include, samples_to_include]
    } else {
      stop('The samples_to_include argument contains samples not in distance matrix')
    }
  }
  
  # Compute PCoA using the ape package
  dist_matrix <- as.dist(dist_matrix)
  PCOA <- pcoa(dist_matrix)
  pcoa_values <- PCOA$vectors[,1:10] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(sample_id_column)
  
  # Format titles for the axes
  ax_x_title <- paste0('Axis ', axis_x, ' (', round(PCOA$values$Relative_eig[axis_x]*100,1) , '% var. explained)')
  ax_y_title <- paste0('Axis ', axis_y, ' (', round(PCOA$values$Relative_eig[axis_y]*100,1) , '% var. explained)')
  
  # Add metadata to pcoa values, for plotting
  pcoa_values <- pcoa_values %>%
    left_join(metadata, by = sample_id_column)
  
  # Is the continuous variable associated with the main PCs? Fit a linear model
  lm_formula <- formula(paste(var_to_color_by, '~', 
                              paste0('Axis.', 
                                     1:num_pcs_for_lm, 
                                     collapse = ' + ')))
  lm_model <- lm(lm_formula, data = pcoa_values)
  lm_r_squared <- round(summary(lm_model)$adj.r.squared, 2)
  lm_p_value <- round(get_lm_p_value(lm_model),3)
  
  # Plot
  p <- ggplot(pcoa_values,
              aes_string(x = paste0('Axis.', axis_x),
                         y = paste0('Axis.', axis_y),
                         fill = var_to_color_by)) +
    geom_point(color = 'black', alpha = 0.8, shape = 21, size = 4) +
    scale_fill_gradient(low = 'lightcyan', high = 'midnightblue',
                        guide = guide_colorbar(barwidth = 0.6, frame.colour = "black", ticks.colour = "black")) +
    ggtitle(paste0(dist_metric_name, ' ~ ', var_to_color_by)) +
    labs(subtitle = paste0('(Linear model R^2: ',lm_r_squared,', P value: ',lm_p_value,')')) +
    xlab(ax_x_title) +
    ylab(ax_y_title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(legend.title = element_text(size = 8))
  
  return(list(pcoa_plot = p, lm_r_squared = lm_r_squared, lm_p_value = lm_p_value))
}

#' Creates a PCoA plot, where samples are colored by some continuous variable.
#'    The function also fits a linear model on the PC's with the variable of 
#'    interest as the dependent variable, to assess whether overall microbiome 
#'    composition associates with the variable. A p value for the linear model 
#'    is returned as well as the plot. FDR is performed for each distance 
#'    metric separately.
#'    
#' @param dist_matrices A named list of microbiome distance matrices (column 
#'    names correspond to sample id's in the first column and are in the same 
#'    order). Names correspond to distance metric name.
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `feat_table` and `metadata`
#' @param vars_to_color_by All continuous features to test (columns in metadata).
#' @param samples_to_include List of samples to include. Null if all samples 
#'    should be used. 
#' @param fdr_threshold FDR threshold used to define significant findings. Plots
#'    will only be generated for significant findings.
#'  
#' @return A list with two elements:
#'    `results` contains p values for the linear models. 
#'    `plots` is a list of PCoA plots for significant findings. 
#' @export
plot_multiple_colored_pcoa <- function(dist_matrices,
                                       metadata,
                                       sample_id_column,
                                       vars_to_color_by,
                                       samples_to_include = NULL,
                                       fdr_threshold = 0.1) {
  require(dplyr)
  require(rlist)
  
  # Validations
  if(!all(vars_to_color_by %in% names(metadata)))
    stop("Invalid 'vars_to_color_by' argument. Some of the variables provided are missing from the metadata table")
  
  # Place holders
  all_stats <- data.frame()
  final_plots <- list()
  
  # Iterate over distance metrics
  for (dist_metric_name in names(dist_matrices)) {
    # Get distance matrix
    dist_matrix <- dist_matrices[[dist_metric_name]]
    
    # Place holders
    all_plots <- list()
    all_ps <- c()
    all_r2s <- c()
    
    # Iterate over measurements
    for (measure in vars_to_color_by) {
      res <- plot_colored_pcoa(dist_matrix,
                               metadata = metadata, 
                               sample_id_column = sample_id_column,
                               var_to_color_by = measure,
                               samples_to_include = samples_to_include,
                               dist_metric_name = dist_metric_name,
                               axis_x = 1,
                               axis_y = 2)
      all_plots[[measure]] <- res$pcoa_plot
      all_ps[measure] <- res$lm_p_value
      all_r2s[measure] <- res$lm_r_squared
    }
    
    # Run FDR
    all_fdrs <- p.adjust(all_ps, 'fdr')
    
    # Save stats
    all_stats <- bind_rows(
      all_stats,
      data.frame(distance_metric = dist_metric_name,
                 variable = vars_to_color_by,
                 p = all_ps %>% unname(),
                 r_squared = all_r2s %>% unname(),
                 fdr = all_fdrs %>% unname())
    )
    
    # Save plots for tests with FDR < fdr_threshold
    if(all(all_fdrs > fdr_threshold)) {
      message('Using ', dist_metric_name, ' - all FDRs are > ', fdr_threshold)
      measures_to_plot <- c()
    } else {
      message('Using ', dist_metric_name, '- the following tests resulted in FDR < ', fdr_threshold, ':')
      
      for (measure in vars_to_color_by) {
        if (all_fdrs[measure] > fdr_threshold) next;
        message(measure)
        
        # Add FDR values to plot subtitle
        p <- all_plots[[measure]]
        curr_subtitle <- p$labels$subtitle
        new_subtitle <- gsub('\\)$', 
                             paste0(', FDR: ', all_fdrs[measure] %>% round(3), ')'), 
                             curr_subtitle)
        p <- p + labs(subtitle = new_subtitle)
        
        # Save plot to final list
        final_plots <- list.append(final_plots, p)
      }
    }
  }
  
  return(list(results = all_stats, plots = final_plots))
}


#' Calculate mantel test between two distance matrices: 
#' 1: Given microbiome distance matrix (such as Bray-Curtis, Weighted UniFrac)
#' 2: The function calculate distance matrix of the metadata metadata_cols_to_analyze.
#'
#' @param dist_matrix A microbiome distance matrix (column names correspond 
#'    to sample id's in the first column and are in the same order).#' 
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param metadata_cols_to_analyze A vector of column names from the metadata 
#'    table for which the analysis should be carried out. If NULL, all available 
#'    columns will be analyzed. 
#' @param sample_id_column The name of the sample_id column.
#' @param metadata_dist_method The method that will be used to calculate 
#'    distances based on metadata features. Methods supported: "euclidean", 
#'    "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param scale TRUE if scaling the metadata[, metadata_cols_to_analyze] before 
#'    calculating the distance is required (default = TRUE).
#' @param center TRUE if centering the metadata[, metadata_cols_to_analyze] before 
#'    calculating the distance is required (default = TRUE).
#'
#' @return mantel p-value.
#'
get_mantel <- function(dist_matrix, 
                       metadata, 
                       metadata_cols_to_analyze, 
                       sample_id_column = "sample_id",
                       metadata_dist_method = "euclidean",
                       scale = TRUE,
                       center = TRUE) {
  
  require(vegan)
  require(dplyr)
  
  # Data verifications 
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in feat_tale")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
    stop("One of metadata_cols_to_analyze not found in metadata")
  if (nrow(metadata) != nrow(dist_matrix))
    warning("Note that nrow(metadata) != nrow(dist_matrix). Only samples with distance and metadata will be used for the following analysis.\n")

  if (is.null(metadata_cols_to_analyze)) {
    metadata_cols_to_analyze <- colnames(metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.")
  } else {
    metadata <- metadata %>% select(all_of(sample_id_column, metadata_cols_to_analyze))
  }
  
  # Check for missing values
  na_counts_per_row <- rowSums(is.na(metadata))
  if (sum(na_counts_per_row) > 0) {
    metadata <- metadata[na_counts_per_row == 0,]
    warning("Removed ", sum(na_counts_per_row>0), " rows with incomplete metadata")
  }
  
  # Merge the two data frames (keep only samples that exists in both original data frames).
  df <- merge(x = dist_matrix, y = metadata, by = sample_id_column, all.x=FALSE)
  
  # Scale and center the metadata parameters.
  metadata_param_scaled_centered <- scale(df[, metadata_cols_to_analyze], center=center, scale=scale)
  
  # Run mantel and return the mantel p-value.
  mantel_res <- mantel(as.dist(df[df[,sample_id_column]]), 
                       dist(metadata_param_scaled_centered, method = metadata_dist_method), 
                       permutations = 999)
  return(mantel_res$signif)
}

#' Get mixed effect model: metadata_col ~ feature + fixed + (1 | random) for all:
#'  metadata_col in `metadata_cols_to_analyze` (from the metadata),
#'  feature in `feat_table` (except `sample_id_column`),
#'  fixed columns in the fixed_cols (from the metadata),
#'  random columns in the random_cols (from the metadata).
#'
#' @param feat_table A feature table with rows as samples and columns as features. 
#'    A 'sample id' column is expected in order to merge the table with the metadata.
##' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected.
#' @param metadata_cols_to_analyze A vector of columns names from the metadata 
#'    table for which the analysis should be carried out. If NULL, all available 
#'    columns will be analyzed.
#' @param sample_id_column The name of the sample_id column.
#' @param fixed_cols columns in the metadata for the fixed effect in the mixed effect model.
#' @param random_cols columns in the metadata for the random effect in the mixed effect model.
#' @param fdr_threshold FDR threshold used to define significant findings. Plots
#'    will only be generated for significant findings.
#' @param color_by_var column in the metadata to color the output scatter plots 
#'    by. By default resorts to the metadata columns being analyzed.
#' @param z_score_normalize Z-normalize all numric columns in feat_table and metadata
#'    
#' @return list of plots and results:
#'    plots- scatter plots of significant results, feature as x and metadata col in y.
#'    results- table of p-values and FDR for all the models.
#'
get_mixed_effect <- function(feat_table, 
                             metadata, 
                             metadata_cols_to_analyze = NULL, 
                             sample_id_column,
                             fixed_cols = NULL,
                             random_cols = NULL,
                             fdr_threshold = 0.1,
                             log_transform_y = FALSE,
                             color_by_var = NULL,
                             plot_output = TRUE,
                             z_score_normalize = FALSE) {
  require(dplyr)
  require(tidyverse)
  require(lme4)
  require(ggpubr)
  require(ggplot2)
  require(ggpmisc)
  
  # Data verifications 
  if (!sample_id_column %in% colnames(feat_table))
    stop("sample_id_column not found in feat_tale")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  if (!is.null(metadata_cols_to_analyze)){
    if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  }
  if (nrow(metadata) != nrow(feat_table))
    warning("Note that nrow(metadata) != nrow(feat_table). Only samples with distance and metadata will be used for the following analysis.\n")
  if (!is.null(color_by_var))
    if (!color_by_var %in% colnames(metadata))
      stop("color_by_var not found in metadata")
  if (!is.null(fixed_cols)) {
    if (any(!fixed_cols %in% colnames(metadata)))
      stop("Some fixed_cols were not found in metadata")
  }
  if (!is.null(random_cols)) {
    if (any(!random_cols %in% colnames(metadata)))
      stop("Some random_cols were not found in metadata")
  }
  
  if (is.null(metadata_cols_to_analyze)){
    metadata_cols_to_analyze <- colnames(metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    
    # remove fixed_cols and random_cols from metadata_cols_to_analyze
    if (length(fixed_cols) > 0 || length(random_cols) > 0) {
      exclude_cols <- unique(c(fixed_cols, random_cols))
      metadata_cols_to_analyze <- setdiff(metadata_cols_to_analyze, exclude_cols)
    }
    
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.")
  }
  
  
  # Convert vector of fixed effect and random effects to string.
  mixed_string = ""
  if (!is.null(fixed_cols)){
    for (fixed_col in fixed_cols)
      mixed_string =  paste(mixed_string, "+", fixed_col)
  }
  if (!is.null(random_cols)){
    for (random_col in random_cols)
      mixed_string =  paste(mixed_string,"+ ( 1 | ", random_col, ")")
  }
  
  # Merge the feat_table and the metadata data frames.
  df <- merge(x = feat_table, y = metadata, by = sample_id_column, all.x=FALSE)
  
  # Conduct z-score normalization if needed
  if (z_score_normalize) {
    df <- df %>%
      mutate(across(where(is.numeric), ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))
  }
  
  # Get the mixed effect model 
  results <- data.frame()
  for (metadata_col in metadata_cols_to_analyze){
    for (feature in colnames(feat_table[, -which(names(feat_table) == sample_id_column)])){
      
      # message("Working on", metadata_col, "<>", feature)
      
      # Organize table - remove rows with missing values
      df_tmp <- df %>% drop_na(metadata_col)
      df_tmp$metadata_col <- df_tmp[,metadata_col]
      if (log_transform_y) df_tmp$metadata_col <- log(df_tmp$metadata_col)
      df_tmp$feature <- df_tmp[,feature]
      
      # Create final LM formula
      formula_string <- paste("metadata_col", "~", "feature", mixed_string)
      if (is.null(random_cols)) {
        # lmm1 <- lm(formula = as.formula(formula_string), data = df_tmp)
        lmm1 <- lm(formula = as.formula(formula_string), data = df_tmp)
        coefs <- summary(lmm1)$coefficients
        pval <- ifelse('feature' %in% rownames(coefs), coefs['feature','Pr(>|t|)'], NA)
        estimate = coefs%>%as_tibble(rownames = "feature")%>%filter(feature == "feature")%>%pull(Estimate)
      } else {
        lmm1 <- lmer(formula = as.formula(formula_string), data = df_tmp)
        for_p_values <- parameters::p_value(lmm1)
        coefs <- summary(lmm1)$coefficients
        pval <- for_p_values$p[for_p_values$Parameter == "feature"] # ifelse('feature' %in% rownames(coefs), coefs['feature','Pr(>|z|)'], NA)
        estimate = coefs%>%as_tibble(rownames = "feature")%>%filter(feature == "feature")%>%pull(Estimate)
      }
      
      results <- bind_rows(
        results,
        data.frame(metadata_col = metadata_col,
                   feature = feature,
                   pval = pval,
                   N_samples = nrow(df_tmp),
                   estimate = estimate,
                   `formula` = formula_string)
      )
    }
  }
  
  results <- results %>%
    group_by(metadata_col) %>%
    mutate(FDR = p.adjust(pval, 'fdr')) %>%
    ungroup()
  
  # Create plots for significant results.
  if(plot_output == TRUE){
    if (is.null(color_by_var)){ # New
      color_by_var <- "tmp_col"
      df$tmp_col <- 0
    }
    
    plots <- list()
    for (metadata_col in metadata_cols_to_analyze){
      for (feature in colnames(feat_table[, -which(names(feat_table) == sample_id_column)])){
        
        # Plot settings
        # color_by_var_ <- ifelse(is.null(color_by_var), metadata_col, color_by_var # OLD 
        feature_ <- feature
        if (grepl(";g__", feature, fixed=TRUE))
          feature_ <- paste("g__", strsplit(feature_, ";g__")[[1]][2], sep="")
        
        fdr <- results$FDR[results$feature == feature & results$metadata_col == metadata_col]
        
        if (!is.na(fdr) & fdr <= fdr_threshold) {
          p <- ggplot(df, aes(x = .data[[metadata_col]], y = .data[[feature]])) + 
            # TODO
            # geom_path(aes(group = .data[[color_by_var_]]),
            #           color = 'darkgrey',
            #           arrow = arrow(length = unit(4, "points"), type = "closed")) +
            geom_point(aes(fill = .data[[color_by_var]]), color = 'black', shape = 21, size = 3, alpha = 0.7) +  # was: color_by_var_
            theme_classic() +
            ylab(feature_) +
            theme(legend.position = "none") + 
            geom_smooth(method = 'lm', formula = y ~ x, se=FALSE, color="darkred") +
            labs(subtitle = paste('FDR:', round(fdr, 4))) 
          plots[[paste(feature, "_", metadata_col)]] <- p
        }
      }
    }
    
    return(list(plots = plots, all_results = results, sig_results = results %>% filter(FDR <= fdr_threshold) %>% arrange(FDR)))
  }
  return(list(all_results = results, sig_results = results %>% filter(FDR <= fdr_threshold) %>% arrange(FDR)))
}


#' get_corr_heatmap
#'
#' The function plots heatmap from statistically significant correlation scores
#' @param  corr_results A results table from get_correlations_to_metadata function
#' @param xlab x axis name, corresponds to the feat_table
#' @param ylab y axis name, corresponds to the metaadata
#' @param fdr_cuttoff fdr cutoff to filter results with higher q values
#' @param rhos_cutoff filtering cutoff to filter results with low rho values
#' @param fill_palette color plate from Rcolorbrewer (for scale_fill_distiller)  
#' @return plot A heatmap correlation plot

get_corr_heatmap <-  function(corr_results, xlab, ylab, fdr_cutoff = 0.1, rhos_cutoff = 0.2, 
                              low_palette = "#4575b4", high_palette = "#d73027", mid_palette = "#ffffbf"){
  
  data_for_plot <- corr_results%>%
    as_tibble(rownames = "rowname")%>%
    separate(rowname, sep = "_vs._", into = c("feature", "metadata"))%>%
    mutate(metadata = str_remove(metadata, "_strat_*."),
           metadata = str_remove(metadata, "_strat"))%>%
    filter(FDR<fdr_cutoff)%>%
    filter(rhos>rhos_cutoff)
  
  
  plot <- ggplot(data_for_plot, aes( x= feature, y = metadata, fill = rhos))+
    geom_tile()+
    theme_classic(base_size = 15)+
    scale_fill_gradient2(low = low_palette, mid = mid_palette, high = high_palette)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = xlab, y = ylab)

  
  return(plot)

}


#' Run PERMANOVA (significance test for group-level differences) with adonis2 
#' (vegan library). 
#' 
#'    
#' @param metadata A table with the metadata features to analyze. A 
#'    'sample id' column is expected.
#' @param dist_matrix Microbiome distance matrix (column names correspond to 
#'    sample id's in the first column and are in the same order).
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `metadata`
#' @param metadata_cols_to_analyze The continuous columns in the metadata 
#'    to run permanova.
#' @param fdr_threshold FDR threshold used to define significant findings. PCoA 
#'    plots will only be generated for significant findings.
#' @param permutations number of permutations in permanova (default 9999)
#' @low_color lower color in the PCoA plot gradient (default "lightcyan")
#' @high_color higher color in the PCoA plot gradient (default "midnightblue")
#' @return A list of two elements:
#'    `results` Table of p-values and FDRs for all `sample_id_column`.
#'    `plots` is a list of relevant plots: PCoA colored by feature for all
#'    significant results.

get_permanova_continuous <- function(metadata, 
                          dist_matrix, 
                          dist_metric_name,
                          sample_id_column = 'sample_id', 
                          metadata_cols_to_analyze = NULL, 
                          metadata_cols_fix_effects = NULL,
                          metadata_cols_rand_effects = NULL,
                          fdr_threshold = 0.1,
                          permutations = 9999,
                          low_color = "lightcyan",
                          high_color = "midnightblue",
                          midpoint = 0) {
  
  # Required libraries loading
  require(dplyr)
  require(tibble)
  require(vegan)
  
  # Define metadata_cols_to_analyze as all columns if not provided, and check if exist if provided.
  if (!is.null(metadata_cols_to_analyze)){
    if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
      stop("One of metadata_cols_to_analyze not found in metadata")
  } else {
    metadata_cols_to_analyze <- colnames(metadata) 
    metadata_cols_to_analyze <- metadata_cols_to_analyze[!metadata_cols_to_analyze %in% c(sample_id_column)]
    warning("metadata_cols_to_analyze not provided, all metadata columns will be taken except sample_id_column.\n")
  }
  
  # Data verification
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  if (!all(metadata_cols_to_analyze %in% colnames(metadata)))
    stop("One of metadata_cols_to_analyze not found in metadata")
  if (nrow(metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the metadata. Only samples with both features and metadata will be used for the following analysis.\n")
  
  # Merge the dataframes
  df <- merge(x = dist_matrix, y = metadata, by = sample_id_column, all.x=FALSE)  %>%
    remove_rownames %>% 
    column_to_rownames(var=sample_id_column)
  
  # Add fixed effects to a string that later will be added to the formula and check that they are not one of the metadata_cols_to_analyze.
  mixed_string = ""
  if (!is.null(metadata_cols_fix_effects)){
    for (fixed_col in metadata_cols_fix_effects){
      # Validate that all metadata_cols_fix_effects not in metadata_cols_to_analyze.
      if (fixed_col %in% metadata_cols_to_analyze)
        stop(paste0("Shared column between metadata_cols_fix_effects and metadata_cols_to_analyze: ", 
                    fixed_col, "\nTry again with different metadata_cols_fix_effects!\n"))
      mixed_string =  paste0(mixed_string, " + tmp_df[,'", fixed_col, "'] ")
    }
  }
  
  # Run adonis2
  pvals <- c()
  r2_col <- c()
  r2_total <- c()
  for (col in metadata_cols_to_analyze) {
    message('Working on variable: ', col)
    tmp_df <- df[!is.na(df[,col]), ]
    
    formula_string <- paste("tmp_df[, rownames(tmp_df)] ~ tmp_df[[col]]", mixed_string)
    permanova <- adonis2(as.formula(formula_string),
                         permutations = permutations)
    
    pvals[col] <- permanova$`Pr(>F)`[1]
    r2_col[col] <- permanova['tmp_df[[col]]','R2']
    r2_total[col] <- permanova['Total','R2']-permanova['Residual','R2']
  }
  
  # Call plot_colored_pcoa_discrete for significant features.
  all_fdrs <- p.adjust(pvals, 'fdr')
  plots <- list()
  if(!all(all_fdrs > fdr_threshold)) {
    for (metadata_feature in names(all_fdrs[all_fdrs <= fdr_threshold])){
      p12 <- plot_colored_pcoa_continuous(dist_matrix, metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 2, low_color, high_color, midpoint)
      plots[[paste0(metadata_feature,"_pc1&pc2")]] <- p12
      p13 <- plot_colored_pcoa_continuous(dist_matrix, metadata, sample_id_column, metadata_feature, dist_metric_name, 1, 3, low_color, high_color, midpoint)
      plots[[paste0(metadata_feature,"_pc1&pc3")]] <- p13
      p23 <- plot_colored_pcoa_continuous(dist_matrix, metadata, sample_id_column, metadata_feature, dist_metric_name, 2, 3, low_color, high_color, midpoint)
      plots[[paste0(metadata_feature,"_pc2&pc3")]] <- p23
    }
  }
  
  # Return a list of a data frame of the p-values and the FDRs, and the PCoA plots.
  results = data.frame(
      pvals = pvals, 
      FDR = all_fdrs,
      R2_feature_only = r2_col,
      R2_permanova_total = r2_total
    ) %>%
    tibble::rownames_to_column("feature")
  
  return(list(results = results, plots = plots))
}



#' Creates a PCoA plot, where samples are colored by some  variable.
#'   Take only samples that exist in the distance matrix and in the metadata.
#'   
#' @param dist_matrix Microbiome distance matrix (column 
#'    names correspond to sample id's in the first column and are in the same 
#'    order).
#' @param metadata A table with the metadata features to analyze. A 'sample id' 
#'    column is expected. 
#' @param sample_id_column The name of the sample_id column. Present in both 
#'    `dist_matrix` and `metadata`
#' @param vars_to_color_by continuous feature to test (column in metadata).
#' @param dist_metric_name Name of distance metric (for plot title). 
#' @param axis_x Which PC should be used for the x axis? Default: 1. 
#' @param axis_y Which PC should be used for the y axis? Default: 2.
#' @param plot_labels PCoA plot labels in case the discrete metadata are numeric values.
#'
#' @return PCoA plot

plot_colored_pcoa_continuous <- function(dist_matrix,
                                         metadata,
                                       sample_id_column,
                                       var_to_color_by,
                                       dist_metric_name,
                                       axis_x = 1,
                                       axis_y = 2,
                                       low_color = "lightcyan",
                                       high_color = "midnightblue", 
                                       midpoint = 0){
  
  # Required libraries loading
  require(dplyr)
  require(ggplot2)
  require(ape)
  
  # Validation
  if(! var_to_color_by %in% names(metadata))
    stop(var_to_color_by, " column is missing from metadata")
  if (!sample_id_column %in% colnames(dist_matrix))
    stop("sample_id_column not found in dist_matrix")
  if (!sample_id_column %in% colnames(metadata))
    stop("sample_id_column not found in metadata")
  if (nrow(metadata) != nrow(dist_matrix))
    warning("# of samples is different between the distance matrix and the metadata Only samples with features and metadata will be used for the following analysis.\n")
  
  df <- merge(x = dist_matrix, y = metadata, by = sample_id_column, all.x=FALSE)
  df <- df %>% drop_na(all_of(var_to_color_by)) # TODO: instead of dropping missing values, give them a default color?
  
  # Compute PCoA using the ape package
  dist_matrix <- as.dist(df[df[,sample_id_column]])
  PCOA <- pcoa(dist_matrix)
  pcoa_values <- PCOA$vectors[,1:10] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(sample_id_column)
  
  # Format titles for the axes
  ax_x_title <- paste0('Axis ', axis_x, ' (', round(PCOA$values$Relative_eig[axis_x]*100,1) , '% var. explained)')
  ax_y_title <- paste0('Axis ', axis_y, ' (', round(PCOA$values$Relative_eig[axis_y]*100,1) , '% var. explained)')
  
  # Add metadata to pcoa values, for plotting
  pcoa_values <- pcoa_values %>%
    left_join(metadata, by = sample_id_column)

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
    theme(legend.title = element_text(size = 8)) +
#    scale_fill_gradient(colours=c("#CA0020", "#F4A582", "#92C5DE", "#0571B0"),
    scale_fill_gradient2(low = low_color, high = high_color, mid = "white", midpoint = midpoint,
                        guide = guide_colorbar(barwidth = 0.6, frame.colour = "black", ticks.colour = "black"))
  return(p)
}




run_continuous_maaslin <- function(feat_table, 
                                       metadata, 
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
                                       quiet = FALSE) {
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
    metadata_cols_to_analyze = colnames(metadata)
  
  # Check if normalization method fits the data (if only one column: normalization method cannot be 'TSS')
  if ((ncol(feat_table) ==  2) & (normalization_method == "TSS")){
    stop("TSS normlaization is not allowed in case there is only one feature in feat_table. 
         Try again with normalization_method = 'NONE'.")
  }
  
  # Create output directory (will be deleted if not required).
  if (!dir.exists(output_dir))
    dir.create(output_dir)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  metadata <- metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(metadata), row.names(feat_table))
  metadata <- metadata[samples_to_include,]
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
    if (((is.numeric(metadata[, fixed_col])) & 
         (n_distinct(metadata[, fixed_col]) < 5)) | 
        (!is.numeric(metadata[, fixed_col]))) {
      
      # Print a message if there are more than 5 categories in the metadata_cols_fix_effects because it might affect the FDR correction.
      if (n_distinct(metadata[[fixed_col]]) > 5)
        message("Too many unique values in ", fixed_col, " column. It might affect the FDR correction and omit some significant results!")
      
      # Perform one-hot encoding for columns if they have more then 2 categories.
      if (n_distinct(metadata[[fixed_col]]) > 2){
        message("Performing one-hot encoding on ", fixed_col)
        fixed_cols_final_names <- c(fixed_cols_final_names, paste(fixed_col, unique(metadata[[fixed_col]]), sep = '_'))
        tmp_metadata_rownames <- rownames(metadata)
        metadata <- metadata %>% mutate(value = 1)  %>% spread(fixed_col, value,  fill = 0, sep = '_')
        rownames(metadata) <- tmp_metadata_rownames
        
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
    
      # Keep only the two selected values (also drop missing values)
      tmp_metadata <- metadata[(!is.na(metadata[[col]])),]# & (metadata[[col]] %in% pair), ]
      
      # Get the final list of samples for which DA will be calculated
      samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
      tmp_metadata <- tmp_metadata[samples_to_include,,drop=F]
      tmp_feat_table <- feat_table[samples_to_include,,drop=F]
      
      # Remove constant columns if exist
      tmp_metadata <- tmp_metadata[, apply(tmp_metadata, 2, n_distinct) > 1]
      
      # Fix column names if needed
      colnames(tmp_metadata) <- make.names(colnames(tmp_metadata))
      if (!is.null(metadata_cols_rand_effects))
        make.names(metadata_cols_rand_effects)
      
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
      
      #tmp_results$results$pair_a <- as.character(pair[[1]])
      #tmp_results$results$pair_b <- as.character(pair[[2]])
      
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
#    }
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
          
          tmp_metadata <- metadata %>% drop_na(col)
          samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table2))
          tmp_metadata <- tmp_metadata[samples_to_include,,drop=FALSE]
          feat_table2 <- feat_table[samples_to_include,,drop=FALSE]
          
          tmp <- data.frame(
            col_vals = tmp_metadata[[col]],
            feat_vals = feat_table2[[feat]]
          )
          
          
          FDR <- unique_signif_results[1,"FDR"][[1]]
          # Convert p-value to scientific style if needed 
          if (FDR < fdr_threshold / 1000){
            plot_FDR <- format(FDR, scientific = TRUE, digits=4)
          } else {
            plot_FDR <- round(FDR, 4)
          }
          
          p <- ggscatter(tmp, x = "col_vals", y = "feat_vals", title="MaAsLin2", subtitle = paste("FDR: ", plot_FDR),
                         color = "black", shape = 19, size = 3, alpha = 0.3, 
                         add = "reg.line",
                         add.params = list(color = "navy", fill = "lightgray"),
                         conf.int = TRUE,
                         cor.coeff.args = list(method = "spearman",  label.sep = "\n"),
                         xlab = col , ylab = feat
          ) + theme(plot.title = element_text(hjust = 0.5)) + font("xlab", size = 10) + font("ylab", size = 10) +
            font("xy.text", size = 10, color = "gray", face = "bold") 
          if (logscale_y_for_plot)
            p <- p + yscale("log10", .format = TRUE)
          p <- p + theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))
          
          # p <- ggplot(tmp, aes(x = col_vals, group = col_vals, y = feat_vals)) +
          #   geom_point(color = 'black', alpha = 0.8, shape = 21, size = 4) +
          #   geom_smooth(method = 'lm', formula = feat_vals ~ col_vals, se=FALSE, color="darkred") +
          #   theme_classic() +
          #   xlab(col) +
          #   ylab(feat) +
          #   labs(subtitle = 'MaAsLin2') +
          #   theme(plot.subtitle = element_text(hjust = 0.5)) 
          # 
          # # Log scale y axis if needed
          # if (logscale_y_for_plot) {
          #   p <- p + yscale("log10", .format = TRUE)
          # }
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