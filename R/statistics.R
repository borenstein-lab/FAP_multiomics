library(glmnet)
library(rsample)
library(pROC)
library(dplyr)
library(permuco)
library(emmeans)

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

create_dummy_vars <- function(df, target_vars){
  # Step 1: Remove factors and find character features
  df <- remove_factors(df)
  char_cols <- sapply(df %>% select(all_of(target_vars)), function(x) is.character(x))
  nonchar_df <- df[!char_cols]
  cat_df <- lapply(df[char_cols], as.factor) %>% as.data.frame()
  
  # Step 2: Create dummy variables (drop intercept so all dummies are separate)
  char_dummies <- model.matrix(~ . - 1, data = cat_df)
  drop_dummies <- lapply(cat_df, function(f) levels(f)[1])
  
  # Step 3: Remove the first level of each factor (the intercept)
  drop_names <- unlist(mapply(function(name, val) paste0(name, val), names(drop_dummies), drop_dummies))
  char_dummies <- char_dummies[, !colnames(char_dummies) %in% drop_names, drop = FALSE]
  
  # Step 4: Combine back with non-categorical columns
  df_binary <- cbind(nonchar_df, char_dummies)
  colnames(df_binary) <- make.names(colnames(df_binary))
  return(df_binary)
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