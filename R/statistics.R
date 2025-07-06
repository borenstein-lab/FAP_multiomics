## Wilcox test between pairs of values in a character column in a df 
add_significance <- function(df, pval_col = "p_value"){
  df$fdr_p_value <- p.adjust(df[[pval_col]], method = "fdr")
  df <- add_sigstart(df,"fdr_p_value")
  return(df)
}

add_sigstart <- function(df, fdr_p_value = "fdr_p_value") {
  df$significance <- case_when(
    df[fdr_p_value] < 0.001 ~ "***",
    df[fdr_p_value] < 0.01  ~ "**",
    df[fdr_p_value] < 0.1  ~ "*",
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

test_wilcox_residuals <- function(df, regressors, x_var, y_var, show_plots = F){
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
  wilcox <- wilcox.test(df[["residuals"]] ~ df[[y_var]], exact = FALSE)
  ggp <- ggplot(df, aes_string(y_var, "residuals", color = y_var)) +
    plot_box_and_jitter() + labs(y = "Residuals")
  return(list(wilcox$p.value, ggp))
}