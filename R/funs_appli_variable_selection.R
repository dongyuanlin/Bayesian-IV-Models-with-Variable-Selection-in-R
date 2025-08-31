# ===============================
# R/funs_appli_variable_selection.R
# ===============================
source("R/funs_models.R")

library(ggplot2)

# -------------------------------
# 1. Data analysis function
# -------------------------------

#' Run Bayesian IV analysis on real data
#' @param y Outcome vector
#' @param X Endogenous variables matrix
#' @param Z Instruments matrix
#' @return List of results for methods M1, M3, M3_Lasso
analyze_real_data <- function(y, X, Z) {
  results <- list()
  results$M1 <- bayesian_iv(y, X, Z, method = "M1", niter = 1000, burn = 200)
  results$M3 <- bayesian_iv(y, X, Z, method = "M3", niter = 1000, burn = 200)
  results$M3_Lasso <- bayesian_iv(y, X, Z, method = "M3_Lasso", niter = 1000, burn = 200)
  return(results)
}

# -------------------------------
# 2. Print results table function
# -------------------------------

#' Print posterior results for all methods
#' @param results Output of analyze_real_data()
#' @param var_names Optional variable names vector
print_results_table <- function(results, var_names = NULL) {
  methods <- names(results)
  if (is.null(var_names)) {
    p <- length(results[[1]]$posterior_means$beta)
    var_names <- paste0("β", 1:p)
  }
  for (method in methods) {
    cat("\n===== Results for", method, "=====\n")
    results_table <- data.frame(
      Parameter = var_names,
      Posterior_Mean = results[[method]]$posterior_means$beta,
      CI_Lower = results[[method]]$credible_intervals$beta[, 1],
      CI_Upper = results[[method]]$credible_intervals$beta[, 2]
    )
    print(results_table)
  }
}

# -------------------------------
# 3. Plotting functions
# -------------------------------

#' Compare all variables across methods
create_all_variables_comparison <- function(results, var_names = NULL) {
  first_method <- names(results)[1]
  n_vars <- length(results[[first_method]]$posterior_means$beta)
  if (is.null(var_names)) var_names <- paste0("β", 1:n_vars)
  
  methods <- names(results)
  
  plot_data <- data.frame(
    Variable = rep(1:n_vars, length(methods)),
    Variable_Name = rep(var_names, length(methods)),
    Method = rep(methods, each = n_vars),
    Estimate = unlist(lapply(results, function(x) x$posterior_means$beta)),
    Lower = unlist(lapply(results, function(x) x$credible_intervals$beta[,1])),
    Upper = unlist(lapply(results, function(x) x$credible_intervals$beta[,2]))
  )
  
  p <- ggplot(plot_data) +
    geom_errorbarh(aes(y = factor(Variable, levels = n_vars:1),
                       xmin = Lower, xmax = Upper,
                       color = Method),
                   height = 0.3,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(y = factor(Variable, levels = n_vars:1),
                   x = Estimate,
                   shape = Method, 
                   color = Method),
               size = 2,
               position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.7) +
    scale_shape_manual(values = c("M1" = 5, "M3" = 1, "M3_Lasso" = 2)) +
    scale_color_manual(values = c("M1" = "red", "M3" = "blue", "M3_Lasso" = "darkgreen")) +
    scale_y_discrete(labels = var_names[n_vars:1]) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    ) +
    labs(x = "Coefficient Estimates",
         y = "Variables",
         title = "All Variables: M1 vs M3 vs M3_Lasso",
         color = "Method",
         shape = "Method")
  
  return(p)
}

#' Plot single method results
create_single_method_all_vars <- function(results, method_name, var_names = NULL) {
  n_vars <- length(results[[method_name]]$posterior_means$beta)
  if (is.null(var_names)) var_names <- paste0("β", 1:n_vars)
  
  estimates <- results[[method_name]]$posterior_means$beta
  ci_lower <- results[[method_name]]$credible_intervals$beta[,1]
  ci_upper <- results[[method_name]]$credible_intervals$beta[,2]
  
  plot_data <- data.frame(
    Variable = 1:n_vars,
    Variable_Name = var_names,
    Estimate = estimates,
    Lower = ci_lower,
    Upper = ci_upper
  )
  
  method_colors <- c("M1" = "red", "M3" = "blue", "M3_Lasso" = "darkgreen")
  method_shapes <- c("M1" = 5, "M3" = 1, "M3_Lasso" = 2)
  method_titles <- c("M1" = "One-Stage (M1)", "M3" = "Bayesian IV (M3)", "M3_Lasso" = "Bayesian IV with Lasso (M3_Lasso)")
  
  color <- method_colors[method_name]
  shape <- method_shapes[method_name]
  title <- method_titles[method_name]
  
  p <- ggplot(plot_data) +
    geom_errorbarh(aes(y = factor(Variable, levels = n_vars:1),
                       xmin = Lower, xmax = Upper),
                   height = 0.3, color = color) +
    geom_point(aes(y = factor(Variable, levels = n_vars:1),
                   x = Estimate),
               shape = shape, size = 2, color = color) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.7) +
    scale_y_discrete(labels = var_names[n_vars:1]) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    labs(x = "Coefficient Estimates",
         y = "Variables",
         title = title)
  
  return(p)
}

# -------------------------------
# 4. Significant variables functions
# -------------------------------

#' Plot significant variables across methods
create_significant_vars_plot <- function(results, var_names = NULL) {
  first_method <- names(results)[1]
  n_vars <- length(results[[first_method]]$posterior_means$beta)
  if (is.null(var_names)) var_names <- paste0("β", 1:n_vars)
  
  methods <- names(results)
  significant_union <- logical(n_vars)
  
  # Identify variables significant in at least one method
  for (method in methods) {
    ci_lower <- results[[method]]$credible_intervals$beta[,1]
    ci_upper <- results[[method]]$credible_intervals$beta[,2]
    significant_union <- significant_union | (ci_lower > 0 | ci_upper < 0)
  }
  
  selected_indices <- which(significant_union)
  full_data <- list()
  
  for (method in methods) {
    ci_lower <- results[[method]]$credible_intervals$beta[,1]
    ci_upper <- results[[method]]$credible_intervals$beta[,2]
    
    full_data[[method]] <- data.frame(
      Variable = selected_indices,
      Variable_Name = var_names[selected_indices],
      Method = method,
      Estimate = results[[method]]$posterior_means$beta[selected_indices],
      Lower = ci_lower[selected_indices],
      Upper = ci_upper[selected_indices]
    )
  }
  
  all_plot_data <- do.call(rbind, full_data)
  
  p <- ggplot(all_plot_data) +
    geom_errorbarh(aes(y = factor(Variable_Name, levels = rev(unique(Variable_Name))),
                       xmin = Lower, xmax = Upper,
                       color = Method),
                   height = 0.3,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(y = factor(Variable_Name, levels = rev(unique(Variable_Name))),
                   x = Estimate,
                   shape = Method, 
                   color = Method),
               size = 3,
               position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.7) +
    scale_shape_manual(values = c("M1" = 5, "M3" = 1, "M3_Lasso" = 2)) +
    scale_color_manual(values = c("M1" = "red", "M3" = "blue", "M3_Lasso" = "darkgreen")) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom") +
    labs(x = "Coefficient Estimates",
         y = "Variables",
         color = "Method",
         shape = "Method")
  
  return(list(plot = p, data = all_plot_data))
}

#' Get summary table of significant variables
get_significant_summary <- function(results, var_names = NULL) {
  first_method <- names(results)[1]
  n_vars <- length(results[[first_method]]$posterior_means$beta)
  if (is.null(var_names)) var_names <- paste0("β", 1:n_vars)
  
  methods <- names(results)
  summary_table <- data.frame()
  
  for (method in methods) {
    ci_lower <- results[[method]]$credible_intervals$beta[,1]
    ci_upper <- results[[method]]$credible_intervals$beta[,2]
    significant_indices <- which(ci_lower > 0 | ci_upper < 0)
    
    for (i in significant_indices) {
      summary_table <- rbind(summary_table, data.frame(
        Method = method,
        Variable = var_names[i],
        Estimate = results[[method]]$posterior_means$beta[i],
        CI_Lower = ci_lower[i],
        CI_Upper = ci_upper[i],
        Significant = ifelse(ci_lower[i] > 0, "Positive", "Negative")
      ))
    }
  }
  
  return(summary_table)
}
