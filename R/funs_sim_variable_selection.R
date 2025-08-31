# ==========================
# Data generation and simulation functions
# ==========================
source("R/funs_models.R")

#' Generate synthetic data for IV regression
#' @param n Sample size
#' @param p Number of endogenous variables
#' @param q Number of instruments
#' @param seed Random seed
#' @return List containing Y, X, Z, true beta
generate_data <- function(n = 500, p = 10, q = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  beta_true <- c(1.5, -0.5, 0.8, rep(0, p-3))
  gamma_true <- matrix(runif(q*p, -0.5, 0.5), q, p)
  
  Z <- matrix(rnorm(n*q), n, q)
  
  A <- matrix(runif((p+1)^2, -1, 1), p+1, p+1)
  Sigma <- A %*% t(A) + 0.01 * diag(p+1)
  
  errors <- mvtnorm::rmvnorm(n, mean = rep(0, p+1), sigma = Sigma)
  
  X <- Z %*% gamma_true + errors[, 1:p]
  Y <- X %*% beta_true + errors[, p+1]
  
  list(Y = Y, X = X, Z = Z, beta_true = beta_true)
}

#' Calculate metrics for one replication
#' @param results List of outputs from bayesian_iv for different methods
#' @param beta_true True beta vector
#' @return List containing bias, MSE, coverage, variable selection metrics
calculate_metrics <- function(results, beta_true) {
  p <- length(beta_true)
  n_methods <- length(results)
  
  bias <- mse <- coverage <- matrix(0, n_methods, p)
  
  for (m in 1:n_methods) {
    if (!is.null(results[[m]])) {
      for (j in 1:p) {
        bias[m,j] <- results[[m]]$posterior_means$beta[j] - beta_true[j]
        mse[m,j] <- (results[[m]]$posterior_means$beta[j] - beta_true[j])^2
        ci <- results[[m]]$credible_intervals$beta[j, ]
        coverage[m,j] <- as.numeric(beta_true[j] >= ci[1] & beta_true[j] <= ci[2])
      }
    } else {
      bias[m,] <- mse[m,] <- coverage[m,] <- NA
    }
  }
  
  # Variable selection metrics
  zero_idx <- which(beta_true == 0)
  nonzero_idx <- which(beta_true != 0)
  vs_metrics <- matrix(NA, n_methods, 6)
  colnames(vs_metrics) <- c("TP", "FP", "FN", "FPR", "FNR", "Precision")
  
  for (m in 1:n_methods) {
    if (!is.null(results[[m]])) {
      intervals <- results[[m]]$credible_intervals$beta
      contains_zero <- sapply(1:p, function(j) intervals[j,1] <= 0 & intervals[j,2] >= 0)
      
      TP <- sum(contains_zero == FALSE & beta_true != 0)
      FP <- sum(contains_zero == FALSE & beta_true == 0)
      FN <- sum(contains_zero == TRUE & beta_true != 0)
      
      FPR <- if (length(zero_idx) > 0) FP / length(zero_idx) else 0
      FNR <- if (length(nonzero_idx) > 0) FN / length(nonzero_idx) else 0
      Precision <- if (TP + FP > 0) TP / (TP + FP) else NA
      
      vs_metrics[m,] <- c(TP, FP, FN, FPR, FNR, Precision)
    }
  }
  
  list(bias = bias, mse = mse, coverage = coverage, vs_metrics = vs_metrics)
}

#' Run a single replication
#' @param rep_id Replication index (used for seed)
#' @param n Sample size
#' @param p Number of endogenous variables
#' @param q Number of instruments
#' @return List containing metrics, true beta, success flags
run_single_replication <- function(rep_id, n = 500, p = 10, q = 10) {
  data <- generate_data(n, p, q, seed = rep_id*100)
  
  methods <- c("M1", "M2", "M3", "M3_Lasso")
  results <- list()
  
  for (method in methods) {
    results[[method]] <- tryCatch(
      bayesian_iv(data$Y, data$X, data$Z, method = method, niter = 2000, burn = 500),
      error = function(e) NULL
    )
  }
  
  metrics <- calculate_metrics(results, data$beta_true)
  
  list(metrics = metrics, beta_true = data$beta_true, success = !sapply(results, is.null))
}

#' Run full simulation
#' @param n_reps Number of replications
#' @param n Sample size
#' @param p Number of endogenous variables
#' @param q Number of instruments
#' @return List with summary statistics and parameters
run_simulation <- function(n_reps = 100, n = 500, p = 10, q = 10) {
  methods <- c("M1", "M2", "M3", "M3_Lasso")
  
  # Storage arrays
  all_bias <- array(NA, dim = c(n_reps, length(methods), p))
  all_mse <- array(NA, dim = c(n_reps, length(methods), p))
  all_coverage <- array(NA, dim = c(n_reps, length(methods), p))
  all_vs <- array(NA, dim = c(n_reps, length(methods), 6))
  success_count <- rep(0, length(methods))
  names(success_count) <- methods
  
  for (rep in 1:n_reps) {
    cat("Replication", rep, "/", n_reps, "\n")
    res <- run_single_replication(rep, n, p, q)
    for (i in 1:length(methods)) {
      if (res$success[i]) {
        all_bias[rep,i,] <- res$metrics$bias[i,]
        all_mse[rep,i,] <- res$metrics$mse[i,]
        all_coverage[rep,i,] <- res$metrics$coverage[i,]
        all_vs[rep,i,] <- res$metrics$vs_metrics[i,]
        success_count[i] <- success_count[i] + 1
      }
    }
  }
  
  summary_stats <- list()
  for (i in 1:length(methods)) {
    method <- methods[i]
    if (success_count[i] > 0) {
      valid_idx <- which(!is.na(all_bias[,i,1]))
      summary_stats[[method]] <- list(
        n_success = success_count[i],
        mean_bias = colMeans(all_bias[valid_idx,i,,drop=FALSE], na.rm = TRUE),
        mean_mse = colMeans(all_mse[valid_idx,i,,drop=FALSE], na.rm = TRUE),
        mean_coverage = colMeans(all_coverage[valid_idx,i,,drop=FALSE], na.rm = TRUE),
        mean_vs = colMeans(all_vs[valid_idx,i,,drop=FALSE], na.rm = TRUE)
      )
    }
  }
  
  list(summary = summary_stats, parameters = list(n_reps = n_reps, n = n, p = p, q = q, beta_true = res$beta_true))
}

#' Print simulation results
#' @param results Output of run_simulation
print_results <- function(results) {
  cat("\n=== SIMULATION RESULTS ===\n")
  params <- results$parameters
  cat("Replications:", params$n_reps, "\n")
  cat("Sample size:", params$n, "\n")
  cat("Parameters: p =", params$p, ", q =", params$q, "\n")
  cat("True beta:", paste(round(params$beta_true,3), collapse=", "), "\n\n")
  
  cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s\n", 
              "Method", "Success", "MSE", "Coverage", "FPR", "FNR", "Precision"))
  cat(paste(rep("-",70), collapse=""), "\n")
  
  for (method in names(results$summary)) {
    stats <- results$summary[[method]]
    overall_mse <- sqrt(sum(stats$mean_mse))
    avg_coverage <- mean(stats$mean_coverage)
    cat(sprintf("%-10s %8d %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                method, stats$n_success, overall_mse, avg_coverage,
                stats$mean_vs[4], stats$mean_vs[5], stats$mean_vs[6]))
  }
}
