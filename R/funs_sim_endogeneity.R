# ===============================
# R/funs_sim_endogeneity.R
# ===============================
source("R/funs_models.R")


# Load required packages
library(parallel)
library(doParallel)
library(foreach)
library(mvtnorm)
library(MASS)
library(MCMCpack)

# -------------------------------
# 1. Single/multivariate simulation data generation
# -------------------------------

sim_gen_multivariate <- function(n, p, q, beta_true, gamma_true, Sigma_u_true, sigma_eps2_true, sigma_u_eps) {
  Z <- matrix(rnorm(n * q), nrow = n, ncol = q)
  Sigma <- rbind(cbind(Sigma_u_true, sigma_u_eps), cbind(t(sigma_u_eps), sigma_eps2_true))
  u_eps <- mvtnorm::rmvnorm(n, sigma = Sigma)
  u_errors <- u_eps[, 1:p, drop = F]
  epsilon <- u_eps[, p+1, drop = F]
  X <- Z %*% gamma_true + u_errors
  Y <- X %*% beta_true + epsilon
  return(list(Y = Y, X = X, Z = Z, true_values = list(beta = beta_true, gamma = gamma_true, Sigma = Sigma)))
}

# -------------------------------
# 2. Run Bayesian IV simulation with multiple configs
# -------------------------------

simulate_bayesian_iv_multiple_configs <- function(n_sim = 100, 
                                                  sample_sizes = c(100, 500),
                                                  configs = list(
                                                    "p<q" = list(p=2, q=3),
                                                    "p=q" = list(p=3, q=3),
                                                    "p>q" = list(p=5, q=3)
                                                  ),
                                                  n_cores = parallel::detectCores() - 1) {
  
  # Source the model functions (bayesian_iv)
  source("R/funs_bayesian_iv_model.R")  # Adjust path if needed
  
  all_results <- list()
  
  for (config_name in names(configs)) {
    cat("\n\n========== Configuration:", config_name, "==========\n")
    p <- configs[[config_name]]$p
    q <- configs[[config_name]]$q
    config_results <- list()
    
    for (n in sample_sizes) {
      cat("\nStarting simulations for n =", n, ", p =", p, ", q =", q, "\n")
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      clusterEvalQ(cl, library(mvtnorm))
      clusterExport(cl, c("sim_gen_multivariate", "bayesian_iv"), envir = environment())
      
      parallel_results <- parLapply(cl, 1:n_sim, function(sim_idx) {
        set.seed(sim_idx)
        beta_true <- rep(1, p)
        gamma_true <- matrix(runif(p*q), q, p)
        Sigma_u_true <- 0.5^abs(outer(1:p, 1:p, "-"))
        sigma_eps2_true <- 1
        sigma_u_eps <- matrix(rep(-0.4, p), ncol = 1)
        
        data <- sim_gen_multivariate(n, p, q, beta_true, gamma_true, Sigma_u_true, sigma_eps2_true, sigma_u_eps)
        y <- data$Y
        X <- data$X
        Z <- data$Z
        
        tryCatch({
          res_m1 <- bayesian_iv(y, X, Z, method="M1", niter=5000, burn=1000)
          res_m2 <- bayesian_iv(y, X, Z, method="M2", niter=5000, burn=1000)
          res_m3 <- bayesian_iv(y, X, Z, method="M3", niter=5000, burn=1000)
          list(res_m1=res_m1, res_m2=res_m2, res_m3=res_m3, error=FALSE, beta_true=beta_true)
        }, error=function(e) list(error=TRUE, message=e$message, sim_idx=sim_idx))
      })
      
      stopCluster(cl)
      config_results[[paste0("n", n)]] <- parallel_results
      cat("Completed simulations for n =", n, ", p =", p, ", q =", q, "\n")
    }
    
    all_results[[config_name]] <- config_results
  }
  
  return(all_results)
}

# -------------------------------
# 3. Result printing function
# -------------------------------

print_simulation_results <- function(all_results) {
  for (config_name in names(all_results)) {
    cat("\n\n========== RESULTS FOR", config_name, "==========\n")
    config_results <- all_results[[config_name]]
    
    first_n_key <- names(config_results)[1]
    n_sim <- length(config_results[[first_n_key]])
    p <- length(config_results[[first_n_key]][[1]]$beta_true)
    
    for (n_key in names(config_results)) {
      cat("\n=== Sample Size:", substr(n_key, 2, nchar(n_key)), "===\n")
      # For brevity, you can implement bias/MSE/coverage calculation here
      # Example: just print number of successful simulations
      n_successful <- sum(sapply(config_results[[n_key]], function(x) !isTRUE(x$error)))
      cat("Number of successful simulations:", n_successful, "out of", n_sim, "\n")
    }
  }
}
