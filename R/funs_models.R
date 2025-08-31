# ==========================
# Bayesian IV Models (with LASSO extension)
# ==========================

# Required packages
library(mvtnorm)
library(MASS)
library(MCMCpack)
library(brms)

#' Bayesian IV Regression
#'
#' Implements several Bayesian IV methods:
#' M1: one-stage regression ignoring instruments
#' M2: full Bayesian using Y,X,Z
#' M3: modified gamma posterior
#' M3_Lasso: Bayesian LASSO variant
#'
#' @param y Response vector
#' @param X Endogenous covariates matrix
#' @param Z Instruments matrix
#' @param niter Number of MCMC iterations
#' @param burn Number of burn-in iterations
#' @param method One of "M1","M2","M3","M3_Lasso"
#' @param lambda Regularization parameter for LASSO
#' @param r Hyperparameter for lambda prior
#' @param delta Hyperparameter for lambda prior
#' @return A list containing posterior samples, means, sds, and credible intervals
#' @export
bayesian_iv <- function(y, X, Z, niter = 5000, burn = 1000, 
                        method = c("M1", "M2", "M3", "M3_Lasso"),
                        lambda = 1, r = 1, delta = 1) {
  
  # Match method argument
  method <- match.arg(method)
  
  # Ensure input data is in matrix/vector format
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (!is.vector(y)) y <- as.vector(y)
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  
  # Center variables
  y_centered <- as.vector(scale(y, center = TRUE, scale = FALSE))
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Z_centered <- scale(Z, center = TRUE, scale = FALSE)
  
  # Prior hyperparameters
  mu_beta <- rep(0, p)
  sigma2_beta <- 10
  mu_gamma <- matrix(0, q, p)
  sigma2_gamma <- 10
  nu0 <- p + 5
  Psi0 <- diag(p + 1)
  epsilon <- 1e-8
  
  # Storage
  beta_samples <- matrix(1, niter, p)
  gamma_samples <- array(0, dim = c(niter, q, p))
  Sigma_samples <- array(0, dim = c(niter, p + 1, p + 1))
  
  if (method == "M3_Lasso") {
    tau2_samples <- matrix(1, niter, p)
    lambda2_samples <- rep(lambda^2, niter)
  }
  
  # -----------------------
  # Implement different methods
  # -----------------------
  
  if (method == "M2") {
    # ===============================================================
    # M2:FULL BAYESIAN IMPLEMENTATION (gamma | X, Y, Z)
    # ===============================================================
    
    # Initialize parameters
    beta_current <- rep(0, p)
    gamma_current <- matrix(0, q, p)
    Sigma_current <- diag(p + 1)
    
    for (i in 1:niter) {
      # 1. Sample beta
      # According to equations (33)-(35) in the PDF
      
      # Calculate conditional variance σ²_{ε|u}
      Sigma_u <- Sigma_current[1:p, 1:p, drop = FALSE]
      Sigma_ue <- Sigma_current[1:p, p+1, drop = FALSE]
      sigma2_e <- Sigma_current[p+1, p+1]
      sigma2_e_u <- sigma2_e - t(Sigma_ue) %*% solve(Sigma_u) %*% Sigma_ue
      
      # Calculate Wi = Yi - Σ'_uε Σ^{-1}_u (Xi - Zi Γ)^T
      W = y_centered - (X_centered - Z_centered %*% gamma_current) %*% solve(Sigma_u, Sigma_ue)
      
      # Calculate posterior variance for beta
      V_beta <- ginv(t(X_centered) %*% X_centered / as.numeric(sigma2_e_u) + 
                       diag(1/sigma2_beta, p))
      
      # Calculate posterior mean for beta
      m_beta <- V_beta %*% (t(X_centered) %*% W / as.numeric(sigma2_e_u) + 
                              mu_beta / sigma2_beta)
      
      # Sample beta from posterior distribution
      beta_current <- mvtnorm::rmvnorm(1, mean = m_beta, sigma = V_beta)[1,]
      
      # 2. Sample gamma
      # Based on equations (36)-(44) in the PDF
      
      # Calculate matrix B1 (equation 41)
      B1 <- t(Z_centered) %*% X_centered %*% solve(Sigma_u) #q \times p
      
      # Calculate C = Σ'_uε Σ^{-1}_u
      C <- t(Sigma_ue) %*% solve(Sigma_u) # 1 \times p
      
      # Calculate residuals
      res_Y <- y_centered - X_centered %*% beta_current
      
      # Calculate Q_i = Z'_i Z_i (equation 39)
      Q <- t(Z_centered) %*% Z_centered
      
      # Calculate R = C'C (equation 40)
      R <- t(C) %*% C # p \times p
      
      # Loop over each equation
      B2 <- 1/(as.numeric(sigma2_e_u)) * t(Z_centered) %*%(W-X_centered %*% t(C)) %*% C
      
      B3 <- mu_gamma / sigma2_gamma
      
      # Calculate Λ (equation 37)
      Lambda <- kronecker(solve(Sigma_u), Q) + kronecker(R,  Q) / as.numeric(sigma2_e_u) +
        1/sigma2_gamma * kronecker(diag(p), diag(q))
      
      vecB = c(B1 + B2 + B3)
      mutildeLambda = solve(Lambda, vecB)
      gamma_current_vec <- mvtnorm::rmvnorm(1, mutildeLambda, solve(Lambda))
      gamma_current <- matrix(gamma_current_vec, q, p) #gamma_true
      
      # 3. Sample Σ
      # According to equation (53) in the PDF
      
      # Calculate residual matrix
      res_X <- X_centered - Z_centered %*% gamma_current
      res_Y <- y_centered - X_centered %*% beta_current
      
      R <- cbind(res_X, res_Y)
      S <- crossprod(R)
      
      # Sample Σ from inverse-Wishart posterior
      Sigma_current <- MCMCpack::riwish(nu0 + n, Psi0 + S)
      
      # Store current samples
      beta_samples[i,] <- beta_current
      gamma_samples[i,,] <- gamma_current
      Sigma_samples[i,,] <- Sigma_current
    }
    
  } else if (method == "M3") {
    # ===============================================================
    # M3:change gamma posterior: gamma | X, Z
    # ===============================================================
    
    # Initialize parameters
    beta_current <- rep(0, p)
    gamma_current <- matrix(0, q, p)
    Sigma_current <- diag(p + 1)
    
    for (i in 1:niter) {
      # 1. Sample beta
      
      # Calculate conditional variance σ²_{ε|u}
      Sigma_u <- Sigma_current[1:p, 1:p, drop = FALSE]
      Sigma_ue <- Sigma_current[1:p, p+1, drop = FALSE]
      sigma2_e <- Sigma_current[p+1, p+1]
      sigma2_e_u <- sigma2_e - t(Sigma_ue) %*% solve(Sigma_u) %*% Sigma_ue
      
      # Calculate Wi = Yi - Σ'_uε Σ^{-1}_u (Xi - Zi Γ)^T
      W = y_centered - (X_centered - Z_centered %*% gamma_current) %*% solve(Sigma_u, Sigma_ue)
      
      # Calculate posterior variance for beta
      V_beta <- ginv(t(X_centered) %*% X_centered / as.numeric(sigma2_e_u) + 
                       diag(1/sigma2_beta, p))
      
      # Calculate posterior mean for beta
      m_beta <- V_beta %*% (t(X_centered) %*% W / as.numeric(sigma2_e_u) + 
                              mu_beta / sigma2_beta)
      
      # Sample beta from posterior distribution
      beta_current <- mvtnorm::rmvnorm(1, mean = m_beta, sigma = V_beta)[1,]
      
      # 2. Sample gamma
      #just sample gamma by x z Sigma:
      Sigma_u <- Sigma_current[1:p, 1:p, drop = FALSE]
      
      V_gamma <- ginv(t(Z_centered) %*% Z_centered + diag(1/sigma2_gamma, q))
      
      mu_gamma_posterior <- V_gamma %*% (t(Z_centered) %*% X_centered +
                                           mu_gamma/sigma2_gamma)
      
      vec_gamma_mean <- as.vector(mu_gamma_posterior)
      vec_gamma_var <- kronecker(Sigma_u, V_gamma)
      
      gamma_vec <- mvtnorm::rmvnorm(1, mean = vec_gamma_mean, sigma = vec_gamma_var)[1,]
      gamma_current <- matrix(gamma_vec, nrow = q, ncol = p)
      
      
      # 3. Sample Σ
      # According to equation (53) in the PDF
      
      # Calculate residual matrix
      res_X <- X_centered - Z_centered %*% gamma_current
      res_Y <- y_centered - X_centered %*% beta_current
      
      R <- cbind(res_X, res_Y)
      S <- crossprod(R)
      
      # Sample Σ from inverse-Wishart posterior
      Sigma_current <- MCMCpack::riwish(nu0 + n, Psi0 + S)
      
      # Store current samples
      beta_samples[i,] <- beta_current
      gamma_samples[i,,] <- gamma_current
      Sigma_samples[i,,] <- Sigma_current
    }
    
  } else if (method == "M3_Lasso") {
    # ===============================================================
    # M3_Lasso: Bayesian LASSO Implementation based on M3
    # ===============================================================
    
    # Initialize parameters
    beta_current <- rep(0, p)
    gamma_current <- matrix(0, q, p)
    Sigma_current <- diag(p + 1)
    tau2_current <- rep(1, p)       # Initialize τ² parameters for LASSO
    lambda2_current <- lambda^2     # Initialize λ² parameter
    
    for (i in 1:niter) {
      # 1. Sample beta with LASSO prior
      D_tau <- diag(tau2_current, p)
      
      Sigma_u <- Sigma_current[1:p, 1:p, drop = FALSE]
      Sigma_ue <- Sigma_current[1:p, p+1, drop = FALSE]
      sigma2_e <- Sigma_current[p+1, p+1]
      sigma2_e_u <- sigma2_e - t(Sigma_ue) %*% solve(Sigma_u) %*% Sigma_ue
      
      W = y_centered - (X_centered - Z_centered %*% gamma_current) %*% solve(Sigma_u, Sigma_ue)
      
      # Calculate posterior variance for beta
      V_beta <- ginv(t(X_centered) %*% X_centered / as.numeric(sigma2_e_u) + 
                       solve(D_tau))
      
      # Calculate posterior mean for beta
      m_beta <- V_beta %*% (t(X_centered) %*% W / as.numeric(sigma2_e_u))
      
      # Sample beta from posterior distribution
      beta_current <- mvtnorm::rmvnorm(1, mean = m_beta, sigma = V_beta)[1,]
      
      # 2. Sample tau^2 parameters (LASSO-specific)
      for (j in 1:p) {
        # Parameters for inverse Gaussian
        mu_param <- sqrt(lambda2_current/(beta_current[j]^2 + epsilon))
        shape_param <- lambda2_current
        
        # Sample from inverse Gaussian distribution
        tau2_current[j] <- 1/brms::rinv_gaussian(1, mu = mu_param, shape = shape_param)
      }
      
      # 3. Sample lambda^2 (optional - if using random λ)
      lambda2_current <- rgamma(1, shape = p + r, 
                                rate = 0.5 * sum(tau2_current) + delta)
      
      # 4. Sample gamma (same as M3)
      Sigma_u <- Sigma_current[1:p, 1:p, drop = FALSE]
      
      V_gamma <- ginv(t(Z_centered) %*% Z_centered + diag(1/sigma2_gamma, q))
      
      mu_gamma_posterior <- V_gamma %*% (t(Z_centered) %*% X_centered +
                                           mu_gamma/sigma2_gamma)
      
      vec_gamma_mean <- as.vector(mu_gamma_posterior)
      vec_gamma_var <- kronecker(Sigma_u, V_gamma)
      
      gamma_vec <- mvtnorm::rmvnorm(1, mean = vec_gamma_mean, sigma = vec_gamma_var)[1,]
      gamma_current <- matrix(gamma_vec, nrow = q, ncol = p)
      
      # 5. Sample Sigma
      # Calculate residual matrix
      res_X <- X_centered - Z_centered %*% gamma_current
      res_Y <- y_centered - X_centered %*% beta_current
      
      R <- cbind(res_X, res_Y)
      S <- crossprod(R)
      
      # Sample Σ from inverse-Wishart posterior
      Sigma_current <- MCMCpack::riwish(nu0 + n, Psi0 + S)
      
      # Store current samples
      beta_samples[i,] <- beta_current
      gamma_samples[i,,] <- gamma_current
      Sigma_samples[i,,] <- Sigma_current
      tau2_samples[i,] <- tau2_current
      lambda2_samples[i] <- lambda2_current
    }
    
  }  else if (method == "M1") {
    # ===============================================================
    # M1: ONE-STAGE BAYESIAN REGRESSION (ignoring instruments)
    # ===============================================================
    
    # Initialize parameters
    beta_current <- rep(0, p)
    sigma_epsilon2_current <- 1
    
    
    for (i in 1:niter) {
      # Sample beta
      V_beta <- solve(t(X_centered) %*% X_centered / sigma_epsilon2_current + 
                        diag(1/sigma2_beta, p))
      m_beta <- V_beta %*% (t(X_centered) %*% y_centered / sigma_epsilon2_current + 
                              mu_beta / sigma2_beta)
      
      # Sample beta from posterior
      beta_current <- mvtnorm::rmvnorm(1, mean = m_beta, sigma = V_beta)[1,]
      
      # Sample sigma_epsilon2
      residual_y <- y_centered - X_centered %*% beta_current
      shape <- (n/2) + 2
      rate <- 2 + sum(residual_y^2)/2
      
      # Sample from inverse-gamma distribution
      sigma_epsilon2_current <- 1/rgamma(1, shape = shape, rate = rate)
      
      # Store current samples
      beta_samples[i,] <- beta_current
      Sigma_samples[i,p+1,p+1] <- sigma_epsilon2_current
      
      # Set other parameters to NA (not used in one-stage)
      gamma_samples[i,,] <- NA
      Sigma_samples[i,1:p,1:p] <- NA
      Sigma_samples[i,1:p,p+1] <- NA
      Sigma_samples[i,p+1,1:p] <- NA
    }
  }
  
  # -----------------------
  # Process posterior samples (after burn-in)
  # -----------------------
  posterior_samples <- list(
    beta = beta_samples[(burn+1):niter,, drop = FALSE],
    gamma = gamma_samples[(burn+1):niter,,, drop = FALSE],
    Sigma = Sigma_samples[(burn+1):niter,,, drop = FALSE]
  )
  
  if (method == "M3_Lasso") {
    posterior_samples$tau2 <- tau2_samples[(burn+1):niter,, drop = FALSE]
    posterior_samples$lambda2 <- lambda2_samples[(burn+1):niter]
  }
  
  # Posterior means and sds
  posterior_means <- list(
    beta = colMeans(posterior_samples$beta),
    gamma = apply(posterior_samples$gamma, c(2,3), mean, na.rm = TRUE),
    Sigma = apply(posterior_samples$Sigma, c(2,3), mean, na.rm = TRUE)
  )
  
  if (method == "M3_Lasso") {
    posterior_means$tau2 <- colMeans(posterior_samples$tau2)
    posterior_means$lambda2 <- mean(posterior_samples$lambda2)
  }
  
  posterior_sds <- list(
    beta = apply(posterior_samples$beta, 2, sd),
    gamma = apply(posterior_samples$gamma, c(2,3), sd, na.rm = TRUE),
    Sigma = apply(posterior_samples$Sigma, c(2,3), sd, na.rm = TRUE)
  )
  
  if (method == "M3_Lasso") {
    posterior_sds$tau2 <- apply(posterior_samples$tau2, 2, sd)
    posterior_sds$lambda2 <- sd(posterior_samples$lambda2)
  }
  
  credible_intervals <- list(
    beta = t(apply(posterior_samples$beta, 2, function(x) quantile(x, c(0.025, 0.975))))
  )
  
  # Return results
  return(list(
    posterior_samples = posterior_samples,
    posterior_means = posterior_means,
    posterior_sds = posterior_sds,
    credible_intervals = credible_intervals,
    method = method
  ))
}
