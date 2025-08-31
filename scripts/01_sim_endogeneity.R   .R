# ===============================
# 01_sim_endogeneity.R
# ===============================

# --- 1. Load functions ---
source("../R/funs_models.R")   
source("../R/funs_sim_endogeneity.R")    # Simulation functions


# --- 2. Set simulation parameters ---
n_sim <- 100
sample_sizes <- c(100, 500)
configs <- list(
  "p<q" = list(p=2, q=3),
  "p=q" = list(p=3, q=3),
  "p>q" = list(p=5, q=3)
)

# --- 3. Run simulations ---
sim_results <- simulate_bayesian_iv_multiple_configs(
  n_sim = n_sim,
  sample_sizes = sample_sizes,
  configs = configs,
  n_cores = parallel::detectCores() - 1
)

# --- 4. Save simulation results ---
saveRDS(sim_results, file = "multivariate_sim_results.Rds")

# --- 5. Print summarized results ---
print_simulation_results(sim_results)
