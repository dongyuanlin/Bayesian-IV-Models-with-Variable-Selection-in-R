# ==========================
# 02_sim_variable_selection.R
# Script to run simulation study for variable selection
# ==========================

# 1️⃣ Load model and simulation functions
source("../R/funs_models.R")                # Bayes IV model functions
source("../R/funs_sim_variable_selection.R") # Simulation helpers

# 2️⃣ Quick test function (optional)
# run_test <- function(n_reps = 5) {
#   results <- run_simulation(n_reps = n_reps, n = 100, p = 5, q = 5)
#   print_results(results)
#   return(results)
# }

# 3️⃣ Full simulation
run_main_simulation <- function() {
  results <- run_simulation(n_reps = 100, n = 500, p = 10, q = 10)
  print_results(results)
  return(results)
}

# 4️⃣ Usage example (commented out if you don't want automatic run)
# main_results <- run_main_simulation()
# saveRDS(main_results, "lassosimulation.rds")

# 5️⃣ Read and print saved results (optional)
# main_results <- readRDS("lassosimulation.rds")
# print_results(main_results)

cat("Available functions:\n")
cat("- run_test(): Quick test (5 reps)\n")
cat("- run_main_simulation(): Full simulation (100 reps)\n")
cat("- run_simulation(n_reps, n, p, q): Custom simulation\n")
