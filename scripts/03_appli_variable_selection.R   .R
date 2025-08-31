# ===============================
# 02_appli_variable_selection.R
# ===============================

# --- 1. Load function library ---
source("../R/funs_appli_variable_selection.R")  # Adjust path according to your project structure

# --- 2. Load data ---
# NOTICE: NEED to find the data by yourself!!!!! if not, all the afetr code can not be run!
dat <- readRDS("dat_IV.Rds")

X <- dat$day1[, !names(dat$day1) %in% c("SEQN", "gender", "race", "bmi")]
Z <- dat$day2[, !names(dat$day2) %in% c("SEQN", "gender", "race", "bmi")]
y <- dat$outcome$bmi

# Optional: Check data
# cat("Dimensions:\n")
# cat("X:", dim(X), "\n")
# cat("Z:", dim(Z), "\n")
# cat("y:", length(y), "\n")
# 
# cat("Missing values:\n")
# cat("y:", sum(is.na(y)), "\n")
# cat("X:", sum(is.na(X)), "\n")
# cat("Z:", sum(is.na(Z)), "\n")

# --- 3. Analyze data ---
results <- analyze_real_data(y, X, Z)

# --- 4. Save results ---
saveRDS(results, file = "lassoapplication.Rds")

# --- 5. Print results table ---
print_results_table(results)

# --- 6. Plotting: comparison of all variables ---
var_names <- paste0("β", 1:ncol(X))  # Generate variable names based on number of columns in X

# Comparison across all methods
p_all <- create_all_variables_comparison(results, var_names)
print(p_all)

# Individual method plots
p_m1_all <- create_single_method_all_vars(results, "M1", var_names)
p_m3_all <- create_single_method_all_vars(results, "M3", var_names)
p_lasso_all <- create_single_method_all_vars(results, "M3_Lasso", var_names)

print(p_m1_all)
print(p_m3_all)
print(p_lasso_all)

# --- 7. Significant variable analysis ---
significant_result <- create_significant_vars_plot(results, var_names)
summary_table <- get_significant_summary(results, var_names)

print(summary_table)
print(significant_result$plot)

# Optional: Save significant variable plot
# ggsave("significant_variables_plot.pdf", plot = significant_result$plot, width = 8, height = 6)
