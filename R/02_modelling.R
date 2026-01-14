# ==============================================================================
# Climate Attitudes UK - Bayesian Hierarchical Model
# Script 02: Model Fitting
#
# Description:
#   Compiles and fits the 3-dimensional hierarchical latent trait model using
#   CmdStan. Runs both prior-predictive sampling and full posterior sampling.
#
# Inputs:
#   - data/stan_data_full.rds (from 01_wrangling.R)
#   - stan/hierarchical_latent_traits.stan
#
# Outputs:
#   - outputs/model/fit_prior_predictive.rds
#   - outputs/model/fit_full.rds
#
# Dependencies:
#   - Run 01_wrangling.R first
#   - Requires CmdStan installation
#
# Author: Henry Baker
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(here)

# ── 1. Load Data ==============================================================

stan_data_path <- here("data", "stan_data_full.rds")
if (!file.exists(stan_data_path)) {
  stop("Stan data not found. Run 01_wrangling.R first.")
}

stan_data <- readRDS(stan_data_path)
message("Loaded stan_data with N = ", stan_data$N, " respondents")

# ── 2. Create output directory ================================================

output_dir <- here("outputs", "model")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ── 3. Sanity checks on stan_data =============================================

if (!is.list(stan_data)) {
  stop("`stan_data` must be a list.")
}

required_names <- c(
  "N", "P", "X", "R", "region_id", "Q", "party_id",
  "J_opt", "N_opt", "i_opt", "j_opt", "y_opt", "lower_opt", "upper_opt",
  "J_env", "N_env", "i_env", "j_env", "y_env", "lower_env", "upper_env",
  "J_rad", "N_rad", "i_rad", "j_rad", "y_rad", "lower_rad", "upper_rad"
)
missing_names <- setdiff(required_names, names(stan_data))
if (length(missing_names) > 0) {
  stop(paste0("`stan_data` is missing these elements: ", paste(missing_names, collapse = ", ")))
}

# Validate dimensions
stopifnot(is.numeric(stan_data$N) && length(stan_data$N) == 1 && stan_data$N >= 1)
stopifnot(is.numeric(stan_data$P) && length(stan_data$P) == 1 && stan_data$P >= 1)
stopifnot(is.matrix(stan_data$X))
stopifnot(nrow(stan_data$X) == stan_data$N && ncol(stan_data$X) == stan_data$P)

# Validate region and party indices
stopifnot(length(stan_data$region_id) == stan_data$N)
stopifnot(all(stan_data$region_id >= 1) && all(stan_data$region_id <= stan_data$R))
stopifnot(length(stan_data$party_id) == stan_data$N)
stopifnot(all(stan_data$party_id >= 1) && all(stan_data$party_id <= stan_data$Q))

# Validate opt block
stopifnot(length(stan_data$i_opt) == stan_data$N_opt)
stopifnot(length(stan_data$j_opt) == stan_data$N_opt)
stopifnot(length(stan_data$y_opt) == stan_data$N_opt)
stopifnot(!any(is.na(stan_data$y_opt)))

# Validate env block
stopifnot(length(stan_data$i_env) == stan_data$N_env)
stopifnot(length(stan_data$j_env) == stan_data$N_env)
stopifnot(length(stan_data$y_env) == stan_data$N_env)
stopifnot(!any(is.na(stan_data$y_env)))

# Validate rad block
stopifnot(length(stan_data$i_rad) == stan_data$N_rad)
stopifnot(length(stan_data$j_rad) == stan_data$N_rad)
stopifnot(length(stan_data$y_rad) == stan_data$N_rad)
stopifnot(!any(is.na(stan_data$y_rad)))

message("All sanity checks passed! stan_data is consistent.")

# ── 4. Compile Stan Model =====================================================

stan_model_path <- here("stan", "hierarchical_latent_traits.stan")
if (!file.exists(stan_model_path)) {
  stop("Stan model not found at: ", stan_model_path)
}

message("Compiling Stan model...")
model <- cmdstan_model(stan_model_path)
message("Model compiled successfully.")

# ── 5. Prior-Predictive Sampling ==============================================

message("Running prior-predictive sampling...")
fit_prior <- model$sample(
  data                = stan_data,
  chains              = 4,
  parallel_chains     = 4,
  iter_warmup         = 2000,
  iter_sampling       = 1000,
  refresh             = 100,
  fixed_param         = TRUE,
  save_cmdstan_config = TRUE
)

prior_path <- here("outputs", "model", "fit_prior_predictive.rds")
fit_prior$save_object(prior_path)
message("Prior-predictive fit saved to: ", prior_path)

# ── 6. Full Posterior Sampling ================================================

message("Running full posterior sampling (this may take a while)...")
fit_full <- model$sample(
  data                = stan_data,
  chains              = 4,
  parallel_chains     = 4,
  iter_warmup         = 2000,
  iter_sampling       = 1000,
  refresh             = 100,
  save_cmdstan_config = TRUE,
  init                = 0.5
)

full_path <- here("outputs", "model", "fit_full.rds")
fit_full$save_object(full_path)
message("Full posterior fit saved to: ", full_path)

# ── 7. Quick Diagnostics Summary ==============================================

message("\n=== Quick Diagnostics Summary ===")

# Sampler diagnostics
diag_summary <- fit_full$diagnostic_summary()
message("Divergent transitions: ", sum(diag_summary$num_divergent))
message("Max treedepth hits: ", sum(diag_summary$num_max_treedepth))

# Summary statistics
summ <- fit_full$summary()
bad_rhat <- sum(summ$rhat > 1.01, na.rm = TRUE)
bad_ess <- sum(summ$ess_bulk < 200, na.rm = TRUE)
message("Parameters with Rhat > 1.01: ", bad_rhat)
message("Parameters with ESS_bulk < 200: ", bad_ess)

message("\n=========================================")
message("Modelling complete!")
message("=========================================")
