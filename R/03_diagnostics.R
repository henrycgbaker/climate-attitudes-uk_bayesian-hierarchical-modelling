# ==============================================================================
# Climate Attitudes UK - Bayesian Hierarchical Model
# Script 03: Model Diagnostics
#
# Description:
#   Comprehensive diagnostic assessment including:
#   - Prior-predictive checks
#   - Convergence diagnostics (Rhat, ESS, divergences)
#   - Posterior-predictive checks
#   - R² estimation
#
# Inputs:
#   - data/stan_data_full.rds
#   - outputs/model/fit_prior_predictive.rds
#   - outputs/model/fit_full.rds
#
# Outputs:
#   - outputs/diagnostic_plots/*.png
#
# Dependencies:
#   - Run 01_wrangling.R and 02_modelling.R first
#
# Author: Henry Baker
# ==============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggdist)
library(here)

set.seed(42)

# ── 1. Load Data & Model Fits =================================================

stan_data <- readRDS(here("data", "stan_data_full.rds"))
fit_prior <- readRDS(here("outputs", "model", "fit_prior_predictive.rds"))
fit_full  <- readRDS(here("outputs", "model", "fit_full.rds"))

# Create output directory
plot_dir <- here("outputs", "diagnostic_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ── 2. Extract Stan Data ======================================================

y_opt   <- stan_data$y_opt
i_opt   <- stan_data$i_opt
j_opt   <- stan_data$j_opt
N_opt   <- stan_data$N_opt
J_opt   <- stan_data$J_opt

y_env   <- stan_data$y_env
i_env   <- stan_data$i_env
j_env   <- stan_data$j_env
N_env   <- stan_data$N_env
J_env   <- stan_data$J_env

y_rad   <- stan_data$y_rad
i_rad   <- stan_data$i_rad
j_rad   <- stan_data$j_rad
N_rad   <- stan_data$N_rad
J_rad   <- stan_data$J_rad

N <- stan_data$N

# Build long-format data frames
opt_long_df <- tibble(respondent = i_opt, item = j_opt, response = y_opt)
env_long_df <- tibble(respondent = i_env, item = j_env, response = y_env)
rad_long_df <- tibble(respondent = i_rad, item = j_rad, response = y_rad)

# ── 3. Prior-Predictive Checks ================================================

message("Generating prior-predictive checks...")

# Extract prior-predictive simulated values
y_opt_sim_pr     <- fit_prior$draws(variables = "y_opt_sim")
y_opt_sim_mat_pr <- as_draws_matrix(y_opt_sim_pr)

y_env_sim_pr     <- fit_prior$draws(variables = "y_env_sim")
y_env_sim_mat_pr <- as_draws_matrix(y_env_sim_pr)

y_rad_sim_pr     <- fit_prior$draws(variables = "y_rad_sim")
y_rad_sim_mat_pr <- as_draws_matrix(y_rad_sim_pr)

# Observed responses
obs_opt <- opt_long_df$response
obs_env <- env_long_df$response
obs_rad <- rad_long_df$response

n_plot_draws_pr <- 200
draws_to_plot_opt_pr <- sample(nrow(y_opt_sim_mat_pr), size = n_plot_draws_pr)
draws_to_plot_env_pr <- sample(nrow(y_env_sim_mat_pr), size = n_plot_draws_pr)
draws_to_plot_rad_pr <- sample(nrow(y_rad_sim_mat_pr), size = n_plot_draws_pr)

# (1) Density overlay: optimism
sim_opt_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_opt_sim_mat_pr)),
  obs_index = rep(1:ncol(y_opt_sim_mat_pr), times = n_plot_draws_pr),
  y_opt_sim = as.vector(y_opt_sim_mat_pr[draws_to_plot_opt_pr, ])
)

ppc_prior_opt <- ggplot() +
  geom_density(
    data  = sim_opt_pr_df,
    aes(x = y_opt_sim, group = draw),
    color = "skyblue", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_opt),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Optimism (Item Densities)",
    subtitle = "Blue = prior draws; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "prior_predictive_optimism_density.png"), ppc_prior_opt,
       width = 6, height = 4, dpi = 300)

# (2) Density overlay: environment
sim_env_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_env_sim_mat_pr)),
  obs_index = rep(1:ncol(y_env_sim_mat_pr), times = n_plot_draws_pr),
  y_env_sim = as.vector(y_env_sim_mat_pr[draws_to_plot_env_pr, ])
)

ppc_prior_env <- ggplot() +
  geom_density(
    data  = sim_env_pr_df,
    aes(x = y_env_sim, group = draw),
    color = "lightgreen", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_env),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Environment (Item Densities)",
    subtitle = "Green = prior draws; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "prior_predictive_environment_density.png"), ppc_prior_env,
       width = 6, height = 4, dpi = 300)

# (3) Density overlay: radical-reform
sim_rad_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_rad_sim_mat_pr)),
  obs_index = rep(1:ncol(y_rad_sim_mat_pr), times = n_plot_draws_pr),
  y_rad_sim = as.vector(y_rad_sim_mat_pr[draws_to_plot_rad_pr, ])
)

ppc_prior_rad <- ggplot() +
  geom_density(
    data  = sim_rad_pr_df,
    aes(x = y_rad_sim, group = draw),
    color = "salmon", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_rad),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Radical-Reform (Item Densities)",
    subtitle = "Red = prior draws; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "prior_predictive_radical_density.png"), ppc_prior_rad,
       width = 6, height = 4, dpi = 300)

# ── 4. Convergence Diagnostics ================================================

message("Running convergence diagnostics...")

# Summary statistics
summ_full <- fit_full$summary()

bad_rhat     <- summ_full %>% filter(rhat > 1.01) %>% arrange(desc(rhat))
bad_ess_bulk <- summ_full %>% filter(ess_bulk < 200) %>% arrange(ess_bulk)
bad_ess_tail <- summ_full %>% filter(ess_tail < 200) %>% arrange(ess_tail)

message("Parameters with Rhat > 1.01: ", nrow(bad_rhat))
message("Parameters with ESS_bulk < 200: ", nrow(bad_ess_bulk))
message("Parameters with ESS_tail < 200: ", nrow(bad_ess_tail))

# Sampler diagnostics
sampler_diag_df <- as_tibble(fit_full$diagnostic_summary())
message("Divergent transitions: ", sum(sampler_diag_df$num_divergent))
message("Max treedepth hits: ", sum(sampler_diag_df$num_max_treedepth))

# Rhat histogram
rhat_df <- tibble(rhat = summ_full$rhat)

rhat_hist <- ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(
    breaks = seq(0.99, 1.05, by = 0.001),
    fill = "steelblue", color = "white"
  ) +
  geom_vline(xintercept = 1.01, linetype = "dotted", color = "red", linewidth = 1) +
  labs(title = "Histogram of Rhat for All Parameters", x = "Rhat", y = "Count") +
  theme_minimal()

ggsave(here(plot_dir, "rhat_histogram.png"), rhat_hist, width = 8, height = 6, dpi = 150)

# Traceplots for selected parameters
draws_arr <- as_draws_array(fit_full$draws())

pars_to_plot <- c(
  "tau_eta[1]", "tau_eta[2]", "tau_eta[3]",
  "sigma_alpha[1]", "sigma_alpha[2]", "sigma_alpha[3]",
  "lambda_opt[1]", "lambda_env[1]", "lambda_rad[1]"
)

traceplot_sel <- mcmc_trace(
  draws_arr,
  pars = pars_to_plot,
  facet_args = list(ncol = 1, strip.position = "right")
) +
  ggtitle("Traceplots for Selected Parameters") +
  theme_minimal()

ggsave(here(plot_dir, "traceplots_selected_params.png"), traceplot_sel,
       width = 6, height = 10, dpi = 300)

mcmc_dens_plot <- mcmc_dens(
  draws_arr,
  pars = pars_to_plot,
  facet_args = list(ncol = 1, strip.position = "right")
) +
  ggtitle("Posterior Densities for Selected Parameters") +
  theme_minimal()

ggsave(here(plot_dir, "densities_selected_params.png"), mcmc_dens_plot,
       width = 6, height = 10, dpi = 300)

# ── 5. Posterior-Predictive Checks ============================================

message("Generating posterior-predictive checks...")

draws_full_df   <- as_draws_df(fit_full$draws())
y_opt_sim_draws <- draws_full_df %>% select(starts_with("y_opt_sim["))
y_env_sim_draws <- draws_full_df %>% select(starts_with("y_env_sim["))
y_rad_sim_draws <- draws_full_df %>% select(starts_with("y_rad_sim["))

n_plot_draws_pp <- 200
draws_to_plot_opt_pp <- sample(nrow(y_opt_sim_draws), size = n_plot_draws_pp)
draws_to_plot_env_pp <- sample(nrow(y_env_sim_draws), size = n_plot_draws_pp)
draws_to_plot_rad_pp <- sample(nrow(y_rad_sim_draws), size = n_plot_draws_pp)

# (1) PPC Density: optimism
sim_opt_df <- y_opt_sim_draws %>%
  slice(draws_to_plot_opt_pp) %>%
  pivot_longer(cols = everything(), names_to = "obs_index", values_to = "y_opt_sim") %>%
  mutate(draw = rep(1:length(draws_to_plot_opt_pp), each = ncol(y_opt_sim_draws)))

ppc_plot_opt <- ggplot() +
  geom_density(
    data = sim_opt_df,
    aes(x = y_opt_sim, group = draw),
    color = "steelblue", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_opt),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Optimism (Item Densities)",
    subtitle = "Blue = simulated; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "posterior_predictive_optimism_density.png"), ppc_plot_opt,
       width = 6, height = 4, dpi = 300)

# (2) PPC Density: environment
sim_env_df <- y_env_sim_draws %>%
  slice(draws_to_plot_env_pp) %>%
  pivot_longer(cols = everything(), names_to = "obs_index", values_to = "y_env_sim") %>%
  mutate(draw = rep(1:length(draws_to_plot_env_pp), each = ncol(y_env_sim_draws)))

ppc_plot_env <- ggplot() +
  geom_density(
    data = sim_env_df,
    aes(x = y_env_sim, group = draw),
    color = "lightgreen", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_env),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Environment (Item Densities)",
    subtitle = "Green = simulated; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "posterior_predictive_environment_density.png"), ppc_plot_env,
       width = 6, height = 4, dpi = 300)

# (3) PPC Density: radical-reform
sim_rad_df <- y_rad_sim_draws %>%
  slice(draws_to_plot_rad_pp) %>%
  pivot_longer(cols = everything(), names_to = "obs_index", values_to = "y_rad_sim") %>%
  mutate(draw = rep(1:length(draws_to_plot_rad_pp), each = ncol(y_rad_sim_draws)))

ppc_plot_rad <- ggplot() +
  geom_density(
    data = sim_rad_df,
    aes(x = y_rad_sim, group = draw),
    color = "salmon", alpha = 0.1, linewidth = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_rad),
    aes(x = response),
    color = "black", linewidth = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Radical-Reform (Item Densities)",
    subtitle = "Red = simulated; black = observed",
    x = "Standardised response", y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "posterior_predictive_radical_density.png"), ppc_plot_rad,
       width = 6, height = 4, dpi = 300)

# ── 6. R² Posterior Densities =================================================

message("Computing R² estimates...")

r2_draws_df <- fit_full$draws(
  variables = c("R2_opt", "R2_env", "R2_rad"),
  format = "draws_df"
)

r2_summary <- r2_draws_df %>%
  pivot_longer(cols = everything(), names_to = "block", values_to = "R2") %>%
  group_by(block) %>%
  summarize(
    mean_R2   = mean(R2),
    median_R2 = median(R2),
    lower95   = quantile(R2, 0.025),
    upper95   = quantile(R2, 0.975),
    .groups   = "drop"
  )

message("\nR² Summary:")
print(r2_summary)

r2_long <- r2_draws_df %>%
  pivot_longer(cols = starts_with("R2_"), names_to = "block", values_to = "R2")

r2_density_plot <- ggplot(r2_long, aes(x = R2, fill = block)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Posterior Densities of R² (Optimism vs. Environment vs. Radical)",
    x = expression(R^2),
    y = "Density"
  ) +
  theme_minimal()

ggsave(here(plot_dir, "r2_posterior_density.png"), r2_density_plot,
       width = 6, height = 4, dpi = 300)

message("\n=========================================")
message("Diagnostics completed successfully!")
message("All plots saved to: ", plot_dir)
message("=========================================")
