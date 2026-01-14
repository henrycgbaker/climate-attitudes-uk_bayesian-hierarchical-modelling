# ==============================================================================
# Climate Attitudes UK - Bayesian Hierarchical Model
# Script 04: Quantities of Interest & Visualisations
#
# Description:
#   Extracts and visualises posterior estimates for publication-ready figures:
#   - Covariate effects (slopes)
#   - Regional and party intercepts
#   - Latent correlations
#   - Measurement reliability (McDonald's ω)
#   - Demographic profile predictions
#
# Inputs:
#   - data/stan_data_full.rds
#   - outputs/model/fit_full.rds
#
# Outputs:
#   - outputs/qoi_plots/*.png
#   - outputs/qoi_plots/*.csv
#
# Dependencies:
#   - Run 01_wrangling.R, 02_modelling.R, 03_diagnostics.R first
#
# Author: Henry Baker
# ==============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggdist)
library(forcats)
library(cowplot)
library(corrr)
library(ggraph)
library(igraph)
library(scales)
library(plotly)
library(ggrepel)
library(ggridges)
library(htmlwidgets)
library(here)

set.seed(42)

# ── 1. Load Data & Model Fit ==================================================

stan_data <- readRDS(here("data", "stan_data_full.rds"))
fit_full  <- readRDS(here("outputs", "model", "fit_full.rds"))

# Create output directories
plot_dir <- here("outputs", "qoi_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(here(plot_dir, "profiles_separate"))) {
  dir.create(here(plot_dir, "profiles_separate"), recursive = TRUE)
}

# Define region and party labels
region_levels <- c(
 "Yorkshire & the Humber", "West Midlands", "Scotland", "Wales",
  "North West", "Eastern", "South West", "East Midlands", "London", "South East"
)

party_levels <- c(
  "Green", "Labour", "Plaid Cymru", "Scottish National Party (SNP)",
  "Liberal Democrat", "Conservative", "Reform UK",
  "Another party", "Don't know", "Won't vote"
)

# ── 2. Extract Posterior Draws ================================================

draws_df <- as_draws_df(fit_full$draws())
n_draws  <- nrow(draws_df)

# ── 3. Posterior Slope Coefficients ===========================================

P <- stan_data$P
slope_names <- paste0("B[", rep(1:3, each = P), ",", rep(1:P, times = 3), "]")

slopes_df <- draws_df %>%
  select(all_of(slope_names)) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  arrange(param)

write_csv(slopes_df, here(plot_dir, "slopes_summary.csv"))

# Covariate names in order
cov_names <- c(
  "gender",
  "age_65plus", "age_55_64", "age_45_54", "age_35_44", "age_25_34",
  "edu_L4plus", "edu_L3", "edu_L2", "edu_L1", "edu_appr", "edu_other",
  "insecurity"
)

# Long format for plotting
slopes_draws_long <- draws_df %>%
  select(all_of(slope_names)) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(
    latent = case_when(
      str_detect(param, "^B\\[1,") ~ "φ (Optimism)",
      str_detect(param, "^B\\[2,") ~ "θ (Environment)",
      str_detect(param, "^B\\[3,") ~ "ψ (Radical-Reform)",
      TRUE ~ NA_character_
    ),
    cov_index = as.integer(str_remove_all(param, "B\\[[123],|\\]"))
  ) %>%
  mutate(covariate = factor(cov_names[cov_index], levels = cov_names)) %>%
  select(latent, covariate, cov_index, value)

# Add nudge for overlapping ridgelines
slopes_draws_long <- slopes_draws_long %>%
  mutate(
    y_position = case_when(
      latent == "φ (Optimism)"       ~ cov_index - 0.15,
      latent == "θ (Environment)"    ~ cov_index,
      latent == "ψ (Radical-Reform)" ~ cov_index + 0.15,
      TRUE ~ cov_index
    )
  )

# Overlapping ridgeline plot
covariate_ridge_overlap <- ggplot(
  slopes_draws_long,
  aes(x = value, y = y_position, fill = latent, color = latent)
) +
  geom_density_ridges(
    aes(height = after_stat(density), group = interaction(latent, covariate)),
    stat = "density",
    scale = 3,
    rel_min_height = 0.01,
    alpha = 0.6,
    trim = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_y_continuous(
    breaks = seq_len(length(levels(slopes_draws_long$covariate))),
    labels = levels(slopes_draws_long$covariate)
  ) +
  scale_fill_discrete(name = "Latent Trait") +
  scale_color_discrete(name = "Latent Trait") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Individual-level Covariate Effects on Latent Traits",
    x = "Posterior Slope Value",
    y = "Covariate"
  ) +
  theme(
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_blank()
  )

ggsave(here(plot_dir, "covariate_effects_overlapping_ridgelines.png"),
       covariate_ridge_overlap, width = 9, height = 15, dpi = 300)

# ── 4. Posterior Latent Scores ================================================

N <- stan_data$N
phi_names   <- paste0("phi[", 1:N, "]")
theta_names <- paste0("theta[", 1:N, "]")
psi_names   <- paste0("psi[", 1:N, "]")

phi_draws   <- draws_df %>% select(all_of(phi_names))
theta_draws <- draws_df %>% select(all_of(theta_names))
psi_draws   <- draws_df %>% select(all_of(psi_names))

phi_means   <- colMeans(phi_draws)
theta_means <- colMeans(theta_draws)
psi_means   <- colMeans(psi_draws)

latent_scores_df <- tibble(
  respondent = 1:N,
  phi_hat    = phi_means,
  theta_hat  = theta_means,
  psi_hat    = psi_means
)

write_csv(latent_scores_df, here(plot_dir, "latent_scores.csv"))

# ── 5. Average Marginal Effects ===============================================

ame_df <- slopes_draws_long %>%
  group_by(latent, covariate) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  arrange(latent, covariate)

write_csv(ame_df, here(plot_dir, "ame_summary.csv"))

# ── 6. Bayesian R² ============================================================

r2_draws_df <- draws_df %>%
  select(starts_with("R2_opt"), starts_with("R2_env"), starts_with("R2_rad")) %>%
  rename_with(~ c("R2_opt", "R2_env", "R2_rad"))

r2_summary <- r2_draws_df %>%
  pivot_longer(cols = everything(), names_to = "block", values_to = "R2") %>%
  group_by(block) %>%
  summarize(mean_R2 = mean(R2), .groups = "drop")

write_csv(r2_summary, here(plot_dir, "r2_summary.csv"))

r2_density <- r2_draws_df %>%
  pivot_longer(cols = everything(), names_to = "block", values_to = "R2") %>%
  ggplot(aes(x = R2, fill = block)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Posterior Densities of Bayesian R²",
    x = expression(R^2),
    y = "Density",
    fill = "Latent Block"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

ggsave(here(plot_dir, "r2_posteriors.png"), r2_density, width = 6, height = 4, dpi = 300)

# ── 7. Group-Level Intercepts =================================================

message("Reconstructing group-level intercepts...")

# Region intercepts
z_alpha_names <- expand.grid(l = 1:3, r = 1:stan_data$R) %>%
  mutate(name = paste0("z_alpha[", l, ",", r, "]")) %>%
  pull(name)

sigma_alpha_names <- paste0("sigma_alpha[", 1:3, "]")

Lcorr_alpha_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_alpha[", i, ",", j, "]")) %>%
  pull(name)

z_alpha_mat     <- draws_df %>% select(all_of(z_alpha_names)) %>% as.matrix()
sigma_alpha_mat <- draws_df %>% select(all_of(sigma_alpha_names)) %>% as.matrix()

Lcorr_alpha_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  L <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_alpha_names) {
    coords <- str_match(nm, "^Lcorr_alpha\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    L[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_alpha_list[[s]] <- L
}

alpha_draws <- array(
  NA_real_,
  dim = c(n_draws, stan_data$R, 3),
  dimnames = list(NULL, NULL, c("phi", "theta", "psi"))
)

for (s in seq_len(n_draws)) {
  D_alpha <- diag(sigma_alpha_mat[s, ])
  L_alpha <- Lcorr_alpha_list[[s]]
  M_alpha <- D_alpha %*% L_alpha
  z_s <- matrix(z_alpha_mat[s, ], nrow = 3, ncol = stan_data$R, byrow = FALSE)
  alpha_draws[s, , ] <- t(M_alpha %*% z_s)
}

# Party intercepts
z_delta_names <- expand.grid(l = 1:3, q = 1:stan_data$Q) %>%
  mutate(name = paste0("z_delta[", l, ",", q, "]")) %>%
  pull(name)

sigma_delta_names <- paste0("sigma_delta[", 1:3, "]")

Lcorr_delta_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_delta[", i, ",", j, "]")) %>%
  pull(name)

z_delta_mat     <- draws_df %>% select(all_of(z_delta_names)) %>% as.matrix()
sigma_delta_mat <- draws_df %>% select(all_of(sigma_delta_names)) %>% as.matrix()

Lcorr_delta_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  Ld <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_delta_names) {
    coords <- str_match(nm, "^Lcorr_delta\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    Ld[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_delta_list[[s]] <- Ld
}

delta_draws <- array(
  NA_real_,
  dim = c(n_draws, stan_data$Q, 3),
  dimnames = list(NULL, NULL, c("phi", "theta", "psi"))
)

for (s in seq_len(n_draws)) {
  D_delta <- diag(sigma_delta_mat[s, ])
  L_delta <- Lcorr_delta_list[[s]]
  M_delta <- D_delta %*% L_delta
  z_s <- matrix(z_delta_mat[s, ], nrow = 3, ncol = stan_data$Q, byrow = FALSE)
  delta_draws[s, , ] <- t(M_delta %*% z_s)
}

# Compute posterior means
alpha_means <- apply(alpha_draws, c(2, 3), mean)
delta_means <- apply(delta_draws, c(2, 3), mean)

region_intercepts_df <- tibble(
  region_id   = 1:stan_data$R,
  region_name = region_levels,
  phi_alpha   = alpha_means[, 1],
  theta_alpha = alpha_means[, 2],
  psi_alpha   = alpha_means[, 3]
)

party_intercepts_df <- tibble(
  party_id    = 1:stan_data$Q,
  party_name  = party_levels,
  phi_delta   = delta_means[, 1],
  theta_delta = delta_means[, 2],
  psi_delta   = delta_means[, 3]
)

write_csv(region_intercepts_df, here(plot_dir, "region_intercepts.csv"))
write_csv(party_intercepts_df, here(plot_dir, "party_intercepts.csv"))

# Region draws long format for plotting
region_draws_long <- map_dfr(
  seq_len(stan_data$R),
  function(r) {
    tibble(
      region_id = r,
      draw = 1:n_draws,
      phi = alpha_draws[, r, 1],
      theta = alpha_draws[, r, 2],
      psi = alpha_draws[, r, 3]
    )
  }
) %>%
  pivot_longer(cols = c(phi, theta, psi), names_to = "latent", values_to = "value") %>%
  mutate(
    region_name = region_levels[region_id],
    latent = case_when(
      latent == "phi"   ~ "φ (Optimism)",
      latent == "theta" ~ "θ (Environment)",
      latent == "psi"   ~ "ψ (Radical-Reform)"
    )
  ) %>%
  select(region_name, latent, value)

# Party draws long format
party_draws_long <- map_dfr(
  seq_len(stan_data$Q),
  function(q) {
    tibble(
      party_id = q,
      draw = 1:n_draws,
      phi = delta_draws[, q, 1],
      theta = delta_draws[, q, 2],
      psi = delta_draws[, q, 3]
    )
  }
) %>%
  pivot_longer(cols = c(phi, theta, psi), names_to = "latent", values_to = "value") %>%
  mutate(
    party_name = party_levels[party_id],
    latent = case_when(
      latent == "phi"   ~ "φ (Optimism)",
      latent == "theta" ~ "θ (Environment)",
      latent == "psi"   ~ "ψ (Radical-Reform)"
    )
  ) %>%
  select(party_name, latent, value)

# Region distribution plot
latent_cols <- c(phi = "#1f78b4", theta = "#33a02c", psi = "#e31a1c")

region_plot_dist <- ggplot(
  region_draws_long,
  aes(x = value, y = fct_rev(region_name), fill = latent, colour = latent)
) +
  stat_halfeye(
    position = position_nudge(y = 0.2),
    slab_alpha = 0.6,
    slab_linewidth = 0.5,
    adjust = 0.8
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  scale_fill_manual(name = "Latent", values = latent_cols) +
  scale_colour_manual(name = "Latent", values = latent_cols) +
  theme_minimal() +
  labs(
    title = "Region-Level Intercepts (Posteriors)",
    subtitle = "Latent traits: φ (Optimism), θ (Environmentalism), ψ (Radical-Reform)",
    x = "Intercept Value",
    y = "Region"
  ) +
  theme(legend.position = "none")

ggsave(here(plot_dir, "region_group_effects_distributions.png"),
       region_plot_dist, width = 10, height = 12, dpi = 300)

# Party distribution plot
party_plot_dist <- ggplot(
  party_draws_long,
  aes(x = value, y = fct_rev(party_name), fill = latent, colour = latent)
) +
  stat_halfeye(
    position = position_nudge(y = 0.2),
    slab_alpha = 0.6,
    slab_linewidth = 0.5,
    adjust = 0.8
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  scale_fill_manual(name = "Latent", values = latent_cols) +
  scale_colour_manual(name = "Latent", values = latent_cols) +
  theme_minimal() +
  labs(
    title = "Party-Level Intercepts (Posteriors)",
    x = "Intercept Value",
    y = "Party"
  ) +
  theme(legend.position = "none")

ggsave(here(plot_dir, "party_group_effects_distributions.png"),
       party_plot_dist, width = 10, height = 12, dpi = 300)

# ── 8. Residual Correlations ==================================================

message("Computing residual correlations...")

Lcorr_eta_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_eta[", i, ",", j, "]")) %>%
  pull(name)

Lcorr_eta_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  Le <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_eta_names) {
    coords <- str_match(nm, "^Lcorr_eta\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    Le[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_eta_list[[s]] <- Le
}

R_eta_list <- lapply(Lcorr_eta_list, function(Le) Le %*% t(Le))
R_eta_array <- simplify2array(R_eta_list)
R_eta_mean <- apply(R_eta_array, c(1, 2), mean)

resid_corr_df <- tibble(
  pair = c("φ–θ", "φ–ψ", "θ–ψ"),
  correlation = c(R_eta_mean[1, 2], R_eta_mean[1, 3], R_eta_mean[2, 3])
)

write_csv(resid_corr_df, here(plot_dir, "residual_correlations.csv"))

# Correlation density plot
corr_long <- map_dfr(R_eta_list, function(Le) {
  tibble(`φ–θ` = Le[1, 2], `φ–ψ` = Le[1, 3], `θ–ψ` = Le[2, 3])
}) %>%
  pivot_longer(cols = everything(), names_to = "pair", values_to = "corr")

point_estimates <- corr_long %>%
  group_by(pair) %>%
  summarize(est = mean(corr)) %>%
  ungroup()

corr_density_with_point <- ggplot(corr_long, aes(x = corr, fill = pair, color = pair)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  geom_vline(
    data = point_estimates,
    aes(xintercept = est, color = pair),
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Latent Correlations: Posterior Densities",
    x = "Correlation",
    y = "Density",
    fill = "Pair",
    color = "Pair"
  ) +
  scale_fill_manual(values = c("φ–θ" = "#666666", "φ–ψ" = "#6a3d9a", "θ–ψ" = "#ff7f00")) +
  scale_color_manual(values = c("φ–θ" = "#666666", "φ–ψ" = "#6a3d9a", "θ–ψ" = "#ff7f00")) +
  theme(legend.position = "top")

ggsave(here(plot_dir, "latent_corr_with_point_estimate.png"),
       corr_density_with_point, width = 8, height = 6, dpi = 300)

# ── 9. Measurement Reliability ================================================

message("Computing measurement reliability...")

lambda_opt_mat <- draws_df %>% select(matches("^lambda_opt\\[")) %>% as.matrix()
sigma_opt_mat  <- draws_df %>% select(matches("^sigma_opt\\[")) %>% as.matrix()
lambda_env_mat <- draws_df %>% select(matches("^lambda_env\\[")) %>% as.matrix()
sigma_env_mat  <- draws_df %>% select(matches("^sigma_env\\[")) %>% as.matrix()
lambda_rad_mat <- draws_df %>% select(matches("^lambda_rad\\[")) %>% as.matrix()
sigma_rad_mat  <- draws_df %>% select(matches("^sigma_rad\\[")) %>% as.matrix()

# McDonald's omega
calc_omega <- function(lambda_mat, sigma_mat) {
  apply(cbind(lambda_mat, sigma_mat), 1, function(row) {
    J <- ncol(lambda_mat)
    lambda <- row[1:J]
    sigma <- row[(J+1):(2*J)]
    num <- sum(lambda)^2
    den <- num + sum(sigma^2)
    num / den
  })
}

omega_opt_draws <- calc_omega(lambda_opt_mat, sigma_opt_mat)
omega_env_draws <- calc_omega(lambda_env_mat, sigma_env_mat)
omega_rad_draws <- calc_omega(lambda_rad_mat, sigma_rad_mat)

omega_summary <- tibble(
  block = c("omega_opt", "omega_env", "omega_rad"),
  mean_omega = c(mean(omega_opt_draws), mean(omega_env_draws), mean(omega_rad_draws))
)

write_csv(omega_summary, here(plot_dir, "reliability_omega.csv"))

# Item discrimination plot
lambda_draws_long <- bind_rows(
  map_dfr(1:stan_data$J_opt, function(j) {
    tibble(item = j, latent = "φ (Optimism)", value = lambda_opt_mat[, j])
  }),
  map_dfr(1:stan_data$J_env, function(j) {
    tibble(item = j, latent = "θ (Environment)", value = lambda_env_mat[, j])
  }),
  map_dfr(1:stan_data$J_rad, function(j) {
    tibble(item = j, latent = "ψ (Radical-Reform)", value = lambda_rad_mat[, j])
  })
)

lambda_dist_plot <- ggplot(
  lambda_draws_long,
  aes(x = value, y = factor(item), fill = latent, colour = latent)
) +
  stat_halfeye(
    position = position_nudge(y = 0.2),
    slab_alpha = 0.6,
    slab_linewidth = 0.5,
    adjust = 0.8
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  scale_y_discrete(drop = TRUE) +
  theme_minimal() +
  labs(
    title = "Item Discrimination (λ): Posterior Distributions",
    x = "Discrimination (λ)",
    y = "Item Index"
  ) +
  theme(legend.position = "none")

ggsave(here(plot_dir, "item_discrimination_lambda_distributions.png"),
       lambda_dist_plot, width = 12, height = 4, dpi = 300)

# ── 10. Regional Portraits ====================================================

message("Creating regional portraits...")

region_counts <- tibble(region_id = stan_data$region_id) %>%
  count(region_id, name = "n_obs")

region_plot_df <- region_intercepts_df %>%
  left_join(region_counts, by = "region_id")

# Pairwise scatterplots
p_psi_phi_col_env <- ggplot(region_plot_df, aes(
  x = psi_alpha, y = phi_alpha, color = theta_alpha, size = n_obs
)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = region_name), size = 3, max.overlaps = 15) +
  scale_color_gradient(low = "grey80", high = "darkgreen", name = "Environmentalism (θ)") +
  scale_size_continuous(range = c(3, 12), name = "# of Respondents") +
  labs(
    title = "Regional profiles in attitude space",
    subtitle = "Radical-Reform (ψ) vs Optimism (φ), coloured by Environmentalism (θ)",
    x = "Radical-Reform ψ (intercept)",
    y = "Optimism φ (intercept)"
  ) +
  theme_minimal(base_size = 13)

ggsave(here(plot_dir, "region_psi_vs_phi_col_env.png"),
       p_psi_phi_col_env, width = 8, height = 6, dpi = 300)

# 3D scatter
fig_3d <- plot_ly(
  data = region_plot_df,
  x = ~phi_alpha,
  y = ~theta_alpha,
  z = ~psi_alpha,
  size = ~n_obs,
  text = ~region_name,
  mode = "markers",
  type = "scatter3d",
  marker = list(
    sizemode = "diameter",
    sizeref = 2 * max(region_plot_df$n_obs) / (20^2),
    opacity = 0.8
  )
) %>%
  layout(
    title = "3D Region Portraits: φ vs θ vs ψ",
    scene = list(
      xaxis = list(title = "Standardised φ (Optimism)"),
      yaxis = list(title = "Standardised θ (Environmentalism)"),
      zaxis = list(title = "Standardised ψ (Radical-Reform)")
    )
  )

saveWidget(fig_3d, here(plot_dir, "region_portraits_3d.html"), selfcontained = TRUE)

# ── 11. Covariate Effects Heatmap =============================================

ame_wide <- slopes_draws_long %>%
  select(latent, covariate, value) %>%
  group_by(latent, covariate) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = latent, values_from = mean)

ame_long <- ame_wide %>%
  pivot_longer(cols = -covariate, names_to = "latent", values_to = "slope")

heatmap_cov <- ggplot(ame_long, aes(x = latent, y = covariate, fill = slope)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name = "Slope") +
  labs(
    title = "Heatmap of Posterior Mean Slopes (AME)",
    x = "Latent Trait",
    y = "Covariate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here(plot_dir, "covariate_effects_heatmap.png"),
       heatmap_cov, width = 8, height = 6, dpi = 300)

# ── 12. Correlation Network ===================================================

corr_mat <- R_eta_mean
colnames(corr_mat) <- c("Optimism", "Environment", "Radical-Reform")
rownames(corr_mat) <- c("Optimism", "Environment", "Radical-Reform")

corr_long_net <- as_tibble(corr_mat, rownames = "trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "corr")

edges <- corr_long_net %>% filter(trait1 != trait2, abs(corr) > 0.2)
graph <- graph_from_data_frame(edges, vertices = tibble(name = rownames(corr_mat)))

corr_network_plot <- ggraph(graph, layout = "circle") +
  geom_edge_link(aes(width = abs(corr), color = corr), alpha = 0.8) +
  geom_node_point(size = 10, color = "#66a61e") +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_color_gradient2(low = "red", mid = "grey80", high = "blue", midpoint = 0, name = "r") +
  scale_edge_width(range = c(0.5, 2), guide = "none") +
  labs(title = "Residual Correlation Network Among Latent Traits") +
  theme_void()

ggsave(here(plot_dir, "latent_correlation_network.png"),
       corr_network_plot, width = 6, height = 6, dpi = 300)

# ── 13. Party Radar Charts ====================================================

party_radar_df <- party_intercepts_df %>%
  select(party_name, phi_delta, theta_delta, psi_delta) %>%
  column_to_rownames(var = "party_name")

radar_plot_list <- map(
  rownames(party_radar_df),
  function(p_name) {
    df <- tibble(
      trait = c("Optimism", "Environment", "Radical-Reform"),
      value = as.numeric(party_radar_df[p_name, ])
    )
    df <- rbind(df, df[1, ])
    ggplot(df, aes(x = trait, y = value, group = 1)) +
      geom_polygon(fill = "#66c2a5", alpha = 0.4) +
      geom_line(color = "#1b9e77", linewidth = 1) +
      geom_point(size = 2) +
      ylim(min(df$value) - 0.5, max(df$value) + 0.5) +
      labs(title = paste0("Party Portrait: ", p_name), y = "Intercept Value") +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")
      )
  }
)

n_parties <- nrow(party_radar_df)
ncols <- 2
nrows <- ceiling(n_parties / ncols)
radar_grid <- plot_grid(plotlist = radar_plot_list, ncol = ncols)

ggsave(here(plot_dir, "party_radar_charts.png"),
       radar_grid, width = 12, height = 6 * nrows / 2, dpi = 300)

message("\n=========================================")
message("Quantities of Interest complete!")
message("All outputs saved to: ", plot_dir)
message("=========================================")
