// =============================================================================
// Hierarchical Latent Trait Model for UK Climate Attitudes
// Three-dimensional measurement model: φ (Optimism), θ (Environment), ψ (Radical-Reform)
// With region and party random intercepts
// =============================================================================

data {
  int<lower=1> N;                // # of individuals
  int<lower=1> P;                // # of numeric covariates
  matrix[N, P] X;                // standardised covariates

  int<lower=1> R;                           // # of regions
  array[N] int<lower=1, upper=R> region_id; // region index

  int<lower=1> Q;                          // # of parties
  array[N] int<lower=1, upper=Q> party_id; // party index

  int<lower=1> J_opt;                          // # of optimism items
  int<lower=1> N_opt;                          // total optimism responses
  array[N_opt] int<lower=1, upper=N>    i_opt;  // respondent index
  array[N_opt] int<lower=1, upper=J_opt> j_opt; // item index
  array[N_opt] real                     y_opt; // observed (standardised)

  int<lower=1> J_env;                          // # of environment items
  int<lower=1> N_env;                          // total environment responses
  array[N_env] int<lower=1, upper=N>    i_env;  // respondent index
  array[N_env] int<lower=1, upper=J_env> j_env; // item index
  array[N_env] real                     y_env; // observed (standardised)

  int<lower=1> J_rad;                          // # of radical-reform items
  int<lower=1> N_rad;                          // total radical-reform responses
  array[N_rad] int<lower=1, upper=N>    i_rad;  // respondent index
  array[N_rad] int<lower=1, upper=J_rad> j_rad; // item index
  array[N_rad] real                     y_rad; // observed (standardised)
}

parameters {
  // A) Latent-factor hierarchy (3-dimensional: φ, θ, ψ)
  cholesky_factor_corr[3] Lcorr_eta;  // corr(φ, θ, ψ)
  vector<lower=0>[3]      tau_eta;    // half-Normal(0, 0.3)
  matrix[3, N]            z_eta;      // non-centred

  // 2) Region intercepts (3-dimensional)
  cholesky_factor_corr[3] Lcorr_alpha; // corr across (α₁, α₂, α₃)
  matrix[3, R]            z_alpha;     // non-centred
  vector<lower=0>[3]      sigma_alpha; // SD ≥ 0

  // 3) Party intercepts (3-dimensional)
  cholesky_factor_corr[3] Lcorr_delta; // corr across (δ₁, δ₂, δ₃)
  matrix[3, Q]            z_delta;     // non-centred
  vector<lower=0>[3]      sigma_delta; // SD ≥ 0

  // 4) Covariate slopes (3 × P)
  matrix[3, P]            B; // Normal(0, 0.5)

  // B) Measurement: optimism items
  vector[J_opt]           beta_opt;    // intercepts
  vector<lower=0>[J_opt]  lambda_opt;  // loadings ≥ 0
  vector<lower=0>[J_opt]  sigma_opt;   // residual SD ≥ 0

  // C) Measurement: environment items
  vector[J_env]           beta_env;    // intercepts
  vector<lower=0>[J_env]  lambda_env;  // loadings ≥ 0
  vector<lower=0>[J_env]  sigma_env;   // residual SD ≥ 0

  // D) Measurement: radical-reform items
  vector[J_rad]           beta_rad;    // intercepts
  vector<lower=0>[J_rad]  lambda_rad;  // loadings ≥ 0
  vector<lower=0>[J_rad]  sigma_rad;   // residual SD ≥ 0
}

transformed parameters {
  // Expand region & party covariance matrices
  cov_matrix[3] Sigma_alpha =
    diag_pre_multiply(sigma_alpha, Lcorr_alpha)
    * diag_pre_multiply(sigma_alpha, Lcorr_alpha)';
  cov_matrix[3] Sigma_delta =
    diag_pre_multiply(sigma_delta, Lcorr_delta)
    * diag_pre_multiply(sigma_delta, Lcorr_delta)';

  // Vectors for each latent dimension
  vector[N] phi;
  vector[N] theta;
  vector[N] psi;

  for (i in 1:N) {
    // region + party intercept contributions (3-vector)
    vector[3] mu_eta_i =
      diag_pre_multiply(sigma_alpha, Lcorr_alpha) * z_alpha[, region_id[i]] +
      diag_pre_multiply(sigma_delta, Lcorr_delta) * z_delta[, party_id[i]] +
      B * to_vector(X[i]);

    // latent noise (3-vector)
    vector[3] noise =
      diag_pre_multiply(tau_eta, Lcorr_eta) * z_eta[, i];

    vector[3] eta_i = mu_eta_i + noise;
    phi[i]   = eta_i[1];
    theta[i] = eta_i[2];
    psi[i]   = eta_i[3];
  }
}

model {
  // 1) Priors on τ_eta (latent SDs)
  // Optimism (φ) and Environment (θ): existing tuned priors
  tau_eta[1] ~ normal(0, 0.3);
  tau_eta[2] ~ normal(0, 0.3);

  // Radical Reform (ψ): tighter prior
  tau_eta[3] ~ normal(0, 0.1);

  Lcorr_eta        ~ lkj_corr_cholesky(2.0);
  to_vector(z_eta) ~ normal(0, 1);

  // 2) Region intercept priors
  Lcorr_alpha      ~ lkj_corr_cholesky(2.0);
  sigma_alpha      ~ normal(0, 0.1) T[0, ];
  to_vector(z_alpha) ~ normal(0, 1);

  // 3) Party intercept priors
  Lcorr_delta      ~ lkj_corr_cholesky(2.0);
  sigma_delta      ~ normal(0, 0.1) T[0, ];
  to_vector(z_delta) ~ normal(0, 1);

  // 4) Covariate slopes
  to_vector(B)     ~ normal(0, 0.5);

  // 5) Measurement: optimism
  beta_opt        ~ normal(0, 0.5);
  lambda_opt      ~ lognormal(log(1), 0.2);
  sigma_opt       ~ normal(1, 0.2);

  // 6) Measurement: environment
  beta_env        ~ normal(0, 0.5);
  lambda_env      ~ lognormal(log(1), 0.2);
  sigma_env       ~ normal(1, 0.2);

  // 7) Measurement: radical-reform
  beta_rad        ~ normal(0, 0.5);
  lambda_rad      ~ lognormal(log(1), 0.2);
  sigma_rad       ~ normal(1, 0.2);

  // 8) Likelihood: y_opt
  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    y_opt[n] ~ normal(mu_opt, sigma_opt[jj]);
  }

  // 9) Likelihood: y_env
  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env = beta_env[jj] + lambda_env[jj] * theta[ii];
    y_env[n] ~ normal(mu_env, sigma_env[jj]);
  }

  // 10) Likelihood: y_rad
  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    y_rad[n] ~ normal(mu_rad, sigma_rad[jj]);
  }
}

generated quantities {
  // A) Posterior-predictive y's for each block
  vector[N_opt] y_opt_sim;
  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    y_opt_sim[n] = normal_rng(mu_opt, sigma_opt[jj]);
  }

  vector[N_env] y_env_sim;
  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env = beta_env[jj] + lambda_env[jj] * theta[ii];
    y_env_sim[n] = normal_rng(mu_env, sigma_env[jj]);
  }

  vector[N_rad] y_rad_sim;
  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    y_rad_sim[n] = normal_rng(mu_rad, sigma_rad[jj]);
  }

  // B) R² computations for each latent block
  // 1) R²_opt
  vector[N_opt] yhat_opt_resp;
  vector[N_opt] var_resid_opt_resp;

  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt_n = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    yhat_opt_resp[n]      = mu_opt_n;
    var_resid_opt_resp[n] = square(sigma_opt[jj]);
  }

  real Var_pred_opt = variance(yhat_opt_resp);
  real E_resid_opt  = mean(var_resid_opt_resp);
  real R2_opt       = Var_pred_opt / (Var_pred_opt + E_resid_opt);

  // 2) R²_env
  vector[N_env] yhat_env_resp;
  vector[N_env] var_resid_env_resp;

  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env_n = beta_env[jj] + lambda_env[jj] * theta[ii];
    yhat_env_resp[n]      = mu_env_n;
    var_resid_env_resp[n] = square(sigma_env[jj]);
  }

  real Var_pred_env = variance(yhat_env_resp);
  real E_resid_env  = mean(var_resid_env_resp);
  real R2_env       = Var_pred_env / (Var_pred_env + E_resid_env);

  // 3) R²_rad
  vector[N_rad] yhat_rad_resp;
  vector[N_rad] var_resid_rad_resp;

  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad_n = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    yhat_rad_resp[n]      = mu_rad_n;
    var_resid_rad_resp[n] = square(sigma_rad[jj]);
  }

  real Var_pred_rad = variance(yhat_rad_resp);
  real E_resid_rad  = mean(var_resid_rad_resp);
  real R2_rad       = Var_pred_rad / (Var_pred_rad + E_resid_rad);
}
