# ==============================================================================
# Climate Attitudes UK - Bayesian Hierarchical Model
# Script 01: Data Wrangling
#
# Description:
#   Loads survey data, recodes items, and prepares Stan-ready data structure.
#   Constructs three latent trait indicators:
#     - φ (Economic Optimism): 6 Likert items
#     - θ (Environmental Priority): 1 binary + 4 Likert items
#     - ψ (Radical-Reform): 8 binary forced-choice items
#   Plus material insecurity composite and demographic covariates.
#
# Inputs:
#   - data/synthetic/synthetic_survey_data.rds (synthetic) OR
#   - data/LfG Data.xlsx (proprietary, not included)
#
# Outputs:
#   - data/stan_data_full.rds
#
# Author: Henry Baker
# ==============================================================================

library(tidyverse)
library(readxl)
library(here)

# ── Helpers ──────────────────────────────────────────────────────────────────
recode_map <- function(x, map) map[x]
z_score    <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# ── 0. Data Loading ===========================================================

# Check for real data first, fall back to synthetic
real_data_path <- here("data", "LfG Data.xlsx")
synthetic_data_path <- here("data", "synthetic", "synthetic_survey_data.rds")

if (file.exists(real_data_path)) {
  message("Loading real survey data from: ", real_data_path)
  data <- read_excel(real_data_path)
} else if (file.exists(synthetic_data_path)) {
  message("Loading synthetic data from: ", synthetic_data_path)
  data <- readRDS(synthetic_data_path)
} else {
  stop("No data found. Run 00_generate_synthetic_data.R first or provide LfG Data.xlsx")
}

data <- data %>% mutate(obs_id = row_number())

# ── 1. Economic Optimism (φᵢ) ─────────────────────────────────────────────────

# Identify the six optimism columns (Q43–Q48)
opt_cols <- grep(
  "^To what extent do you feel optimistic or pessimistic about the following",
  names(data),
  value = TRUE
)

# Recode text → numeric on {−2,−1,0,1,2}
opt_map <- c(
  "Very pessimistic"                   = -2L,
  "Fairly pessimistic"                 = -1L,
  "Neither optimistic nor pessimistic" =  0L,
  "Fairly optimistic"                  =  1L,
  "Very optimistic"                    =  2L
)

# Step 1: Create raw numeric columns
data <- data %>%
  mutate(across(all_of(opt_cols),
                ~ recode_map(., opt_map),
                .names = "num_{.col}"
  ))

# Step 2: Z-score each of those six
num_opt_cols <- paste0("num_", opt_cols)
data <- data %>%
  mutate(across(all_of(num_opt_cols),
                ~ z_score(.),
                .names = "z_{.col}"
  ))

# Now gather into long format
z_opt_cols <- paste0("z_", num_opt_cols)

opt_long_df <- data %>%
  select(obs_id, all_of(z_opt_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_opt",
    values_to = "y_opt"
  ) %>%
  filter(!is.na(y_opt)) %>%
  mutate(
    j_opt = as.integer(factor(item_opt, levels = z_opt_cols))
  ) %>%
  select(obs_id, j_opt, y_opt)

# Compute lower and upper bounds for each of the six optimism items
lower_opt <- sapply(opt_cols, function(col) {
  zname <- paste0("z_num_", col)
  min(data[[zname]], na.rm = TRUE)
})
upper_opt <- sapply(opt_cols, function(col) {
  zname <- paste0("z_num_", col)
  max(data[[zname]], na.rm = TRUE)
})

# ── 2. Environmental Priority (θᵢ) ────────────────────────────────────────────

#  2a) Forced choice Q121
env_forced_col <- grep("\\.\\.\\.121$", names(data), value = TRUE)

#  2b) Four Likert items Q122–Q125
env_likert_cols <- grep(
  "^To what extent do you agree or disagree with the following statements",
  names(data),
  value = TRUE
)

# Recode Likert {Strongly disagree → -2, …, Strongly agree → +2}
env_map <- c(
  "Strongly disagree"           = -2L,
  "Somewhat disagree"           = -1L,
  "Neither agree nor disagree"  =  0L,
  "Somewhat agree"              =  1L,
  "Strongly agree"              =  2L
)

# Which of those four needs reversing? (the "cost-first" statement is the 2nd)
anti_col <- env_likert_cols[2]

data <- data %>%
  # 1) Recode all four to {−2,…,+2}
  mutate(across(all_of(env_likert_cols),
                ~ recode_map(., env_map),
                .names = "num_{.col}"
  )) %>%
  # 2) Reverse-score the anti-environment item
  mutate(
    !!paste0("num_", anti_col) := -get(paste0("num_", anti_col))
  ) %>%
  # 3) Recode the forced-choice into {0,1}
  mutate(
    num_forced = case_when(
      .data[[env_forced_col]] ==
        "We should prioritise protecting the natural environment, even if that sometimes prevents action to reduce the cost of living" ~ 1L,
      .data[[env_forced_col]] ==
        "We should prioritise action to reduce the cost of living, even if that sometimes comes at the cost of the environment"      ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  # 4) Z-score all five (forced + four recoded Likert)
  mutate(across(
    c("num_forced", paste0("num_", env_likert_cols)),
    ~ z_score(.),
    .names = "z_{.col}"
  ))

# Build long format for environment
num_env_cols <- c("num_forced", paste0("num_", env_likert_cols))
z_env_cols   <- paste0("z_", num_env_cols)

env_long_df <- data %>%
  select(obs_id, all_of(z_env_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_env",
    values_to = "y_env"
  ) %>%
  filter(!is.na(y_env)) %>%
  mutate(
    j_env = as.integer(factor(item_env, levels = z_env_cols))
  ) %>%
  select(obs_id, j_env, y_env)

# Compute lower and upper bounds for environment items
lower_env <- sapply(z_env_cols, function(col) {
  min(data[[col]], na.rm = TRUE)
})
upper_env <- sapply(z_env_cols, function(col) {
  max(data[[col]], na.rm = TRUE)
})

# ── 3. Radical-Reform (ψᵢ) ────────────────────────────────────────────────────

# The eight forced-choice questions Q72–Q75, Q78–Q79, Q82–Q83
rad_cols <- colnames(data)[c(72, 73, 74, 75, 78, 79, 82, 83)]

# Mapping: 1 = radical statement, 0 = status-quo statement
radical_map_list <- list(
  Q72 = c(
    "Westminster politicians have moved too slowly on initiatives to help the economy" = 1L,
    "Westminster politicians have generally moved at a good pace to get things done to help the economy" = 0L
  ),
  Q73 = c(
    "I am looking for politicians who show that they can get things done" = 1L,
    "I am looking for politicians who take their time and work for long term goals even if that means there is less getting done in the short term" = 0L
  ),
  Q74 = c(
    "I would be more likely to support politicians who show they are prepared to take radical action to improve everyday people's lives" = 1L,
    "I would be more likely to support politicians who do not want to change the system and keep things broadly the same with no surprises" = 0L
  ),
  Q75 = c(
    "Britain has had a rough few years but is broadly doing fine and just needs to have a steady few years" = 0L,
    "Britain is broken and needs radical action to fix it" = 1L
  ),
  Q78 = c(
    "The UK is not taking the right steps to achieve more economic growth and raise living standards" = 1L,
    "The UK is taking the right steps to achieve more economic growth and raise living standards" = 0L
  ),
  Q79 = c(
    "The UK is broadly speaking on the wrong track and needs radical reform to deliver a good quality of life for its citizens" = 1L,
    "The UK is broadly speaking headed in the right direction and should continue as normal to deliver a good quality of life for its citizens" = 0L
  ),
  Q82 = c(
    "It is acceptable for unelected bodies and not elected representatives to make the final decisions about things that affect the UK" = 1L,
    "Ultimately final decisions about things that affect the UK should always lie with an elected representative" = 0L
  ),
  Q83 = c(
    "Politicians have allowed unelected bodies and officials to make too many decisions to shield them from accountability" = 1L,
    "Politicians have taken too many decisions and should trust experts and unelected officials to do what is best for the country" = 0L
  )
)

# Map to column names
qs <- c("Q72", "Q73", "Q74", "Q75", "Q78", "Q79", "Q82", "Q83")
radical_map_list2 <- setNames(radical_map_list[qs], rad_cols)

# Recode each column
for (col in rad_cols) {
  recode_vec <- radical_map_list2[[col]]
  data <- data %>%
    mutate(
      !!paste0("num_", col) := recode_map(.data[[col]], recode_vec)
    )
}

num_rad_cols <- paste0("num_", rad_cols)

# Z-score each radical item
data <- data %>%
  mutate(across(
    all_of(num_rad_cols),
    ~ z_score(.),
    .names = "z_{.col}"
  ))

z_rad_cols <- paste0("z_num_", rad_cols)

# Assemble into long format
rad_long_df <- data %>%
  select(obs_id, all_of(z_rad_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_rad",
    values_to = "y_rad"
  ) %>%
  filter(!is.na(y_rad)) %>%
  mutate(
    j_rad = as.integer(factor(item_rad, levels = z_rad_cols))
  ) %>%
  select(obs_id, j_rad, y_rad)

# Compute lower and upper bounds for radical items
lower_rad <- sapply(z_rad_cols, function(col) {
  min(data[[col]], na.rm = TRUE)
})
upper_rad <- sapply(z_rad_cols, function(col) {
  max(data[[col]], na.rm = TRUE)
})

# ── 4. Material-Insecurity (Mᵢ) ───────────────────────────────────────────────

# Q90–Q97 frequency items
insec_cols <- grep(
  "^For each of the following, say how they apply or do not apply to you:",
  names(data),
  value = TRUE
)

insec_map <- c(
  "Never"        = 0L,
  "Occasionally" = 1L,
  "Fairly often" = 2L,
  "Very often"   = 3L
)

# 1) Recode each to 0–3
data <- data %>%
  mutate(across(
    all_of(insec_cols),
    ~ recode_map(., insec_map),
    .names = "num_{.col}"
  )) %>%
  # 2) Build composite (row mean) and then standardise
  mutate(
    M_full     = rowMeans(select(., paste0("num_", insec_cols)), na.rm = TRUE),
    insecurity = z_score(M_full)
  )

# ── 5. Demographics & Covariates (Xᵢ) ─────────────────────────────────────────

#  5a) Combine "vote tomorrow" / "vote recent" into one
vote_tomorrow <- grep(
  "Which party would you vote for if there were a general election tomorrow",
  names(data), value = TRUE, ignore.case = TRUE
)
vote_recent <- grep(
  "Which party did you vote for in the most recent general election",
  names(data), value = TRUE, ignore.case = TRUE
)

data <- data %>%
  mutate(
    vote_tomorrow = na_if(.data[[vote_tomorrow]], ""),
    vote_recent   = na_if(.data[[vote_recent]], ""),
    vote          = coalesce(vote_tomorrow, vote_recent)
  )

#  5b) Define party categories
party_levels <- c(
  "Green", "Labour", "Plaid Cymru",
  "Scottish National Party (SNP)",
  "Liberal Democrat", "Conservative", "Reform UK",
  "Another party", "Dont know", "Won't vote"
)
data <- data %>%
  mutate(
    party = factor(vote, levels = party_levels),
    party_id = as.integer(party)
  )

#  5c) Region → integer 1..10
region_col <- grep("^Region$", names(data), value = TRUE, ignore.case = TRUE)
region_levels <- c(
  "Yorkshire & the Humber", "West Midlands", "Scotland", "Wales",
  "North West", "Eastern", "South West", "East Midlands", "London", "South East"
)
data <- data %>%
  mutate(
    region = factor(.data[[region_col]], levels = region_levels),
    region_id = as.integer(region)
  )

#  5d) Gender → binary dummy (0 = Male, 1 = Female)
gender_col <- grep("^Gender$", names(data), value = TRUE, ignore.case = TRUE)
data <- data %>%
  mutate(
    gender = case_when(
      .data[[gender_col]] == "Male"   ~ 0L,
      .data[[gender_col]] == "Female" ~ 1L,
      TRUE                            ~ NA_integer_
    )
  )

#  5e) Age bracket (6 levels → 5 dummies, reference = "18–24")
age_bracket_cols_both <- grep("^Age.*\\d+$", names(data), value = TRUE)
age_bracket_col <- age_bracket_cols_both[2]
data <- data %>%
  mutate(
    age_raw = factor(.data[[age_bracket_col]],
                     levels = c("65+", "55-64", "45-54", "25-34", "35-44", "18-24"))
  ) %>%
  mutate(
    age_65plus  = as.integer(age_raw == "65+"),
    age_55_64   = as.integer(age_raw == "55-64"),
    age_45_54   = as.integer(age_raw == "45-54"),
    age_25_34   = as.integer(age_raw == "25-34"),
    age_35_44   = as.integer(age_raw == "35-44")
  )

#  5f) Education level (7 categories → 6 dummies, reference = "No qualifications")
edu_col <- grep("educational.*qualification", names(data), value = TRUE, ignore.case = TRUE)
edu_levels <- c(
  "Level 4 qualifications or above, for example a bachelor's degree or above",
  "Level 2 qualifications, for example 5 or more GCSE passes (formerly O levels)",
  "Level 1 and entry level qualifications, for example 1 to 4 GCSE passes (formerly O levels)",
  "Level 3 qualifications, for example 2 or more A levels",
  "Apprenticeship",
  "No qualifications",
  "Other qualifications"
)
data <- data %>%
  mutate(
    edu_raw = factor(.data[[edu_col]], levels = edu_levels)
  ) %>%
  mutate(
    edu_L4plus  = as.integer(edu_raw == edu_levels[1]),
    edu_L2      = as.integer(edu_raw == edu_levels[2]),
    edu_L1      = as.integer(edu_raw == edu_levels[3]),
    edu_L3      = as.integer(edu_raw == edu_levels[4]),
    edu_appr    = as.integer(edu_raw == edu_levels[5]),
    edu_other   = as.integer(edu_raw == edu_levels[7])
  )

# ── 6. Assemble Covariate DataFrame (cov_df) ──────────────────────────────────

cov_df <- data %>%
  select(
    obs_id,
    party_id,
    region_id,
    gender,
    age_65plus, age_55_64, age_45_54, age_35_44, age_25_34,
    edu_L4plus, edu_L3, edu_L2, edu_L1, edu_appr, edu_other,
    insecurity
  )

# Drop rows where any covariate is NA
cov_df <- cov_df %>% drop_na()

# ── 7. Re-index respondents and merge with long responses ─────────────────────

# 7a) Determine which obs_id remain
valid_ids <- cov_df$obs_id

# 7b) Re-label them to 1..N
cov_df <- cov_df %>%
  arrange(obs_id) %>%
  mutate(i = row_number())

# 7c) Build a lookup table: old obs_id → new i
id_map <- tibble(old = cov_df$obs_id, new = cov_df$i)

# 7d) Filter and remap opt_long_df
opt_long_df2 <- opt_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_opt, y_opt)

# 7e) Filter and remap env_long_df
env_long_df2 <- env_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_env, y_env)

# 7f) Filter and remap rad_long_df
rad_long_df2 <- rad_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_rad, y_rad)

# Check how many subjects remain
N_final <- nrow(cov_df)
message("Final sample size: N = ", N_final)

# ── 8. Extract region_id and party_id vectors ─────────────────────────────────

resp_region <- cov_df$region_id
resp_party  <- cov_df$party_id

# ── 9. Build design matrix X (N_final × 13) ───────────────────────────────────

Xmat <- cov_df %>%
  select(
    gender,
    age_65plus, age_55_64, age_45_54, age_35_44, age_25_34,
    edu_L4plus, edu_L3, edu_L2, edu_L1, edu_appr, edu_other,
    insecurity
  ) %>%
  as.matrix()

P_cov <- ncol(Xmat)
message("Number of covariates: P = ", P_cov)

# ── 10. Final sanity checks ───────────────────────────────────────────────────

stopifnot(range(opt_long_df2$j_opt) == c(1, 6))
stopifnot(range(env_long_df2$j_env) == c(1, 5))
stopifnot(range(rad_long_df2$j_rad) == c(1, 8))
stopifnot(range(opt_long_df2$i) == c(1, N_final))
stopifnot(range(env_long_df2$i) == c(1, N_final))
stopifnot(range(rad_long_df2$i) == c(1, N_final))

message("All sanity checks passed!")

# ── 11. Prepare the stan_data list ────────────────────────────────────────────

stan_data <- list(
  N         = N_final,
  P         = P_cov,
  X         = Xmat,
  R         = length(unique(resp_region)),
  region_id = resp_region,
  Q         = max(resp_party),
  party_id  = resp_party,

  J_opt     = length(unique(opt_long_df2$j_opt)),
  N_opt     = nrow(opt_long_df2),
  i_opt     = opt_long_df2$i,
  j_opt     = opt_long_df2$j_opt,
  y_opt     = opt_long_df2$y_opt,
  lower_opt = as.vector(lower_opt),
  upper_opt = as.vector(upper_opt),

  J_env     = length(unique(env_long_df2$j_env)),
  N_env     = nrow(env_long_df2),
  i_env     = env_long_df2$i,
  j_env     = env_long_df2$j_env,
  y_env     = env_long_df2$y_env,
  lower_env = as.vector(lower_env),
  upper_env = as.vector(upper_env),

  J_rad     = length(unique(rad_long_df2$j_rad)),
  N_rad     = nrow(rad_long_df2),
  i_rad     = rad_long_df2$i,
  j_rad     = rad_long_df2$j_rad,
  y_rad     = rad_long_df2$y_rad,
  lower_rad = as.vector(lower_rad),
  upper_rad = as.vector(upper_rad)
)

# Save data for Stan
output_path <- here("data", "stan_data_full.rds")
saveRDS(stan_data, output_path)
message("Stan data saved to: ", output_path)
