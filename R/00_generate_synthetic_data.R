# ==============================================================================
# Climate Attitudes UK - Bayesian Hierarchical Model
# Script 00: Generate Synthetic Data
#
# Description:
#   Generates synthetic survey data matching the structure expected by
#   01_wrangling.R. This allows the analysis pipeline to be tested and
#   demonstrated without access to the proprietary survey data.
#
#   The original data is from the Looking for Growth tracker
#   (https://tracker.lookingforgrowth.uk/) and is not publicly available.
#
# Outputs:
#   - data/synthetic/synthetic_survey_data.rds
#
# Author: Henry Baker
# ==============================================================================

library(tidyverse)
library(here)
library(MASS)  # For mvrnorm

set.seed(42)

# ── Configuration =============================================================

N_SYNTH <- 500  # Number of synthetic respondents

# Response option labels
opt_labels <- c(
  "Very pessimistic",
  "Fairly pessimistic",
  "Neither optimistic nor pessimistic",
  "Fairly optimistic",
  "Very optimistic"
)

agree_labels <- c(
  "Strongly disagree",
  "Somewhat disagree",
  "Neither agree nor disagree",
  "Somewhat agree",
  "Strongly agree"
)

insec_labels <- c("Never", "Occasionally", "Fairly often", "Very often")

# Demographics
region_levels <- c(
  "Yorkshire & the Humber", "West Midlands", "Scotland", "Wales",
  "North West", "Eastern", "South West", "East Midlands", "London", "South East"
)

party_levels <- c(
  "Green", "Labour", "Plaid Cymru", "Scottish National Party (SNP)",
  "Liberal Democrat", "Conservative", "Reform UK",
  "Another party", "Dont know", "Won't vote"
)

age_levels <- c("18-24", "25-34", "35-44", "45-54", "55-64", "65+")

edu_levels <- c(
  "Level 4 qualifications or above, for example a bachelor's degree or above",
  "Level 2 qualifications, for example 5 or more GCSE passes (formerly O levels)",
  "Level 1 and entry level qualifications, for example 1 to 4 GCSE passes (formerly O levels)",
  "Level 3 qualifications, for example 2 or more A levels",
  "Apprenticeship",
  "No qualifications",
  "Other qualifications"
)

# ── Helper Functions ==========================================================

# Convert continuous latent to ordinal response
latent_to_ordinal <- function(latent, n_cats, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- qnorm(seq(1/n_cats, (n_cats-1)/n_cats, length.out = n_cats - 1))
  }
  cut_vals <- c(-Inf, thresholds, Inf)
  as.integer(cut(latent, breaks = cut_vals, labels = FALSE))
}

# ── Generate Latent Traits ====================================================

message("Generating correlated latent traits...")

# Population correlation structure for latent traits
# φ (Optimism), θ (Environment), ψ (Radical-Reform)
latent_corr <- matrix(c(
  1.0,  0.1, -0.2,
  0.1,  1.0,  0.3,
 -0.2,  0.3,  1.0
), nrow = 3, byrow = TRUE)

# Generate latent trait scores
latent_traits <- mvrnorm(n = N_SYNTH, mu = c(0, 0, 0), Sigma = latent_corr)
colnames(latent_traits) <- c("phi", "theta", "psi")
latent_traits <- as_tibble(latent_traits)

# Add material insecurity (correlated with traits)
latent_traits$insecurity <- -0.3 * latent_traits$phi + 0.1 * latent_traits$psi + rnorm(N_SYNTH, 0, 0.8)

# ── Generate Demographics =====================================================

message("Generating demographics...")

# Region (weighted by UK population)
region_probs <- c(0.08, 0.09, 0.08, 0.05, 0.11, 0.10, 0.09, 0.07, 0.13, 0.15)
region_probs <- region_probs / sum(region_probs)
regions <- sample(region_levels, N_SYNTH, replace = TRUE, prob = region_probs)

# Party (synthetic distribution)
party_probs <- c(0.08, 0.30, 0.02, 0.04, 0.10, 0.22, 0.08, 0.03, 0.08, 0.05)
party_probs <- party_probs / sum(party_probs)
parties <- sample(party_levels, N_SYNTH, replace = TRUE, prob = party_probs)

# Gender
gender <- sample(c("Male", "Female"), N_SYNTH, replace = TRUE)

# Age brackets
age_probs <- c(0.10, 0.15, 0.18, 0.18, 0.17, 0.22)
age <- sample(age_levels, N_SYNTH, replace = TRUE, prob = age_probs)

# Education
edu_probs <- c(0.35, 0.15, 0.10, 0.15, 0.08, 0.10, 0.07)
education <- sample(edu_levels, N_SYNTH, replace = TRUE, prob = edu_probs)

# ── Generate Item Responses ===================================================

message("Generating item responses...")

# Create base data frame
synth_data <- tibble(
  obs_id = 1:N_SYNTH,
  Region = regions,
  Gender = gender
)

# Add age columns (the script expects specific column positions)
synth_data$`Age...1` <- age  # First age column
synth_data$`Age...2` <- age  # Second age column (script uses this one)

# Add education
synth_data$`What is your highest educational qualification?` <- education

# Add party vote columns
synth_data$`Which party would you vote for if there were a general election tomorrow?` <- parties
synth_data$`Which party did you vote for in the most recent general election?` <- parties

# ── Optimism Items (Q43-Q48) ==================================================

opt_item_names <- paste0(
  "To what extent do you feel optimistic or pessimistic about the following_Q",
  43:48
)

# Generate 6 optimism items based on φ latent trait
for (i in 1:6) {
  item_latent <- latent_traits$phi + rnorm(N_SYNTH, 0, 0.7)  # Add measurement error
  item_cat <- latent_to_ordinal(item_latent, 5)
  synth_data[[opt_item_names[i]]] <- opt_labels[item_cat]
}

# ── Environment Items (Q121-Q125) =============================================

# Q121: Forced choice
env_forced_col <- "Which of the following statements comes closest to your view?...121"
env_forced_vals <- ifelse(
  latent_traits$theta + rnorm(N_SYNTH, 0, 0.5) > 0,
  "We should prioritise protecting the natural environment, even if that sometimes prevents action to reduce the cost of living",
  "We should prioritise action to reduce the cost of living, even if that sometimes comes at the cost of the environment"
)
synth_data[[env_forced_col]] <- env_forced_vals

# Q122-Q125: Likert agreement items
env_likert_names <- paste0(
  "To what extent do you agree or disagree with the following statements_Q",
  122:125
)

# Generate 4 environment Likert items (Q122 is reverse-coded in original)
for (i in 1:4) {
  if (i == 2) {
    # This item is reverse-coded (anti-environment statement)
    item_latent <- -latent_traits$theta + rnorm(N_SYNTH, 0, 0.7)
  } else {
    item_latent <- latent_traits$theta + rnorm(N_SYNTH, 0, 0.7)
  }
  item_cat <- latent_to_ordinal(item_latent, 5)
  synth_data[[env_likert_names[i]]] <- agree_labels[item_cat]
}

# ── Radical-Reform Items (Q72-Q83) ============================================

# Define the forced-choice pairs (Q72, Q73, Q74, Q75, Q78, Q79, Q82, Q83)
radical_items <- list(
  list(
    name = "For the following pairs...Q72",
    radical = "Westminster politicians have moved too slowly on initiatives to help the economy",
    status_quo = "Westminster politicians have generally moved at a good pace to get things done to help the economy"
  ),
  list(
    name = "For the following pairs...Q73",
    radical = "I am looking for politicians who show that they can get things done",
    status_quo = "I am looking for politicians who take their time and work for long term goals even if that means there is less getting done in the short term"
  ),
  list(
    name = "For the following pairs...Q74",
    radical = "I would be more likely to support politicians who show they are prepared to take radical action to improve everyday people's lives",
    status_quo = "I would be more likely to support politicians who do not want to change the system and keep things broadly the same with no surprises"
  ),
  list(
    name = "For the following pairs...Q75",
    radical = "Britain is broken and needs radical action to fix it",
    status_quo = "Britain has had a rough few years but is broadly doing fine and just needs to have a steady few years"
  ),
  list(
    name = "For the following pairs...Q78",
    radical = "The UK is not taking the right steps to achieve more economic growth and raise living standards",
    status_quo = "The UK is taking the right steps to achieve more economic growth and raise living standards"
  ),
  list(
    name = "For the following pairs...Q79",
    radical = "The UK is broadly speaking on the wrong track and needs radical reform to deliver a good quality of life for its citizens",
    status_quo = "The UK is broadly speaking headed in the right direction and should continue as normal to deliver a good quality of life for its citizens"
  ),
  list(
    name = "For the following pairs...Q82",
    radical = "It is acceptable for unelected bodies and not elected representatives to make the final decisions about things that affect the UK",
    status_quo = "Ultimately final decisions about things that affect the UK should always lie with an elected representative"
  ),
  list(
    name = "For the following pairs...Q83",
    radical = "Politicians have allowed unelected bodies and officials to make too many decisions to shield them from accountability",
    status_quo = "Politicians have taken too many decisions and should trust experts and unelected officials to do what is best for the country"
  )
)

# Generate radical-reform responses based on ψ latent trait
for (item in radical_items) {
  item_latent <- latent_traits$psi + rnorm(N_SYNTH, 0, 0.6)
  item_choice <- ifelse(item_latent > 0, item$radical, item$status_quo)
  synth_data[[item$name]] <- item_choice
}

# ── Insecurity Items (Q90-Q97) ================================================

insec_item_names <- paste0(
  "For each of the following, say how they apply or do not apply to you:_Q",
  90:97
)

# Generate 8 insecurity items based on insecurity latent
for (i in 1:8) {
  item_latent <- latent_traits$insecurity + rnorm(N_SYNTH, 0, 0.6)
  # Insecurity items use 4-point scale (Never to Very often)
  item_cat <- latent_to_ordinal(item_latent, 4)
  synth_data[[insec_item_names[i]]] <- insec_labels[item_cat]
}

# ── Ensure column positions match expected indices ============================

# The original script uses column positions [72, 73, 74, 75, 78, 79, 82, 83]
# for radical items. We need to pad the data frame to match.

# Current columns
current_cols <- names(synth_data)

# Create a data frame with placeholder columns to get radical items at right positions
placeholder_df <- synth_data

# Add placeholder columns to push radical items to positions 72-83
# First, rename radical columns to simple names
radical_cols_old <- sapply(radical_items, function(x) x$name)

# We need columns 1-71 to be non-radical, then radical items at 72, 73, 74, 75
# then some columns, then 78, 79, then more, then 82, 83

# Simpler approach: just ensure the column names match the grep pattern
# The script uses: colnames(data)[c(72, 73, 74, 75, 78, 79, 82, 83)]
# So we need to ensure these positions have the radical columns

# Create ordered column list
base_cols <- c("obs_id", "Region", "Gender", "Age...1", "Age...2",
               "What is your highest educational qualification?",
               "Which party would you vote for if there were a general election tomorrow?",
               "Which party did you vote for in the most recent general election?")

# Add placeholder columns to reach position 72
n_placeholders_before <- 71 - length(base_cols) - 6 - 1 - 4  # minus opt, env_forced, env_likert
placeholder_names_before <- paste0("Placeholder_", 1:max(1, n_placeholders_before))

# Reorder columns to match expected structure
# Structure: base cols, opt cols (6), placeholders, env cols, more placeholders, insecurity, radical at 72-83

# For simplicity, let's create a new data frame with correct column order
ordered_data <- synth_data[, c("obs_id", "Region", "Gender", "Age...1", "Age...2",
                                "What is your highest educational qualification?",
                                "Which party would you vote for if there were a general election tomorrow?",
                                "Which party did you vote for in the most recent general election?")]

# Add optimism columns
for (col in opt_item_names) {
  ordered_data[[col]] <- synth_data[[col]]
}

# Add environment forced choice
ordered_data[[env_forced_col]] <- synth_data[[env_forced_col]]

# Add environment likert
for (col in env_likert_names) {
  ordered_data[[col]] <- synth_data[[col]]
}

# Add placeholders to reach position 72 for first radical item
n_cols_so_far <- ncol(ordered_data)
n_placeholders_needed <- 71 - n_cols_so_far
if (n_placeholders_needed > 0) {
  for (i in 1:n_placeholders_needed) {
    ordered_data[[paste0("Placeholder_", i)]] <- NA
  }
}

# Add radical items at positions 72-75
for (i in 1:4) {
  ordered_data[[radical_items[[i]]$name]] <- synth_data[[radical_items[[i]]$name]]
}

# Positions 76-77 placeholders
ordered_data$Placeholder_76 <- NA
ordered_data$Placeholder_77 <- NA

# Add radical items at positions 78-79
for (i in 5:6) {
  ordered_data[[radical_items[[i]]$name]] <- synth_data[[radical_items[[i]]$name]]
}

# Positions 80-81 placeholders
ordered_data$Placeholder_80 <- NA
ordered_data$Placeholder_81 <- NA

# Add radical items at positions 82-83
for (i in 7:8) {
  ordered_data[[radical_items[[i]]$name]] <- synth_data[[radical_items[[i]]$name]]
}

# Add insecurity items
for (col in insec_item_names) {
  ordered_data[[col]] <- synth_data[[col]]
}

# ── Save Synthetic Data =======================================================

output_dir <- here("data", "synthetic")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_path <- here(output_dir, "synthetic_survey_data.rds")
saveRDS(ordered_data, output_path)

message("\n=========================================")
message("Synthetic data generated successfully!")
message("N = ", N_SYNTH, " respondents")
message("Saved to: ", output_path)
message("=========================================")

# Print column positions for verification
rad_positions <- which(names(ordered_data) %in% sapply(radical_items, function(x) x$name))
message("\nRadical item column positions: ", paste(rad_positions, collapse = ", "))
