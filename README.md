# Mapping Multidimensional Climate Attitudes in Britain

[![Paper - Substantive](https://img.shields.io/badge/Paper-Substantive_Findings-2ea44f)](paper/main_substantive.pdf)
[![Paper - Technical](https://img.shields.io/badge/Paper-Full_Methodology-0969da)](paper/main_technical.pdf)
![Submitted](https://img.shields.io/badge/Submitted-June_2025-grey)

A Bayesian hierarchical latent trait model for analysing UK public attitudes toward climate policy.

## Overview

This project develops a three-dimensional measurement model capturing:

- **φ (Economic Optimism)**: Confidence in economic prospects
- **θ (Environmentalism)**: Prioritisation of environmental protection over economic growth
- **ψ (Support for Radical Reform)**: Preference for systemic change vs. status quo

The hierarchical structure accounts for variation across:
- Individual demographics (age, gender, education, material insecurity)
- Political party affiliation (10 categories)
- UK region (10 regions)

## Key Findings

- Material insecurity shows strong negative association with economic optimism
- Age effects differ across latent dimensions
- Significant party-level variation, particularly between Green/Labour and Conservative/Reform UK
- Positive residual correlation between environmentalism and economic optimism

## Repository Structure

```
├── R/                          # Analysis scripts
│   ├── 00_generate_synthetic_data.R
│   ├── 01_wrangling.R
│   ├── 02_modelling.R
│   ├── 03_diagnostics.R
│   └── 04_quantities_of_interest.R
├── stan/                       # Stan model
│   └── hierarchical_latent_traits.stan
├── data/                       # Data directory
│   └── synthetic/              # Synthetic data for testing
├── outputs/                    # Generated outputs (gitignored)
│   ├── model/
│   ├── diagnostic_plots/
│   └── qoi_plots/
└── paper/                      # Research note (LaTeX)
    ├── main_substantive.pdf    # Substantive findings
    ├── main_technical.pdf      # Full methodology
    └── figures/
```

## Quick Start

### Prerequisites

- R 4.0+
- CmdStan (for Stan model compilation)
- Required R packages: `tidyverse`, `cmdstanr`, `posterior`, `bayesplot`, `ggdist`, `here`

### Installation

```bash
# Clone repository
git clone https://github.com/henrybaker/climate-attitudes-uk_bayesian-hierarchical-modelling.git
cd climate-attitudes-uk_bayesian-hierarchical-modelling

# Install R dependencies (using renv)
Rscript -e "install.packages('renv'); renv::restore()"
```

### Running the Analysis

```bash
# 1. Generate synthetic data (if real data not available)
Rscript R/00_generate_synthetic_data.R

# 2. Prepare data for Stan
Rscript R/01_wrangling.R

# 3. Fit the model (takes ~30-60 minutes)
Rscript R/02_modelling.R

# 4. Generate diagnostics
Rscript R/03_diagnostics.R

# 5. Extract quantities of interest
Rscript R/04_quantities_of_interest.R
```

## Data Availability

The original survey data is from the [Looking for Growth](https://tracker.lookingforgrowth.uk/) academic partnership and is **proprietary**. It cannot be included in this repository.

For demonstration and code testing, synthetic data can be generated using `R/00_generate_synthetic_data.R`. The synthetic data matches the expected column structure but uses simulated responses based on correlated latent traits.

Researchers interested in accessing the real data should contact the Looking for Growth project directly.

## Research Note

Two versions of the research note are available:

| Version | Description |
|---------|-------------|
| **Substantive** (`paper/main_substantive.pdf`) | Domain-focused findings for social scientists and policy researchers |
| **Technical** (`paper/main_technical.pdf`) | Full methodology including model specification, priors, and diagnostics |

## Model Summary

The model is a 3-dimensional hierarchical latent trait model with:

- **Measurement models**: Factor analysis structure linking observed items to latent traits
- **Structural model**: Latent traits predicted by demographics + random intercepts for region/party
- **Non-centred parameterisation**: For efficient MCMC sampling
- **Generated quantities**: Posterior predictive checks and Bayesian R²

See `stan/hierarchical_latent_traits.stan` for the full model specification.

## Citation

If you use this code or methodology, please cite:

```
Baker, H. (2025). Mapping Multidimensional Climate Attitudes in Britain:
A Bayesian Hierarchical Latent Trait Approach.
GitHub: https://github.com/henrybaker/climate-attitudes-uk_bayesian-hierarchical-modelling
```

## Author

Henry Baker

## Licence

Code: MIT Licence
Paper content: CC-BY-4.0
