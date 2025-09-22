# Optimizing Dementia Risk Reduction

This repository contains the analysis code for the paper Sleep and Physical Activity Trade-offs and Dementia Risk: A Prospective Cohort Study in UK Biobank Participants by Yiallourou et al. This project made use of the UK Biobank resource. Data analysis was performed by Beaudan Campbell-Brown and Lachlan Cribb.

## Getting Started

### Prerequisites

This project uses [Nix](https://nixos.org/) for reproducible dependency management. To get started:

1. Install Nix
2. Enable flakes (if not already enabled)
3. Run `nix develop` to enter the development environment

All R packages and system dependencies are managed through the `flake.nix` configuration, ensuring a reproducible analysis environment.

### Running the Analysis

The analysis pipeline is implemented using the R [targets](https://books.ropensci.org/targets/) package for reproducible workflow management. To execute the full pipeline:

```r
# Load the targets library
library(targets)

# Run the entire pipeline
tar_make()

# View the pipeline structure
tar_visnetwork(targets_only = TRUE)
```

## Project Structure

The pipeline is organized into several target files:
- `data_targets.R` - Data preparation targets
- `primary_targets.R` - Primary analysis targets
- `mri_targets.R` - MRI analysis targets
- `cum_targets.R` - Cumulative incidence (risk) plot targets
- `sensitivity_*.R` - Various sensitivity analysis targets

## Dependencies

This project uses Nix flakes for dependency management, providing:
- R environment with all required packages
- System-level dependencies (Quarto, formatters)
- Reproducible package versions across different systems