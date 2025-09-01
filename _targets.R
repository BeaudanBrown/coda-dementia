library(targets)
library(autometric)
library(tarchetypes)
library(crew)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")
ncpus <- future::availableCores() - 1
# ncpus <- 2

Sys.setenv(R_DATATABLE_NUM_THREADS = 1)
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)

# set target configs
tar_config_set(store = cache_dir)
unlink("./logs/*", recursive = FALSE)

controller <- crew_controller_local(
  options_local = crew_options_local(log_directory = "./logs"),
  options_metrics = crew_options_metrics(
    path = "/dev/stdout",
    seconds_interval = 10
  ),
  workers = ncpus
)

# Set target options:
tar_option_set(
  packages = c(
    "utils",
    "mice",
    "tidyverse",
    "survival",
    "rms",
    "compositions",
    "parallel",
    "boot",
    "rlang",
    "data.table",
    "cowplot",
    "extrafont",
    "fastglm",
    "dotenv",
    "lme4",
    "emmeans",
    "lubridate",
    "mvtnorm",
    "lmtp",
    "extrafont",
    "RhpcBLASctl",
    "ks",
    "grid",
    "patchwork"
  ),
  format = "qs",
  controller = controller,
  workspace_on_error = TRUE,
  seed = 5678
)

# Run the R scripts in the R/ folder
tar_source()

## pipeline

source("R/constants.R")
source("data_targets.R")
source("primary_targets.R")
source("mri_targets.R")
source("cum_targets.R")
source("cum_test.R")
source("sensitivity_reverse_causation_targets.R")
source("sensitivity_reverse_causation_test.R")
source("sensitivity_covar_targets.R")
source("sensitivity_representative_targets.R")

list(
  #### DATA PREPARATION ####
  data_targets,
  #### PRIMARY ANALYSIS ####
  primary_targets,
  #### MRI ANALYSIS ####
  mri_targets,
  #### CUMULATIVE ANALYSIS ####
  cum_targets,
  #### SENSITIVITIES ####
  reverse_causation_targets,
  reverse_causation_targets_test,
  covar_sensitivity_targets,
  representative_targets,

  tar_target(
    sub_tables,
    list(
      shorts = process_sub_table(primary_risk_ratios_short_sleeper),
      avgs = process_sub_table(primary_risk_ratios_avg_sleeper),
      longs = process_sub_table(primary_risk_ratios_long_sleeper)
    )
  )
)
