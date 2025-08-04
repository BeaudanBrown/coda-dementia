library(targets)
library(autometric)
library(tarchetypes)
library(crew)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")
ncpus <- future::availableCores() - 1
data.table::setDTthreads(1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

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
  covar_sensitivity_targets,
  representative_targets,
  cum_test_targets
)
