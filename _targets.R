library(targets)
library(tarchetypes)
library(crew)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")
ncpus <- as.integer(Sys.getenv("NCPUS"))

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
mins_in_day <- 1440
mins_in_hour <- 60
sub_steps <- 4
sub_step_mins <- mins_in_hour / sub_steps
m <- 10 # Number of imputed datasets
maxit <- 2 # Number of MICE iterations
bootstrap_iterations <- 16 # For ideal/worst plots

# set target configs
tar_config_set(store = cache_dir)
unlink("./logs/*", recursive = FALSE)

# Set target options:
tar_option_set(
  packages = c(
    "utils",
    "mice",
    "tidyverse",
    "survival",
    "rms",
    "compositions",
    "data.table",
    "parallel",
    "boot",
    "rlang",
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
    "targets",
    "tarchetypes",
    "RhpcBLASctl"
  ),
  format = "qs",
  controller = crew_controller_local(
    options_local = crew_options_local(log_directory = "./logs"),
    workers = ncpus
  ),
  seed = 5678
)

source("sub_targets.R")

# Run the R scripts in the R/ folder
tar_source()

# Set data table cores to 1

data.table::setDTthreads(1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

## pipeline
list(
  tar_target(
    core_file,
    file.path(data_dir, Sys.getenv("CORE_FILE")),
    format = "file"
  ),
  tar_target(
    demdeath_file,
    file.path(data_dir, Sys.getenv("DEMDEATH_FILE")),
    format = "file"
  ),
  tar_target(
    snp_file,
    file.path(data_dir, Sys.getenv("SNP_FILE")),
    format = "file"
  ),
  tar_target(
    diet_file,
    file.path(data_dir, Sys.getenv("DIET_FILE")),
    format = "file"
  ),
  tar_target(
    accel_file,
    file.path(data_dir, Sys.getenv("ACCEL_FILE")),
    format = "file"
  ),
  tar_target(
    sleep_file,
    file.path(data_dir, Sys.getenv("SLEEP_FILE")),
    format = "file"
  ),
  tar_target(
    sri_file,
    file.path(data_dir, Sys.getenv("SRI_FILE")),
    format = "file"
  ),
  tar_target(
    mri_file,
    file.path(data_dir, Sys.getenv("MRI_FILE")),
    format = "file"
  ),
  tar_target(
    df_raw,
    create_data(
      core_file,
      demdeath_file,
      snp_file,
      diet_file,
      accel_file,
      sleep_file,
      sri_file
    )
  ),
  tar_target(
    disease_file,
    file.path(data_dir, Sys.getenv("DISEASE_FILE")),
    format = "file"
  ),
  tar_target(df, prepare_dataset(df_raw, disease_file)),
  tar_target(test_df, sample_frac(df, 0.03)),
  tar_target(imp, impute_data(df, m, maxit), iteration = "list"),
  tar_target(
    timegroup_cuts,
    make_cuts(df)
  ),
  tar_target(
    imp_wide,
    widen_data(imp, timegroup_cuts),
    iteration = "list"
  ),
  sub_targets
)
