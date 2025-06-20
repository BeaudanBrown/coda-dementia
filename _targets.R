library(targets)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")
seed_val <- Sys.getenv("SEED")

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
mins_in_day <- 1440
mins_in_hour <- 60
sub_steps <- 4
sub_step_mins <- mins_in_hour / sub_steps
ncpus <- as.integer(Sys.getenv("NCPUS"))
bootstrap_iterations <- as.integer(Sys.getenv("BOOT_ITRS"))
maxit <- as.integer(Sys.getenv("MAXIT"))

# set target configs
tar_config_set(store = cache_dir)

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
    "renv",
    "extrafont"
  ),
  format = "qs"
)

# Run the R scripts in the R/ folder
tar_source()

# Set data table cores to 1

data.table::setDTthreads(1)

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

  # Run bootstrap for sensitivity analysis model 1
  tar_target(
    s1,
    run_bootstrap(
      data = df,
      create_formula_fn = get_s1_formula
    )
  ),
  # Run bootstrap for sensitivity analysis model 2
  tar_target(
    s2,
    run_bootstrap(
      data = df,
      create_formula_fn = get_s2_formula,
      intervals = intervals
    )
  ),
  ## Run bootstrap for sensitivity analysis model 3
  # standardising to pseudo pop of Schoeler et al.
  tar_target(
    s3,
    run_bootstrap(
      data = df,
      create_formula_fn = get_s3_formula,
      intervals = intervals,
      empirical = FALSE
    )
  ),
  tar_target(main_plots, produce_plots(primary, s1, s2, s3)),
  tar_target(mri_vars, prepare_mri(df_raw)),
  tar_target(mri_df, make_mri_df(df, mri_vars)),
  tar_target(
    mri_boot,
    run_mri_bootstrap(mri_df, full_best_worst, get_mri_formula)
  ),
  tar_target(mri_output, process_mri_output(mri_boot)),
  tar_target(full_best_worst, get_best_and_worst_comp(df))
)
