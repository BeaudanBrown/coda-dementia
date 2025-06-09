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
  tar_target(df_raw, create_data()),
  tar_target(df, prepare_dataset(df_raw)),
  # Run bootstrap for primary model
  tar_target(
    primary,
    run_bootstrap(
      data = df,
      create_formula_fn = get_primary_formula
    )
  ),

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
  tar_target(mri_df, make_mri_df(mri_vars, df)),
  tar_target(full_best_worst, get_best_and_worst_comp(df)),
)
