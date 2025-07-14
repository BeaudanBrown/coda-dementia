library(targets)
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
maxit <- 5 # Number of MICE iterations
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
    "extrafont"
  ),
  format = "qs",
  controller = crew_controller_local(
    options_local = crew_options_local(log_directory = "./logs"),
    workers = ncpus
  ),
  seed = 5678
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
  tar_target(test_df, sample_frac(df, 0.1)),
  tar_target(imp, impute_data(test_df, m, maxit)),
  tar_target(imp2, back_to_factor(imp), pattern = map(imp)),
  tar_target(
    timegroup_cuts,
    make_cuts(df)
  ),
  tar_target(
    imp_wide,
    widen_data(imp2, timegroup_cuts),
    pattern = map(imp2)
  ),
  tar_target(
    sub_names,
    c("avg_sleep", "avg_mvpa", "avg_light", "avg_inactivity")
  ),
  tar_target(sub_durations, seq(from = 15, to = 60, by = 15)),
  tar_target(
    substitutions,
    {
      expand.grid(
        from_var = sub_names,
        to_var = sub_names,
        stringsAsFactors = FALSE
      ) |>
        dplyr::filter(!from_var == to_var)
    }
  ),
  tar_target(
    substitution_results,
    apply_substitution(
      imp_wide,
      substitutions$from_var,
      substitutions$to_var,
      sub_durations,
      timegroup_cuts
    ),
    pattern = cross(substitutions, sub_durations),
    iteration = "vector"
  ),
  tar_target(
    sub_test,
    process_substitution(
      imp_wide,
      substitution_results,
      baseline = c(
        "fruit_veg",
        "alc_freq",
        "sex",
        "retired",
        "shift",
        "apoe_e4",
        "highest_qual",
        "townsend_deprivation_index",
        "psych_meds",
        "ethnicity",
        "avg_total_household_income",
        "smok_status"
      )
    ),
    pattern = map(substitution_results),
    iteration = "list"
  )
)
