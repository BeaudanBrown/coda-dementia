library(targets)
library(autometric)
library(tarchetypes)
library(crew)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")
ncpus <- future::availableCores() - 1

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
mins_in_day <- 1440
mins_in_hour <- 60
sub_steps <- 4
sub_step_mins <- mins_in_hour / sub_steps
m <- 1 # Number of imputed datasets
maxit <- 5 # Number of MICE iterations
n_boots <- 3

# set target configs
tar_config_set(store = cache_dir)
unlink("./logs/*", recursive = FALSE)
unlink("./log.txt", recursive = FALSE)

controller <- crew_controller_local(
  options_local = crew_options_local(log_directory = "./logs"),
  options_metrics = crew_options_metrics(
    path = "/dev/stdout",
    seconds_interval = 10
  ),
  workers = ncpus
)

if (tar_active()) {
  controller$start()
  log_start(
    path = "log.txt",
    seconds = 30,
    pids = controller$pids()
  )
}

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
    "RhpcBLASctl"
  ),
  format = "qs",
  controller = controller,
  seed = 5678
)

# Run the R scripts in the R/ folder
tar_source()

# Set data table cores to 1

data.table::setDTthreads(1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

## pipeline

list(
  #### FILES ####
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

  #### PREPARE DATA ####
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

  #### DEFINE ANALYSIS PARAMETERS ####
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
        dplyr::filter(
          !from_var == to_var &
            (from_var == "avg_sleep" | to_var == "avg_sleep")
        )
    }
  ),
  tar_target(timegroup_cuts, make_cuts(df)),
  tar_target(final_time, length(timegroup_cuts) - 1),

  #### IMPUTATION ####
  tar_rep(boots, bootstrap_sample(df), batches = n_boots),
  tar_target(imp, impute_data(boots, m, maxit), pattern = map(boots)),

  #### PRIMARY ANALYSIS ####
  tar_target(
    primary_models,
    fit_models(imp, timegroup_cuts, get_primary_formula),
    pattern = map(imp)
  ),
  tar_target(
    primary_ref_risk,
    get_ref_risk(imp, primary_models, final_time),
    pattern = map(imp, primary_models),
    iteration = "list"
  ),
  tar_target(
    primary_sub_risk,
    get_sub_risk(
      imp,
      substitutions$from_var,
      substitutions$to_var,
      sub_durations,
      primary_models,
      final_time
    ),
    pattern = cross(
      map(imp, primary_models),
      cross(substitutions, sub_durations)
    ),
    iteration = "list"
  ),
  tar_target(primary_results, intervals(primary_ref_risk, primary_sub_risk)),
  tar_target(
    test_sub,
    get_sub_risk(
      imp,
      "avg_mvpa",
      "avg_sleep",
      60,
      primary_models,
      final_time
    ),
    map(imp, primary_models)
  )
)
