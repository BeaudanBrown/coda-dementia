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
long_sleep_min_hours <- 8
long_sleep_hours <- 9
hrs_in_day <- 24
mins_in_day <- 1440
mins_in_hour <- 60
sub_steps <- 4
sub_step_mins <- mins_in_hour / sub_steps
m <- 1 # Number of imputed datasets
maxit <- 10 # Number of MICE iterations
n_boots <- 500

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
  seed = 5678
)

# Run the R scripts in the R/ folder
tar_source()

# Set data table cores to 1

data.table::setDTthreads(1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

short_sleeper_filter_fn <- function(df) {
  df |> filter(avg_sleep < short_sleep_hours * 60)
}

avg_sleeper_filter_fn <- function(df) {
  df |>
    filter(
      avg_sleep >= short_sleep_hours * 60 & avg_sleep <= long_sleep_hours * 60
    )
}

long_sleeper_filter_fn <- function(df) {
  df |> filter(avg_sleep > long_sleep_min_hours * 60)
}

substitutions <- expand.grid(
  from_var = c("avg_mvpa", "avg_light", "avg_inactivity"),
  to_var = "avg_sleep",
  stringsAsFactors = FALSE
)
durations <- data.frame(
  duration = seq(from = -60, to = 60, by = 15)[
    seq(from = -60, to = 60, by = 15) != 0
  ]
)
all_subs <- merge(substitutions, durations)

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
  tar_target(
    mri_qc_file,
    file.path(data_dir, Sys.getenv("MRI_QC_FILE")),
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
    pattern = map(imp, primary_models)
  ),
  tar_target(
    primary_sub_risk,
    bind_rows(apply(all_subs, 1, function(sub) {
      get_sub_risk(
        imp,
        sub["from_var"],
        sub["to_var"],
        as.numeric(sub["duration"]),
        primary_models,
        final_time
      )
    })),
    pattern = map(imp, primary_models),
  ),
  tar_map(
    values = list(
      filter_fn = rlang::syms(c(
        "long_sleeper_filter_fn",
        "avg_sleeper_filter_fn",
        "short_sleeper_filter_fn"
      ))
    ),
    tar_target(
      primary_ref_avg_risks,
      average_risks(primary_ref_risk, df, filter_fn)
    ),
    tar_target(
      primary_sub_avg_risks,
      average_risks(primary_sub_risk, df, filter_fn),
      pattern = map(primary_sub_risk)
    ),
    tar_target(
      primary_risk_ratios,
      merge_risks(primary_sub_avg_risks, primary_ref_avg_risks)
    ),
    tar_map(
      values = list(
        name = c("inactivity", "light_activity", "mvpa"),
        from = c("avg_inactivity", "avg_light", "avg_mvpa"),
        colour = c("#ff747b", "#6ed853", "#708ff9")
      ),
      names = name,
      tar_target(
        primary_plots,
        make_plot(primary_risk_ratios, from, colour)
      )
    )
  ),
  tar_target(
    all_primary_plots,
    list(
      short_inactive = primary_plots_inactivity_short_sleeper_filter_fn,
      short_light = primary_plots_light_activity_short_sleeper_filter_fn,
      short_mvpa = primary_plots_mvpa_short_sleeper_filter_fn,
      avg_inactive = primary_plots_inactivity_avg_sleeper_filter_fn,
      avg_light = primary_plots_light_activity_avg_sleeper_filter_fn,
      avg_mvpa = primary_plots_mvpa_avg_sleeper_filter_fn,
      long_inactive = primary_plots_inactivity_long_sleeper_filter_fn,
      long_light = primary_plots_light_activity_long_sleeper_filter_fn,
      long_mvpa = primary_plots_mvpa_long_sleeper_filter_fn
    )
  ),
  tar_target(
    primary_grid,
    {
      make_plot_grid(all_primary_plots)
    }
  ),

  #### FIND IDEAL/WORST COMPOSITION ####

  tar_target(synth_comps, generate_compositions(df)),
  tar_target(synth_comps_dens, add_density(df, synth_comps, 0.05)),
  tar_target(
    synth_comps_filtered,
    synth_comps_dens[dens > dens_threshold & !is.na(dens), ]
  ),
  tar_target(
    train_indices,
    sample(seq_len(nrow(df)), size = floor(0.5 * nrow(df)))
  ),
  tar_target(df_train, df[train_indices]),
  tar_target(df_test, df[-train_indices]),
  tar_target(imp_train, impute_data(df_train, m, maxit)),
  tar_target(
    train_models,
    train_model(imp_train, timegroup_cuts, get_primary_formula)
  ),
  tar_target(covar_data, get_covar_data(imp_train)),
  tar_target(
    synth_split,
    split(
      synth_comps_filtered[1:100],
      seq_len(nrow(synth_comps_filtered[1:100]))
    ),
    iteration = "list"
  ),
  tar_target(
    synth_comp_risk,
    get_synth_risk(
      covar_data,
      synth_comps_filtered,
      train_models,
      final_time
    ),
    pattern = map(synth_comps_filtered)
  ),
  tar_target(
    reference_comps,
    {
      list(
        best = synth_comp_risk[order(risk), ][1],
        worst = synth_comp_risk[order(-dens), ][1],
        typical = synth_comp_risk[order(-risk), ][1]
      )
    }
  ),

  #### ESTIMATE RISK FOR IDEAL/WORST COMPS ####

  #### MRI ANALYSIS ####
  tar_target(mri_raw, prepare_mri(df_raw, mri_file, mri_qc_file)),
  tar_target(mri_df, make_mri_df(mri_raw, df)),
  tar_rep(mri_boots, bootstrap_sample(mri_df), batches = n_boots),
  tar_target(
    mri_imp,
    impute_mri_data(mri_boots, m, maxit),
    pattern = map(mri_boots)
  ),
  tar_target(mri_outcomes, c("tbv", "wmv", "gmv", "hip", "log_wmh")),
  tar_map(
    values = list(
      outcome = c("tbv", "wmv", "gmv", "hip", "log_wmh")
    ),
    names = outcome,
    tar_target(
      mri_models,
      get_mri_model(mri_imp, outcome),
      pattern = map(mri_imp),
      iteration = "list"
    ),
    tar_target(
      mri_ref_results,
      get_mri_ref(mri_imp, outcome, mri_models),
      pattern = map(mri_imp, mri_models)
    ),
    tar_target(
      mri_sub_results,
      bind_rows(apply(all_subs, 1, function(sub) {
        get_mri_subs(
          mri_imp,
          outcome,
          mri_models,
          sub["from_var"],
          sub["to_var"],
          as.numeric(sub["duration"])
        )
      })),
      pattern = map(mri_imp, mri_models)
    ),
    tar_map(
      values = list(
        filter_fn = rlang::syms(c(
          "avg_sleeper_filter_fn",
          "short_sleeper_filter_fn"
        ))
      ),
      tar_target(
        mri_ref_avg_estimate,
        average_estimates(mri_ref_results, df, filter_fn)
      ),
      tar_target(
        mri_sub_avg_estimate,
        average_estimates(mri_sub_results, df, filter_fn),
        pattern = map(mri_sub_results)
      ),
      tar_target(
        mri_mean_diffs,
        merge_estimates(mri_sub_avg_estimate, mri_ref_avg_estimate)
      )
    )
  )
  # tar_map(
  #   values = list(
  #     wmh = list(
  #       ylabel = "",
  #       ylabel = "",
  #       ylabel = "",
  #     ),
  #   ),
  #   tar_target(
  #     plot,
  #     make_mri_plot(mri_results, ylabel, sub_name, colour, mri_sub_results)
  #   )
  # )
)
