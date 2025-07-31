data_targets <- list(
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
  tar_target(
    disease_file,
    file.path(data_dir, Sys.getenv("DISEASE_FILE")),
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
  tar_target(df, prepare_dataset(df_raw, disease_file)),

  #### DEFINE ANALYSIS PARAMETERS ####
  tar_target(timegroup_cuts, make_cuts(df)),
  tar_target(final_time, length(timegroup_cuts) - 1),
  tar_target(comp_limits, {
    sleep_quants <- quantile(df$avg_sleep, probs = c(0.01, 0.99), na.rm = TRUE)
    inactivity_quants <- quantile(
      df$avg_inactivity,
      probs = c(0.01, 0.99),
      na.rm = TRUE
    )
    light_quants <- quantile(df$avg_light, probs = c(0.01, 0.99), na.rm = TRUE)
    mvpa_quants <- quantile(df$avg_mvpa, probs = c(0.01, 0.99), na.rm = TRUE)
    list(
      avg_sleep = list(
        lower = sleep_quants[1],
        upper = sleep_quants[2]
      ),
      avg_inactivity = list(
        lower = inactivity_quants[1],
        upper = inactivity_quants[2]
      ),
      avg_light = list(
        lower = light_quants[1],
        upper = light_quants[2]
      ),
      avg_mvpa = list(
        lower = mvpa_quants[1],
        upper = mvpa_quants[2]
      )
    )
  }),

  #### IMPUTATION ####
  tar_rep(boots, bootstrap_sample(df), batches = n_boots),
  tar_target(imp, impute_data(boots, m, maxit), pattern = map(boots)),
  tar_target(full_imp, {
    full_imp <- impute_data(df, m, maxit)[, id := .I]
    full_imp[, tar_batch := 1]
  }),

  #### MRI ANALYSIS ####
  tar_target(mri_raw, prepare_mri(df_raw, mri_file, mri_qc_file)),
  tar_target(mri_df, make_mri_df(mri_raw, df)),
  tar_rep(mri_boots, bootstrap_sample(mri_df), batches = n_boots),
  tar_target(
    mri_imp,
    impute_mri_data(mri_boots, m, maxit),
    pattern = map(mri_boots)
  )
)
