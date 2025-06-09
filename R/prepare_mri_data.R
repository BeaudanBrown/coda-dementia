### Prepare MRI data for analysis ###
prepare_mri <- function(df_raw) {
  ## read data

  mri <- fread(file.path(
    data_dir,
    "../../../Generic_data/MRI/mri_full_trimmed.csv"
  ))

  ## Filter to MRI data available

  mri$head_scale <-
    mri$volumetric_scaling_from_t1_head_image_to_standard_space

  mri <- mri |> filter(!is.na(head_scale))

  # create variables of interest

  mri$tot_wmh <- mri$total_vol_white_matter_hyperintensities_from_t1_and_t2_flair_images

  mri$tbv <- (mri$vol_grey_matter_norm + mri$vol_white_matter_norm) / 1000

  mri$wmv <- mri$vol_white_matter_norm / 1000

  mri$gmv <- mri$vol_grey_matter_norm / 1000

  mri$hip <-
    ((mri$vol_hippocampus_l + mri$vol_hippocampus_r) *
      mri$head_scale) /
    1000

  mri$log_wmh <- log(mri$tot_wmh * mri$head_scale)

  ## select relevant vars

  mri <- mri |>
    select(
      eid,
      age_assessment_mri1,
      date_mri1,
      assessment_centre_mri1,
      height_mri1,
      head_scale,
      tot_wmh,
      tbv,
      wmv,
      gmv,
      hip,
      log_wmh
    )

  ## add in QC variables

  qc <- fread(file.path(data_dir, "../../../Generic_data/MRI/mri_qc.csv"))

  qc <- qc |>
    select(eid, ends_with("-2.0")) |>
    as_tibble()

  qc <- qc |>
    rename(
      "discrep_t1_stand_lin" = `25731-2.0`,
      "discrep_t1_stand_nonlin" = `25732-2.0`,
      "warping_t1" = `25733-2.0`,
      "inv_sig_noise_t1" = `25734-2.0`,
      "inv_contr_noise_t1" = `25735-2.0`,
      "discrep_t2_t1" = `25736-2.0`,
      "discrep_swi_t1" = `25738-2.0`,
      "discrep_rfmri_t1" = `25739-2.0`,
      "discrep_tfmri_t1" = `25740-2.0`,
      "mean_tfmri_headmot" = `25742-2.0`,
      "inv_temp_signoise_pp_rfmri" = `25743-2.0`,
      "inv_temp_signoise_art_rfmri" = `25744-2.0`,
      "num_dmri_outslices" = `25746-2.0`,
      "scan_lat_bpos" = `25756-2.0`,
      "scan_trans_bpos" = `25757-2.0`,
      "scan_long_bpos" = `25758-2.0`,
      "scan_tabpos" = `25759-2.0`,
      "acq_prot_phase" = `25780-2.0`
    )

  qc <- qc |> select(-starts_with("2"))

  ## Merge

  mri <- left_join(mri, qc, by = "eid")

  ## add main data

  mri2 <- left_join(df_raw, mri, by = "eid")

  ## MRI data available

  mri2 <- filter(mri2, !is.na(head_scale))

  ## Accelerometry before MRI?

  mri2 <- mri2 |>
    mutate(
      date_mri1 = as_date(date_mri1),
      calendar_date = as_date(calendar_date)
    ) |>
    mutate(mri_accel_time_dif = as.integer(date_mri1 - calendar_date)) |>
    mutate(mri_before_accel = ifelse(mri_accel_time_dif < 0, 1, 0))

  ## Remove participants whose mri was before actigraphy

  mri2 <- mri2 |> filter(mri_before_accel == 0)

  ## select analysis variables

  mri_model_data <- select(
    mri2,
    calendar_date,
    date_mri1,
    eid,
    assessment_centre_mri1,
    tbv,
    wmv,
    gmv,
    hip,
    log_wmh,
    head_scale,
    mean_tfmri_headmot,
    scan_lat_bpos,
    scan_trans_bpos,
    scan_long_bpos,
    scan_tabpos
  )

  ## write it out

  return(mri_model_data)
}

## Add in main dataset predictors

make_mri_df <- function(mri_vars, df) {
  mri_df <- mri_vars

  mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)

  ## Read main data
  mri_df <- left_join(mri_df, df, by = "eid")

  # Age at MRI scan
  mri_df <- mri_df |>
    mutate(
      age_mri = (as.numeric(
        as.Date(date_mri1) - as.Date(calendar_date)
      ) /
        365) +
        age_accel
    )

  return(mri_df)
}
