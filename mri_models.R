source("ideal_comp.R")

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)
v <- gsi.buildilrBase(t(sbp))

make_ilr <- function(comp) {
  return(ilr(acomp(comp), V = v))
}
best_and_worst <- get_best_and_worst_comp()
best_and_worst <- as.data.frame(apply(best_and_worst, 2, make_ilr))

mri_df <-
  fread(file.path(data_dir, "../MRI/mri_full_trimmed_v3.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
mri_df <- mri_df |> rename("insomnia_scale_sr" = "insomnia_sr")
mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)
mri_df <- mri_df |> select(-starts_with("dur_day_total_"))

# MRI QC data
#########################################################################################
mri_qc_df <-
  fread(file.path(data_dir, "../MRI/mri_qc.csv"))

mri_qc_df <- mri_qc_df |>
  select(eid, ends_with("-2.0")) |>
  as_tibble()

mri_qc_df <- mri_qc_df |>
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

mri_qc_df <- mri_qc_df |> select(-starts_with("2"))

boot_df <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

mri_df <- mri_df |>
  select(-any_of(names(boot_df)[names(boot_df) != "eid"])) |>
  left_join(boot_df, by = "eid")

mri_df <- mri_df |>
  left_join(mri_qc_df |> select(-any_of(names(mri_df)[names(boot_df) != "eid"])), by = "eid")

# Age at MRI scan
mri_df <- mri_df %>%
  mutate(
    date_mri1 = parse_date_time(date_mri1, "ymd"),
    date_birth = paste(year_birth, month_birth, sep = "-")
  ) %>%
  mutate(date_birth = parse_date_time(date_birth, "ym")) %>%
  mutate(age_mri = as.numeric(date_mri1 - date_birth) / 365)

## Create derived variables
mri_df$tbv <-
  (mri_df$vol_grey_matter_norm + mri_df$vol_white_matter_norm) / 1000
mri_df$wmv <-
  mri_df$vol_white_matter_norm / 1000
mri_df$gmv <-
  mri_df$vol_grey_matter_norm / 1000
mri_df$hip <-
  (mri_df$vol_hippocampus_l + mri_df$vol_hippocampus_r) *
  mri_df$volumetric_scaling_from_t1_head_image_to_standard_space / 1000
mri_df$tot_wmh <- mri_df$total_vol_white_matter_hyperintensities_from_t1_and_t2_flair_images
mri_df$log_wmh <-
  log(mri_df$tot_wmh * mri_df$volumetric_scaling_from_t1_head_image_to_standard_space)

# Model data for mri model
mri_model_data <- select(
  mri_df,
  R1, R2, R3,
  avg_sleep, avg_inactivity, avg_light, avg_mvpa,
  assessment_centre_mri1,
  tbv, wmv, gmv, hip, log_wmh,
  mean_tfmri_headmot,
  scan_lat_bpos,
  scan_trans_bpos,
  scan_long_bpos,
  scan_tabpos,
  smok_status,
  dem,
  time_to_dem,
  bp_syst_avg,
  avg_sri,
  age_assessment,
  bp_med,
  any_cvd,
  highest_qual,
  apoe_e4,
  BMI,
  antidepressant_med,
  antipsychotic_med,
  avg_WASO,
  retired,
  shift,
  pa_modvig,
  insomnia_med,
  ethnicity,
  age_mri,
  sex,
  parent_dementa,
  avg_total_household_income,
  townsend_deprivation_index,
  vol_white_matter_norm,
  vol_grey_matter_norm,
  vol_hippocampus_l,
  vol_hippocampus_r,
  volumetric_scaling_from_t1_head_image_to_standard_space,
  mri_accel_time_dif,
  tot_wmh,
  dem,
  avg_sleep_duration,
  time_to_dem,
  overall_health_rating
)

options(datadist = datadist(mri_model_data))

mri_model_data <- na.omit(mri_model_data)

predict_mri_outcome <- function(outcome_var, model_data, best_and_worst) {
  knots_avg_sri_str <-
    paste(quantile(model_data[["avg_sri"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_age_mri_str <-
    paste(quantile(model_data[["age_mri"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_deprivation_str <-
    paste(quantile(model_data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_tfmri_str <-
    paste(quantile(model_data[["mean_tfmri_headmot"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_lat_bpos_str <-
    paste(quantile(model_data[["scan_lat_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_trans_bpos_str <-
    paste(quantile(model_data[["scan_trans_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_long_bpos_str <-
    paste(quantile(model_data[["scan_long_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")

  model_formula <- as.formula(paste(outcome_var, " ~
      rcs(R1, 3) + rcs(R2, 3) + rcs(R3, 3) +
      rcs(avg_sri, c(", knots_avg_sri_str, ")) +
      rcs(age_mri, c(", knots_age_mri_str, ")) +
      retired +
      rcs(townsend_deprivation_index, c(", knots_deprivation_str, ")) +
      sex +
      antidepressant_med +
      antipsychotic_med +
      insomnia_med +
      ethnicity +
      avg_total_household_income +
      highest_qual +
      apoe_e4 +
      smok_status +
      rcs(mean_tfmri_headmot, c(", knots_tfmri_str, ")) +
      rcs(scan_lat_bpos, c(", knots_lat_bpos_str, ")) +
      rcs(scan_trans_bpos, c(", knots_trans_bpos_str, ")) +
      rcs(scan_long_bpos, c(", knots_long_bpos_str, ")) +
      scan_tabpos +
      assessment_centre_mri1"))

  tbv_model <- ols(model_formula,
    data = model_data
  )

  ref_row <- as.data.frame(model_data[1, ])
  ref_row[c("R1", "R2", "R3")] <- best_and_worst$best
  best <- predict(tbv_model, newdata = ref_row)

  ref_row[c("R1", "R2", "R3")] <- best_and_worst$worst
  worst <- predict(tbv_model, newdata = ref_row)

  ref_row[c("R1", "R2", "R3")] <- best_and_worst$most_common
  common <- predict(tbv_model, newdata = ref_row)

  return(data.frame(
    reference = c("worst", "common", "best"),
    value = c(worst, common, best)
  ))
}

results <- lapply(c("tbv", "wmv", "gmv", "tot_wmh", "log_wmh"), function(outcome) {
  df <- predict_mri_outcome(outcome, mri_model_data, best_and_worst)
  df$outcome <- outcome
  return(df)
})

final_df <- do.call(rbind, results)

ggplot(final_df, aes(x = reference, y = value, color = outcome)) +
  geom_point() +
  facet_wrap(~outcome, scales = "free") +
  xlab("") +
  ylab("Value") +
  theme_bw()
