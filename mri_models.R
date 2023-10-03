mri_df <-
  fread(file.path(data_dir, "MRI/mri_full_trimmed_v3.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
mri_df <- mri_df |> rename("insomnia_scale_sr" = "insomnia_sr")
mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)
mri_df <- mri_df |> select(-starts_with("dur_day_total_"))

# MRI QC data
#########################################################################################
mri_qc_df <-
  fread(file.path(data_dir, "MRI/mri_qc.csv"))

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

# Accelerometery data
accel_df <-
  fread(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_personsumMM_output.csv"))

# Sleep disorders data
sleep_dis_df <- fread(file.path(data_dir, "Sleep_disorders/sleep_disorders_selfreport_primarycare_hosp.csv"))

mri_df <- left_join(mri_df, accel_df, by = "eid")
mri_df <- left_join(mri_df, mri_qc_df, by = "eid")
mri_df <- left_join(mri_df, sleep_dis_df, by = "eid")

mri_df$OSA_dx <- ifelse(mri_df$OSA_sr == 1 | mri_df$date_osa_dx <= mri_df$accel_date, 1, 0)
mri_df$insomnia_dx <-
  ifelse(mri_df$insomnia_sr == 1 | mri_df$date_insomnia_dx <= mri_df$accel_date, 1, 0)
mri_df$sleep_disorder_dx <-
  ifelse(mri_df$sleep_disorder_sr == 1 | mri_df$date_any_sleep_dx <= mri_df$accel_date, 1, 0)

# replace missing with zero
mri_df <- mri_df |> mutate(across(c(OSA_dx, insomnia_dx, sleep_disorder_dx), ~ if_else(is.na(.), 0, .)))

# set pack years to zero for non-smokers
mri_df$smok_pckyrs <- ifelse(mri_df$smok_status == 0, 0, mri_df$smok_pckyrs)

# set reference categories for factors
mri_df$diagnosed_diabetes <- as.factor(mri_df$diagnosed_diabetes)
levels(mri_df$diagnosed_diabetes) <- c("prefer not answer", "dont know", "no", "yes")
mri_df$smok_status <- as.factor(mri_df$smok_status)
levels(mri_df$smok_status) <- c("prefer not answer", "never", "former", "current")
mri_df$apoe_e4 <- as.factor(mri_df$apoe_e4)
mri_df$highest_qual <- fct_relevel(mri_df$highest_qual, "Grad")
mri_df$ethnicity <- fct_relevel(mri_df$ethnicity, "white")
mri_df$insomnia_scale_sr <- fct_relevel(mri_df$insomnia_scale_sr, "never")
mri_df$chronotype <- fct_relevel(mri_df$chronotype, "morning")
mri_df$freq_depressed_twoweeks <- fct_relevel(mri_df$freq_depressed_twoweeks, "not at all")
mri_df$avg_total_household_income <- fct_relevel(mri_df$avg_total_household_income, "31-50")
mri_df$sick_disabled <- ifelse(mri_df$employment == "sick or disabled", 1, 0)
mri_df$retired <- ifelse(mri_df$employment == "retired", 1, 0)
mri_df$overall_health_rating <- as.factor(mri_df$overall_health_rating)
levels(mri_df$overall_health_rating) <- c("prefer not answer", "dont know", "excellent", "good", "fair", "poor")
mri_df$shift <- ifelse(mri_df$job_night_shift %in% c(3, 4) | mri_df$job_shift_work %in% c(3, 4), 1, 0)

# mark prefer not answer as missing
mri_df <- mri_df |>
  mutate(across(where(is.factor), as.character)) |>
  mutate(across(where(is.character), ~ na_if(., "prefer not answer"))) |>
  mutate(across(where(is.character), ~ na_if(., "dont know"))) |>
  mutate(across(where(is.character), ~ na_if(., ""))) |>
  mutate(across(where(is.character), as.factor)) |>
  mutate(across(where(is.factor), fct_drop))

# TODO: Consult data dictionary to figure out why totals we are calculating don't match existing totals
mri_df$awake_sleep <-
  mri_df$dur_spt_wake_IN_min_pla +
  mri_df$dur_spt_wake_LIG_min_pla +
  mri_df$dur_spt_wake_MOD_min_pla +
  mri_df$dur_spt_wake_VIG_min_pla

# TODO: Deal with these columns being duplicated from join
mri_df$mins_worn <-
  mri_df$dur_spt_sleep_min_pla +
  mri_df$dur_day_total_IN_min_pla +
  mri_df$dur_day_total_LIG_min_pla +
  mri_df$dur_day_total_MOD_min_pla +
  mri_df$dur_day_total_VIG_min_pla +
  mri_df$awake_sleep

# Variables normalised to 1440 relative to their proportion of total wear time
mri_df$sleep_n <-
  (mri_df$dur_spt_sleep_min_pla / mri_df$mins_worn) * mins_in_day
mri_df$inactive_n <-
  ((mri_df$dur_day_total_IN_min_pla + mri_df$dur_spt_wake_IN_min_pla) / mri_df$mins_worn) * mins_in_day
mri_df$light_n <-
  ((mri_df$dur_day_total_LIG_min_pla + mri_df$dur_spt_wake_LIG_min_pla) / mri_df$mins_worn) * mins_in_day
mri_df$moderate_n <-
  ((mri_df$dur_day_total_MOD_min_pla + mri_df$dur_spt_wake_MOD_min_pla) / mri_df$mins_worn) * mins_in_day
mri_df$vigorous_n <-
  ((mri_df$dur_day_total_VIG_min_pla + mri_df$dur_spt_wake_VIG_min_pla) / mri_df$mins_worn) * mins_in_day

mri_df$mvpa_n <- mri_df$moderate_n + mri_df$vigorous_n

# Update age variable to age at accelerometry study
mri_df$age_accel <-
  mri_df$age_assessment + ((as.Date(mri_df$accel_date) - as.Date(mri_df$date_baseline)) / 365)
mri_df$age_accel <-
  as.numeric(mri_df$age_accel)

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

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)
v <- gsi.buildilrBase(t(sbp))

mri_comp <- acomp(data.frame(mri_df$sleep_n, mri_df$inactive_n, mri_df$light_n, mri_df$mvpa_n))
mri_base_ilr <-
  ilr(mri_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

# Model data for mri model
mri_model_data <- select(
  mri_df,
  assessment_centre_mri1,
  tbv, wmv, gmv, hip, log_wmh,
  mean_tfmri_headmot,
  scan_lat_bpos,
  scan_trans_bpos,
  scan_long_bpos,
  scan_tabpos,
  sleep_n, inactive_n, light_n, mvpa_n,
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
  overall_health_rating,
  OSA_dx,
  insomnia_dx,
  sleep_disorder_dx
) |>
  cbind(mri_base_ilr)

options(datadist = datadist(mri_model_data))

tbv_model <- ols(
  tbv ~
    rcs(R1, 3) + rcs(R2, 3) + rcs(R3, 3) +
    rcs(avg_sri, c(47, 61, 72)) +
    rcs(age_mri, c(54, 65, 75)) +
    retired +
    rcs(townsend_deprivation_index, c(-5, -2, 2.5)) +
    sex +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    highest_qual +
    apoe_e4 +
    smok_status +
    rcs(mean_tfmri_headmot, c(0.08, 0.13, 0.22)) +
    rcs(scan_lat_bpos, c(-2.8, 0.6, 4.3)) +
    rcs(scan_trans_bpos, c(59, 66, 74)) +
    rcs(scan_long_bpos, c(-53.6, -34, 16.7)) +
    scan_tabpos +
    assessment_centre_mri1,
  data = mri_model_data
)

