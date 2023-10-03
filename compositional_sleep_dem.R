# Basic Linear Regression and isotemporal substitution with Compositional Explanatory Variables

# Packages and installs
library(tidyverse)
library(compositions)
library(data.table)
library(here)
library(dotenv)
library(rms)
library(survival)
library(ggplot2)
library(mvtnorm)
library(marginaleffects)

sleep_idx <- 1
inactive_idx <- 2
light_idx <- 3
vig_idx <- 4
sub_amount <- 120
mins_in_day <- 1440
hrs_in_day <- 24
short_sleep_hours <- 6
long_sleep_hours <- 9

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

# MRI and sri data
mri_df <-
  fread(file.path(data_dir, "MRI/mri_full_trimmed_v3.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
mri_df <- mri_df |> rename("insomnia_scale_sr" = "insomnia_sr")
mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)
mri_df <- mri_df |> select(-starts_with("dur_day_total_"))

sri_df <-
  fread(file.path(data_dir, "SRI/sri_data_may2023.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
sri_df <- sri_df |> rename("insomnia_scale_sr" = "insomnia_sr")
sri_df <- sri_df |> select(-starts_with("dur_day_total_"))

# Accelerometery data
accel_df <-
  fread(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_personsumMM_output.csv"))

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
#########################################################################################

# Sleep disorders data
sleep_dis_df <- fread(file.path(data_dir, "Sleep_disorders/sleep_disorders_selfreport_primarycare_hosp.csv"))

# Merge All dataframes
sri_df <- left_join(sri_df, accel_df, by = "eid")
sri_df <- left_join(sri_df, sleep_dis_df, by = "eid")

mri_df <- left_join(mri_df, accel_df, by = "eid")
mri_df <- left_join(mri_df, mri_qc_df, by = "eid")
mri_df <- left_join(mri_df, sleep_dis_df, by = "eid")

clean_base_df <- function(df) {
  # sleep disorder prior to accelerometry?
  df$OSA_dx <- ifelse(df$OSA_sr == 1 | df$date_osa_dx <= df$accel_date, 1, 0)
  df$insomnia_dx <-
    ifelse(df$insomnia_sr == 1 | df$date_insomnia_dx <= df$accel_date, 1, 0)
  df$sleep_disorder_dx <-
    ifelse(df$sleep_disorder_sr == 1 | df$date_any_sleep_dx <= df$accel_date, 1, 0)

  # replace missing with zero
  df <- df |> mutate(across(c(OSA_dx, insomnia_dx, sleep_disorder_dx), ~ if_else(is.na(.), 0, .)))

  # Death from other causes censoring event
  df <- df |>
    mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |>
    mutate(
      time_to_dem =
        case_when(
          dem == 1 ~ difftime(date_acdem2, accel_date),
          competing == 1 ~ difftime(date_of_death, accel_date),
          TRUE ~ difftime("2022-01-01", accel_date)
        )
    ) |>
    mutate(time_to_dem = as.integer(time_to_dem))

  # set pack years to zero for non-smokers
  df$smok_pckyrs <- ifelse(df$smok_status == 0, 0, df$smok_pckyrs)

  # set reference categories for factors
  df$diagnosed_diabetes <- as.factor(df$diagnosed_diabetes)
  levels(df$diagnosed_diabetes) <- c("prefer not answer", "dont know", "no", "yes")
  df$smok_status <- as.factor(df$smok_status)
  levels(df$smok_status) <- c("prefer not answer", "never", "former", "current")
  df$apoe_e4 <- as.factor(df$apoe_e4)
  df$highest_qual <- fct_relevel(df$highest_qual, "Grad")
  df$ethnicity <- fct_relevel(df$ethnicity, "white")
  df$insomnia_scale_sr <- fct_relevel(df$insomnia_scale_sr, "never")
  df$chronotype <- fct_relevel(df$chronotype, "morning")
  df$freq_depressed_twoweeks <- fct_relevel(df$freq_depressed_twoweeks, "not at all")
  df$avg_total_household_income <- fct_relevel(df$avg_total_household_income, "31-50")
  df$sick_disabled <- ifelse(df$employment == "sick or disabled", 1, 0)
  df$retired <- ifelse(df$employment == "retired", 1, 0)
  df$overall_health_rating <- as.factor(df$overall_health_rating)
  levels(df$overall_health_rating) <- c("prefer not answer", "dont know", "excellent", "good", "fair", "poor")
  df$shift <- ifelse(df$job_night_shift %in% c(3, 4) | df$job_shift_work %in% c(3, 4), 1, 0)

  # mark prefer not answer as missing
  df <- df |>
    mutate(across(where(is.factor), as.character)) |>
    mutate(across(where(is.character), ~ na_if(., "prefer not answer"))) |>
    mutate(across(where(is.character), ~ na_if(., "dont know"))) |>
    mutate(across(where(is.character), ~ na_if(., ""))) |>
    mutate(across(where(is.character), as.factor)) |>
    mutate(across(where(is.factor), fct_drop))

  # TODO: Consult data dictionary to figure out why totals we are calculating don't match existing totals
  df$awake_sleep <-
    df$dur_spt_wake_IN_min_pla +
    df$dur_spt_wake_LIG_min_pla +
    df$dur_spt_wake_MOD_min_pla +
    df$dur_spt_wake_VIG_min_pla

  # TODO: Deal with these columns being duplicated from join
  df$mins_worn <-
    df$dur_spt_sleep_min_pla +
    df$dur_day_total_IN_min_pla +
    df$dur_day_total_LIG_min_pla +
    df$dur_day_total_MOD_min_pla +
    df$dur_day_total_VIG_min_pla +
    df$awake_sleep

  # Variables normalised to 1440 relative to their proportion of total wear time
  df$sleep_n <-
    (df$dur_spt_sleep_min_pla / df$mins_worn) * mins_in_day
  df$inactive_n <-
    ((df$dur_day_total_IN_min_pla + df$dur_spt_wake_IN_min_pla) / df$mins_worn) * mins_in_day
  df$light_n <-
    ((df$dur_day_total_LIG_min_pla + df$dur_spt_wake_LIG_min_pla) / df$mins_worn) * mins_in_day
  df$moderate_n <-
    ((df$dur_day_total_MOD_min_pla + df$dur_spt_wake_MOD_min_pla) / df$mins_worn) * mins_in_day
  df$vigorous_n <-
    ((df$dur_day_total_VIG_min_pla + df$dur_spt_wake_VIG_min_pla) / df$mins_worn) * mins_in_day

  df$mvpa_n <- df$moderate_n + df$vigorous_n

# Update age variable to age at accelerometry study
  df$age_accel <-
    df$age_assessment + ((as.Date(df$accel_date) - as.Date(df$date_baseline)) / 365)
  df$age_accel <-
    as.numeric(df$age_accel)
  return(df)
}

mri_df <- clean_base_df(mri_df)
sri_df <- clean_base_df(sri_df)

# MRI df prep
###################################################
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

sri_comp <- acomp(data.frame(sri_df$sleep_n, sri_df$inactive_n, sri_df$light_n, sri_df$mvpa_n))
mri_comp <- acomp(data.frame(mri_df$sleep_n, mri_df$inactive_n, mri_df$light_n, mri_df$mvpa_n))
mri_base_ilr <-
  ilr(mri_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))
sri_base_ilr <-
  ilr(sri_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

# Model data for time_to_dem model
sri_model_data <- select(
  sri_df,
  sleep_n, inactive_n, light_n, mvpa_n,
  dem,
  time_to_dem,
  bp_syst_avg,
  avg_sri,
  age_assessment,
  retired,
  shift,
  townsend_deprivation_index,
  sex,
  antidepressant_med,
  antipsychotic_med,
  insomnia_med,
  avg_sleep_duration,
  sleep_duration_sr,
  ethnicity,
  avg_total_household_income,
  bp_med,
  any_cvd,
  parent_dementa,
  chronotype,
  sick_disabled,
  highest_qual,
  pa_modvig,
  smok_pckyrs,
  apoe_e4,
  freq_depressed_twoweeks,
  diagnosed_diabetes,
  BMI,
  smok_status,
) |> cbind(sri_base_ilr)

# Only run this once and then save the output

imp <- aregImpute(
  ~ R1 + R2 + R3 + dem + time_to_dem + bp_syst_avg + parent_dementa + bp_med + any_cvd + highest_qual +
    avg_sleep_duration + apoe_e4 + BMI + antidepressant_med + antipsychotic_med + insomnia_med + I(sleep_duration_sr) +
    ethnicity + age_assessment + sex + chronotype + freq_depressed_twoweeks + shift +
    avg_total_household_income + townsend_deprivation_index + retired + sick_disabled + avg_sri + smok_status +
    smok_pckyrs + pa_modvig + diagnosed_diabetes,
  data = model_data, n.impute = 10, type = "pmm"
)

# save imputed datasets
saveRDS(imp, "imp_sri.rds")

imp <- readRDS("imp_sri.rds")


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
  mean_tfmri_headmot,
  scan_lat_bpos,
  scan_trans_bpos,
  scan_long_bpos,
  scan_tabpos,
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

tbv_model <- lm(
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

model <- fit.mult.impute(
  Surv(time_to_dem, dem) ~
    rcs(R1, 3) + rcs(R2, 3) + rcs(R3, 3) +
    rcs(avg_sri, c(47, 61, 72)) +
    rcs(age_assessment, c(44, 57, 66)) +
    retired +
    shift +
    rcs(townsend_deprivation_index, c(-5, -2, 2.5)) +
    sex +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    highest_qual +
    apoe_e4 +
    smok_status,
  fitter = cph, xtrans = imp, data = sri_model_data
)

model <- coxph(
  Surv(time_to_dem, dem) ~
    rcs(R1, c(0.86, 1.29, 1.80)) + rcs(R2, c(-0.6809912, -0.4121788, -0.2104333)) + rcs(R3, c(0.11, 0.45, 0.84)) +
    rcs(avg_sri, c(47, 61, 72)) +
    rcs(age_assessment, c(44, 57, 66)) +
    retired +
    shift +
    rcs(townsend_deprivation_index, c(-5, -2, 2.5)) +
    sex +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    highest_qual +
    apoe_e4 +
    smok_status,
  data = model_data
)

sri_comp_df <- as.data.frame(sri_comp)
sri_comp_model <- lm(cbind(sri_df.sleep_n, sri_df.inactive_n, sri_df.light_n) ~ 1, data = sri_comp_df)
vcv_mat <- vcov(sri_comp_model)
val <- dmvnorm(sri_comp_df[1, -4], mean = coef(sri_comp_model), sigma = vcv_mat, log = TRUE)
sri_comp_df$dens <- apply(sri_comp_df[, -4], 1, dmvnorm, mean = coef(sri_comp_model), sigma = vcv_mat, log = TRUE)
sorted_df <- sri_comp_df[order(sri_comp_df$dens, decreasing = TRUE), ]
sorted_df[, 1:4] <- sorted_df[, 1:4] * 24

threshold <- quantile(sri_comp_df$dens, probs = c(0.025), na.rm = TRUE)

step_size <- 15

quantiles <- apply(sri_comp, 2, function(column) {
  quantile(column, probs = c(0.025, 0.975), na.rm = TRUE)
})

# Convert into minutes
quantiles_in_minutes <- quantiles * 24 * 60

# Store in "lower" and "upper" vectors
lower <- quantiles_in_minutes[1, ]
upper <- quantiles_in_minutes[2, ]

generate_compositions <- function(lower, upper) {
  # Round down for the lower bounds
  lower <- floor(lower / step_size) * step_size
  # Round up for the upper bounds
  upper <- ceiling(upper / step_size) * step_size

  # Create an empty data frame
  df <- expand.grid(
    sleep = seq(lower[1], upper[1], by = step_size),
    inactive = seq(lower[2], upper[2], by = step_size),
    light = seq(lower[3], upper[3], by = step_size),
    modvig = seq(lower[4], upper[4], by = step_size)
  )

  # Calculate total times
  df$total <- rowSums(df)

  # Filter out rows where total time is not 24 * 60
  df <- df[df$total == mins_in_day, ]
  rownames(df) <- NULL

  return(df)
}

generated_comps <- generate_compositions(lower, upper)
generated_comps <- generated_comps / mins_in_day
generated_comps$dens <- apply(generated_comps[, -c(4, 5)], 1, dmvnorm, mean = coef(sri_comp_model), sigma = vcv_mat, log = TRUE)
generated_comps <- generated_comps[generated_comps$dens > threshold, ]

generate_hazards <- function(comps, base_comp, base_ilr) {
  ilrs <- t(apply(comps, 1, function(comp) ilr(acomp(comp), V = v)))

  # Convert contrasts into a data frame
  contrasts <- t(apply(ilrs, 1, function(new_ilr) {
    contrast_out <- contrast(
      model,
      list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
      list(R1 = base_ilr[1], R2 = base_ilr[2], R3 = base_ilr[3])
    )
    return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
  }))
  contrasts_df <- as.data.frame(contrasts)
  colnames(contrasts_df) <- c("Contrast", "Lower", "Upper")

  # Add contrasts to comps data frame
  comps$Contrast <- contrasts_df$Contrast
  comps$Lower <- contrasts_df$Lower
  comps$Upper <- contrasts_df$Upper

  return(comps)
}

contrasts <- generate_hazards(generated_comps[, c(1, 2, 3, 4)], avg_sleep_geo_mean, sri_base_ilr)
sorted_contrasts <- contrasts[order(-contrasts$Contrast), ]
best_comp <- sorted_contrasts[nrow(sorted_contrasts), c(1, 2, 3, 4)]
worst_comp <- sorted_contrasts[1, c(1, 2, 3, 4)]
best_ilr <- ilr(acomp(best_comp))
worst_ilr <- ilr(acomp(worst_comp))

# Predicting average tbv for best and worst compositions
best_df <- mri_model_data
best_df$R1 <- rep(best_ilr[1], times = nrow(best_df))
best_df$R2 <- rep(best_ilr[2], times = nrow(best_df))
best_df$R3 <- rep(best_ilr[3], times = nrow(best_df))
best_results <- predict(tbv_model, newdata = best_df)
best_avg <- mean(best_results, na.rm = TRUE)

worst_df <- mri_model_data
worst_df$R1 <- rep(worst_ilr[1], times = nrow(worst_df))
worst_df$R2 <- rep(worst_ilr[2], times = nrow(worst_df))
worst_df$R3 <- rep(worst_ilr[3], times = nrow(worst_df))
worst_results <- predict(tbv_model, newdata = worst_df)
worst_avg <- mean(worst_results, na.rm = TRUE)

mri_avg_sleep_comp <-
  mri_comp[mri_comp$mri_df.sleep_n >= short_sleep_hours / hrs_in_day & mri_comp$mri_df.sleep_n <= long_sleep_hours / hrs_in_day]

avg_sleep_comp <-
  sri_comp[sri_comp$sri_df.sleep_n >= short_sleep_hours / hrs_in_day & sri_comp$sri_df.sleep_n <= long_sleep_hours / hrs_in_day]
short_sleep_comp <-
  sri_comp[sri_comp$sri_df.sleep_n < short_sleep_hours / hrs_in_day]
long_sleep_comp <-
  sri_comp[sri_comp$sri_df.sleep_n > long_sleep_hours / hrs_in_day]

geo_mean <-
  apply(sri_comp, 2, function(x) exp(mean(log(x))))

avg_sleep_geo_mean <-
  apply(avg_sleep_comp, 2, function(x) exp(mean(log(x))))
mri_avg_sleep_geo_mean <-
  apply(mri_avg_sleep_comp, 2, function(x) exp(mean(log(x))))
short_sleep_geo_mean <-
  apply(short_sleep_comp, 2, function(x) exp(mean(log(x))))
long_sleep_geo_mean <-
  apply(long_sleep_comp, 2, function(x) exp(mean(log(x))))

avg_ilr <-
  ilr(acomp(avg_sleep_geo_mean), V = v)
mri_avg_ilr <-
  ilr(acomp(mri_avg_sleep_geo_mean), V = v)
short_sleep_mean_ilr <-
  ilr(acomp(short_sleep_geo_mean), V = v)
long_sleep_mean_ilr <-
  ilr(acomp(long_sleep_geo_mean), V = v)

inc <- -sub_amount:sub_amount / mins_in_day
avg_offset <-
  matrix(rep(avg_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
mri_avg_offset <-
  matrix(rep(mri_avg_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
short_sleep_base_offset <-
  matrix(rep(short_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
long_sleep_base_offset <-
  matrix(rep(long_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)

# List of replacement prefixes and their respective pairs of indexes
# These are then used to generate and plot the relevant contrasts

sub_sleep_conditions <- list(
  inactive_sleep = c(inactive_idx, sleep_idx),
  light_sleep = c(light_idx, sleep_idx),
  vig_sleep = c(vig_idx, sleep_idx)
)

sub_inactive_conditions <- list(
  sleep_inactive = c(sleep_idx, inactive_idx),
  light_inactive = c(light_idx, inactive_idx),
  vig_inactive = c(vig_idx, inactive_idx)
)

sub_light_conditions <- list(
  inactive_light = c(inactive_idx, light_idx),
  sleep_light = c(sleep_idx, light_idx),
  vig_light = c(vig_idx, light_idx)
)

sub_vig_conditions <- list(
  inactive_vig = c(inactive_idx, vig_idx),
  light_vig = c(light_idx, vig_idx),
  sleep_vig = c(sleep_idx, vig_idx)
)

# generate_plot <- function(conditions, base_offset, mean_ilr) {
generate_plot <- function(title, substitution, cohorts, model) {
  # Create an empty data frame
  df <- data.frame()
  scaled_x <- inc * mins_in_day

  # Process each condition
  for (cohort in cohorts) {
    # Create offset matrix
    adjusted_comp <- cohort$base_comp
    mean_ilr <- cohort$ilr

    adjusted_comp[, substitution[1]] <- adjusted_comp[, substitution[1]] + inc

    other_idx <- setdiff(1:ncol(cohort$base_comp), substitution[1])

    for (i in other_idx) {
      adjusted_comp[, i] <-
        adjusted_comp[, i] - (inc * (adjusted_comp[, i] / rowSums(adjusted_comp[, -substitution[1]])))
    }
    # adjusted_comp[, substitution[2]] <- adjusted_comp[, substitution[2]] - inc

    # filtered_indices <- any(adjusted_comp <= 0)
    # adjusted_comp <- adjusted_comp[filtered_indices,]
    # scaled_x <- scaled_x[filtered_indices] * mins_in_day

    # Calculate ilrs
    ilrs <- t(apply(adjusted_comp, 1, function(comp) ilr(acomp(comp), V = v)))

    # Get contrasts
    contrasts <- apply(ilrs, 1, function(new_ilr) {
      contrast_out <- contrast(
        model,
        list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
        list(R1 = mean_ilr[1], R2 = mean_ilr[2], R3 = mean_ilr[3])
      )
      return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
    })

    df <- rbind(df, data.frame(
      X = scaled_x,
      Y = contrasts[1, ],
      Y_lower = contrasts[2, ],
      Y_upper = contrasts[3, ],
      name = cohort$name
    ))
  }

  return(ggplot(df, aes(x = X)) +
    geom_line(aes(y = Y, color = name)) +
    geom_ribbon(aes(ymin = Y_lower, ymax = Y_upper, fill = name), alpha = 0.2) +
    facet_wrap(~name) +
    xlab("min/day realocated") +
    ylab("Delta log hazard") +
    geom_hline(yintercept = 0) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = "darkblue")))
}

sub_sleep_flat_plot <-
  flat_sub(avg_offset, avg_ilr)
short_sleep_flat_plot <-
  flat_sub(short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_flat_plot <-
  flat_sub(long_sleep_base_offset, long_sleep_mean_ilr)

cohorts <- list(
  list(
    name = "1 Short Sleepers",
    base_comp = short_sleep_base_offset,
    ilr = short_sleep_mean_ilr
  ),
  list(
    name = "2 Average Sleepers",
    base_comp = avg_offset,
    ilr = avg_ilr
  ),
  list(
    name = "3 Long Sleepers",
    base_comp = long_sleep_base_offset,
    ilr = long_sleep_mean_ilr
  )
)

mri_cohorts <- list(
  list(
    name = "Average Sleepers",
    base_comp = mri_avg_offset,
    ilr = mri_avg_ilr
  )
)

sleep_inactive_plot <-
  generate_plot("Substituting sleep for inactivity", c(sleep_idx, 0), cohorts)
sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), cohorts)
sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), cohorts)

mri_sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), mri_cohorts, tbv_model)
mri_sleep_inactive_plot <-
  generate_plot("Substituting sleep for inactivity", c(sleep_idx, inactive_idx), mri_cohorts, tbv_model)
mri_sleep_light_plot <-
  generate_plot("Substituting sleep for light", c(sleep_idx, light_idx), mri_cohorts, tbv_model)

sub_sleep_avg_plot <-
  generate_plot(sub_sleep_conditions, avg_offset, avg_ilr)
short_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_inactive_avg_plot <-
  generate_plot(sub_inactive_conditions, avg_offset, avg_ilr)
short_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_light_avg_plot <-
  generate_plot(sub_light_conditions, avg_offset, avg_ilr)
short_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_vig_avg_plot <-
  generate_plot(sub_vig_conditions, avg_offset, avg_ilr)
short_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

ggsave(filename = "plots/flat_avg_plot.png", plot = sub_sleep_flat_plot, bg = "white")
ggsave(filename = "plots/flat_short_sleep_plot.png", plot = short_sleep_flat_plot, bg = "white")
ggsave(filename = "plots/flat_long_sleep_plot.png", plot = long_sleep_flat_plot, bg = "white")

ggsave(filename = "plots/avg_sub_sleep_plot.png", plot = sub_sleep_avg_plot)
ggsave(filename = "plots/short_sleep_sub_sleep_plot.png", plot = short_sleep_sub_sleep_plot)
ggsave(filename = "plots/long_sleep_sub_sleep_plot.png", plot = long_sleep_sub_sleep_plot)

ggsave(filename = "plots/avg_sub_inactive_plot.png", plot = sub_inactive_avg_plot)
ggsave(filename = "plots/short_sleep_sub_inactive_plot.png", plot = short_sleep_sub_inactive_plot)
ggsave(filename = "plots/long_sleep_sub_inactive_plot.png", plot = long_sleep_sub_inactive_plot)

ggsave(filename = "plots/avg_sub_light_plot.png", plot = sub_light_avg_plot)
ggsave(filename = "plots/short_sleep_sub_light_plot.png", plot = short_sleep_sub_light_plot)
ggsave(filename = "plots/long_sleep_sub_light_plot.png", plot = long_sleep_sub_light_plot)

ggsave(filename = "plots/avg_sub_vig_plot.png", plot = sub_vig_avg_plot)
ggsave(filename = "plots/short_sleep_sub_vig_plot.png", plot = short_sleep_sub_vig_plot)
ggsave(filename = "plots/long_sleep_sub_vig_plot.png", plot = long_sleep_sub_vig_plot)

