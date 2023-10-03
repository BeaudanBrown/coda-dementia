library(compositions)
library(tidyverse)
library(dotenv)
library(data.table)
library(rms)

mins_in_day <- 1440

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

dem_df <-
  fread(file.path(data_dir, "SRI/sri_data_may2023.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
dem_df <- dem_df |> rename("insomnia_scale_sr" = "insomnia_sr")
dem_df <- dem_df |> select(-starts_with("dur_day_total_"))

# Accelerometery data
accel_df <-
  fread(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_personsumMM_output.csv"))

# Sleep disorders data
sleep_dis_df <- fread(file.path(data_dir, "Sleep_disorders/sleep_disorders_selfreport_primarycare_hosp.csv"))

# Merge All dataframes
dem_df <- left_join(dem_df, accel_df, by = "eid")
dem_df <- left_join(dem_df, sleep_dis_df, by = "eid")

# sleep disorder prior to accelerometry?
dem_df$OSA_dx <- ifelse(dem_df$OSA_sr == 1 | dem_df$date_osa_dx <= dem_df$accel_date, 1, 0)
dem_df$insomnia_dx <-
  ifelse(dem_df$insomnia_sr == 1 | dem_df$date_insomnia_dx <= dem_df$accel_date, 1, 0)
dem_df$sleep_disorder_dx <-
  ifelse(dem_df$sleep_disorder_sr == 1 | dem_df$date_any_sleep_dx <= dem_df$accel_date, 1, 0)

# replace missing with zero
dem_df <- dem_df |> mutate(across(c(OSA_dx, insomnia_dx, sleep_disorder_dx), ~ if_else(is.na(.), 0, .)))

# Death from other causes censoring event
dem_df <- dem_df |>
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
dem_df$smok_pckyrs <- ifelse(dem_df$smok_status == 0, 0, dem_df$smok_pckyrs)

# set reference categories for factors
dem_df$diagnosed_diabetes <- as.factor(dem_df$diagnosed_diabetes)
levels(dem_df$diagnosed_diabetes) <- c("prefer not answer", "dont know", "no", "yes")
dem_df$smok_status <- as.factor(dem_df$smok_status)
levels(dem_df$smok_status) <- c("prefer not answer", "never", "former", "current")
dem_df$apoe_e4 <- as.factor(dem_df$apoe_e4)
dem_df$highest_qual <- fct_relevel(dem_df$highest_qual, "Grad")
dem_df$ethnicity <- fct_relevel(dem_df$ethnicity, "white")
dem_df$insomnia_scale_sr <- fct_relevel(dem_df$insomnia_scale_sr, "never")
dem_df$chronotype <- fct_relevel(dem_df$chronotype, "morning")
dem_df$freq_depressed_twoweeks <- fct_relevel(dem_df$freq_depressed_twoweeks, "not at all")
dem_df$avg_total_household_income <- fct_relevel(dem_df$avg_total_household_income, "31-50")
dem_df$sick_disabled <- ifelse(dem_df$employment == "sick or disabled", 1, 0)
dem_df$retired <- ifelse(dem_df$employment == "retired", 1, 0)
dem_df$overall_health_rating <- as.factor(dem_df$overall_health_rating)
levels(dem_df$overall_health_rating) <- c("prefer not answer", "dont know", "excellent", "good", "fair", "poor")
dem_df$shift <- ifelse(dem_df$job_night_shift %in% c(3, 4) | dem_df$job_shift_work %in% c(3, 4), 1, 0)

# mark prefer not answer as missing
dem_df <- dem_df |>
  mutate(across(where(is.factor), as.character)) |>
  mutate(across(where(is.character), ~ na_if(., "prefer not answer"))) |>
  mutate(across(where(is.character), ~ na_if(., "dont know"))) |>
  mutate(across(where(is.character), ~ na_if(., ""))) |>
  mutate(across(where(is.character), as.factor)) |>
  mutate(across(where(is.factor), fct_drop))

# TODO: Consult data dictionary to figure out why totals we are calculating don't match existing totals
dem_df$awake_sleep <-
  dem_df$dur_spt_wake_IN_min_pla +
  dem_df$dur_spt_wake_LIG_min_pla +
  dem_df$dur_spt_wake_MOD_min_pla +
  dem_df$dur_spt_wake_VIG_min_pla

# TODO: Deal with these columns being duplicated from join
dem_df$mins_worn <-
  dem_df$dur_spt_sleep_min_pla +
  dem_df$dur_day_total_IN_min_pla +
  dem_df$dur_day_total_LIG_min_pla +
  dem_df$dur_day_total_MOD_min_pla +
  dem_df$dur_day_total_VIG_min_pla +
  dem_df$awake_sleep

# Variables normalised to 1440 relative to their proportion of total wear time
dem_df$sleep_n <-
  (dem_df$dur_spt_sleep_min_pla / dem_df$mins_worn) * mins_in_day
dem_df$inactive_n <-
  ((dem_df$dur_day_total_IN_min_pla + dem_df$dur_spt_wake_IN_min_pla) / dem_df$mins_worn) * mins_in_day
dem_df$light_n <-
  ((dem_df$dur_day_total_LIG_min_pla + dem_df$dur_spt_wake_LIG_min_pla) / dem_df$mins_worn) * mins_in_day
dem_df$moderate_n <-
  ((dem_df$dur_day_total_MOD_min_pla + dem_df$dur_spt_wake_MOD_min_pla) / dem_df$mins_worn) * mins_in_day
dem_df$vigorous_n <-
  ((dem_df$dur_day_total_VIG_min_pla + dem_df$dur_spt_wake_VIG_min_pla) / dem_df$mins_worn) * mins_in_day

dem_df$mvpa_n <- dem_df$moderate_n + dem_df$vigorous_n

# Update age variable to age at accelerometry study
dem_df$age_accel <-
  dem_df$age_assessment + ((as.Date(dem_df$accel_date) - as.Date(dem_df$date_baseline)) / 365)
dem_df$age_accel <-
  as.numeric(dem_df$age_accel)

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)
v <- gsi.buildilrBase(t(sbp))

dem_comp <- acomp(data.frame(dem_df$sleep_n, dem_df$inactive_n, dem_df$light_n, dem_df$mvpa_n))

dem_base_ilr <-
  ilr(dem_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

# Model data for time_to_dem model
dem_model_data <- select(
  dem_df,
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
) |> cbind(dem_base_ilr)

# Only run this once and then save the output
# imp <- aregImpute(
#   ~ R1 + R2 + R3 + dem + time_to_dem + bp_syst_avg + parent_dementa + bp_med + any_cvd + highest_qual +
#     avg_sleep_duration + apoe_e4 + BMI + antidepressant_med + antipsychotic_med + insomnia_med + I(sleep_duration_sr) +
#     ethnicity + age_assessment + sex + chronotype + freq_depressed_twoweeks + shift +
#     avg_total_household_income + townsend_deprivation_index + retired + sick_disabled + avg_sri + smok_status +
#     smok_pckyrs + pa_modvig + diagnosed_diabetes,
#   data = model_data, n.impute = 10, type = "pmm"
# )

# # save imputed datasets
# saveRDS(imp, "imp_sri.rds")

imp <- readRDS("imp_sri.rds")

dem_model <- fit.mult.impute(
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
  fitter = cph, xtrans = imp, data = dem_model_data
)

# dem_model <- coxph(
#   Surv(time_to_dem, dem) ~
#     rcs(R1, c(0.86, 1.29, 1.80)) + rcs(R2, c(-0.6809912, -0.4121788, -0.2104333)) + rcs(R3, c(0.11, 0.45, 0.84)) +
#     rcs(avg_sri, c(47, 61, 72)) +
#     rcs(age_assessment, c(44, 57, 66)) +
#     retired +
#     shift +
#     rcs(townsend_deprivation_index, c(-5, -2, 2.5)) +
#     sex +
#     antidepressant_med +
#     antipsychotic_med +
#     insomnia_med +
#     ethnicity +
#     avg_total_household_income +
#     highest_qual +
#     apoe_e4 +
#     smok_status,
#   data = model_data
# )
