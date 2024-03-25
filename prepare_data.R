source("utils.R")

# load packages
library(tidyverse)
library(dotenv)
library(data.table)
library(compositions)
library(mice)
library(survival)

# read main dataset

dem_df <-
  fread(file.path(data_dir, "24hr_behaviours_25_03_24.csv"),
    stringsAsFactors = TRUE
  ) |>
  as_tibble()

# rename and remove some variables
dem_df <- dem_df |> rename("insomnia_scale_sr" = "insomnia_sr")
dem_df <- dem_df |> select(-starts_with("dur_day_total_"))

# Sleep disorders data
sleep_dis_df <-
  fread(file.path(
    data_dir,
    "../Sleep_disorders/sleep_disorders_selfreport_primarycare_hosp.csv"
  ))

# Merge dataframes
dem_df <- left_join(dem_df, sleep_dis_df, by = "eid")

# sleep disorder prior to accelerometry?
dem_df$OSA_dx <-
  ifelse(dem_df$OSA_sr == 1 | dem_df$date_osa_dx <= dem_df$calendar_date, 1, 0)
dem_df$insomnia_dx <-
  ifelse(dem_df$insomnia_sr == 1 | dem_df$date_insomnia_dx <= dem_df$calendar_date, 1, 0)
dem_df$sleep_disorder_dx <-
  ifelse(dem_df$sleep_disorder_sr == 1 | dem_df$date_any_sleep_dx <= dem_df$calendar_date, 1, 0)

# replace missing with zero
dem_df <-
  dem_df |>
  mutate(across(c(OSA_dx, insomnia_dx, sleep_disorder_dx), ~ if_else(is.na(.), 0, .)))

### Prepare confounding variables ###
# set pack years to zero for non-smokers
dem_df$smok_pckyrs <- ifelse(dem_df$smok_status == 0, 0, dem_df$smok_pckyrs)

# set alcohol prefer not to answer to missing
dem_df$alc_freq <- na_if(dem_df$alc_freq, -3)

# set do not know and prefer not to answer as missing for diet variables
dem_df <- dem_df |> mutate(across(c("veg_cooked","veg_raw","fruit_fresh","fruit_dry"), ~ na_if(., -1)))
dem_df <- dem_df |> mutate(across(c("veg_cooked","veg_raw","fruit_fresh","fruit_dry"), ~ na_if(., -3)))

# set "less than one" to 0 for diet variables
dem_df <- dem_df |> mutate(across(c("veg_cooked","veg_raw","fruit_fresh","fruit_dry"), ~if_else(. == -10, 0, .)))

# create total fruit and veg variable 
dem_df$fruit_veg <- rowSums(dem_df[,c(c("veg_cooked","veg_raw","fruit_fresh","fruit_dry"))])

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

# Update age variable to age at accelerometry study
dem_df$age_accel <-
  dem_df$age_assessment + ((as.Date(dem_df$calendar_date) - as.Date(dem_df$date_baseline)) / 365)

dem_df$age_accel <- as.numeric(dem_df$age_accel)

dem_df <- dem_df |> rename("date_accel" = "calendar_date")

### Create isometric log ratio coordinates

# close time use variables

dem_comp <-
  acomp(data.frame(
    dem_df$avg_sleep, dem_df$avg_inactivity, dem_df$avg_light, dem_df$avg_mvpa
  ))

dem_base_ilr <-
  ilr(dem_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

## Construct risk set
# Death from other causes remain in risk set until end of FU

dem_df <- dem_df |>
  mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |>
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, date_accel),
    competing == 1 ~ difftime("2022-01-01", date_accel),
    TRUE ~ difftime("2022-01-01", date_accel)
  )) |>
  mutate(time_to_dem = as.integer(time_to_dem))

# Create age at dementia or censoring/end of follow-up variable

dem_df$age_dem <- dem_df$age_accel + (dem_df$time_to_dem / 365)

## Add WASO and SRI to dataset ##

waso_dat <- fread(file.path(data_dir, "../SRI/sri_data_may2023.csv"))

waso_dat <- select(waso_dat, eid, avg_WASO, avg_sri)

dem_df <- left_join(dem_df, waso_dat, by = "eid")

### Add prevalent disease variables

### Add in prevalent illness data ###

prev <- fread(file.path(data_dir, "../SRI/disease_dates.csv"),
              stringsAsFactors = T)

# add in date of actigraphy and number of cancers

adate <- dem_df |> select(eid, date_accel, num_sr_cancers)

prev <- left_join(prev, adate, by = "eid")

prev <- as_tibble(prev)

# create prevalent illness variables

prev <- prev |> mutate(across(starts_with("sr_"), ~ifelse(is.na(.), 0, .)))

prev$date_accel <- as_date(prev$date_accel)

prev <- prev |>
  mutate(prev_diabetes = case_when(sr_diabetes == 1 ~ 1,
                                   sr_diabetes == 0 & date_diabetes < date_accel ~ 1,
                                   TRUE ~ 0)) |>
  mutate(prev_cancer = case_when(num_sr_cancers > 0 ~ 1,
                                 num_sr_cancers == 0 & date_neoplasm < date_accel ~ 1,
                                 TRUE ~ 0)) |>
  mutate(prev_mental_disorder = case_when(sr_mental_disorder == 1 ~ 1,
                                          sr_mental_disorder == 0 & date_mental_disorder < date_accel ~ 1,
                                          TRUE ~ 0)) |>
  mutate(prev_nervous_system = case_when(sr_nervous_system == 1 ~ 1,
                                         sr_nervous_system == 0 & date_nervous_system_disorder < date_accel ~ 1,
                                         TRUE ~ 0)) |>
  mutate(prev_cvd = case_when(sr_cvd == 1 ~ 1,
                              sr_cvd == 0 & date_CVD < date_accel ~ 1,
                              TRUE ~ 0))

prev <- prev |> select(eid, starts_with("prev_"))

dem_df <- left_join(dem_df, prev, by = "eid")

### Select model data

dem_model_data <- select(
  dem_df,
  eid,
  fruit_veg,
  alc_freq,
  avg_sleep, avg_inactivity, avg_light, avg_mvpa,
  dem, time_to_dem,
  avg_WASO, avg_sri,
  bp_syst_avg,
  age_accel,
  retired,
  shift,
  townsend_deprivation_index,
  sex,
  antidepressant_med,
  antipsychotic_med,
  insomnia_med,
  ethnicity,
  avg_total_household_income,
  bp_med,
  any_cvd,
  parent_dementa,
  chronotype,
  sick_disabled,
  highest_qual,
  smok_pckyrs,
  apoe_e4,
  freq_depressed_twoweeks,
  diagnosed_diabetes,
  BMI,
  smok_status,
  date_of_death,
  age_at_death,
  date_acdem2,
  date_accel,
  age_dem,
  sick_disabled,
  prev_diabetes,
  prev_cancer,
  prev_mental_disorder,
  prev_nervous_system,
  prev_cvd,
  bp_med,
  BMI,
  bp_syst_avg
) |> cbind(dem_base_ilr)

## Output RDS file ##

write_rds(dem_model_data, file.path(data_dir, "bootstrap_data_25_03_24.rds"))

