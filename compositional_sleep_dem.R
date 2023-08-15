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

sleep_idx <- 1
inactive_idx <- 2
light_idx <- 3
vig_idx <- 4
sub_amount <- 60
mins_in_day <- 1440
hours_in_day <- 24

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

# Variables and source data
source_df <-
  fread(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_personsumMM_output.csv"))

# Merge in outcome data
outcome_df <-
  fread(file.path(data_dir, "SRI/sri_data_may2023.csv"), stringsAsFactors = TRUE) |>
  as_tibble() |>
  select(
    eid, dem, time_to_dem, date_of_death, sex, age_assessment, accel_date, date_acdem2, smok_pckyrs, smok_status,
    bp_syst_avg, bp_med, any_cvd, diagnosed_diabetes, highest_qual, avg_sleep_duration,
    apoe_e4, BMI, antidepressant_med, antipsychotic_med, pa_modvig,
    insomnia_med, sleep_duration_sr, ethnicity,
    insomnia_sr, chronotype, freq_depressed_twoweeks,
    parent_dementa, avg_total_household_income, townsend_deprivation_index,
    employment, avg_sri, avg_sleep_duration, neurological_exclude_first2fu,
    job_shift_work, job_night_shift, overall_health_rating
  )

outcome_df <- outcome_df |> rename("insomnia_scale_sr" = "insomnia_sr")

sleep_dis <- fread(file.path(data_dir, "Sleep_disorders/sleep_disorders_selfreport_primarycare_hosp.csv"))

outcome_df <- left_join(outcome_df, sleep_dis, by = "eid")

# sleep disorder prior to accelerometry?

outcome_df$OSA_dx <- ifelse(outcome_df$OSA_sr == 1 | outcome_df$date_osa_dx <= outcome_df$accel_date, 1, 0)

outcome_df$insomnia_dx <-
  ifelse(outcome_df$insomnia_sr == 1 | outcome_df$date_insomnia_dx <= outcome_df$accel_date, 1, 0)

outcome_df$sleep_disorder_dx <-
  ifelse(outcome_df$sleep_disorder_sr == 1 | outcome_df$date_any_sleep_dx <= outcome_df$accel_date, 1, 0)

# replace missing with zero

outcome_df <- outcome_df |> mutate(across(c(OSA_dx, insomnia_dx, sleep_disorder_dx), ~ if_else(is.na(.), 0, .)))

### Death from other causes censoring event ###

outcome_df <- outcome_df |>
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

###################################################

# set pack years to zero for non-smokers

outcome_df$smok_pckyrs <- ifelse(outcome_df$smok_status == 0, 0, outcome_df$smok_pckyrs)

# set reference categories for factors

outcome_df$diagnosed_diabetes <- as.factor(outcome_df$diagnosed_diabetes)
levels(outcome_df$diagnosed_diabetes) <- c("prefer not answer", "dont know", "no", "yes")
outcome_df$smok_status <- as.factor(outcome_df$smok_status)
levels(outcome_df$smok_status) <- c("prefer not answer", "never", "former", "current")
outcome_df$apoe_e4 <- as.factor(outcome_df$apoe_e4)
outcome_df$highest_qual <- fct_relevel(outcome_df$highest_qual, "Grad")
outcome_df$ethnicity <- fct_relevel(outcome_df$ethnicity, "white")
outcome_df$insomnia_scale_sr <- fct_relevel(outcome_df$insomnia_scale_sr, "never")
outcome_df$chronotype <- fct_relevel(outcome_df$chronotype, "morning")
outcome_df$freq_depressed_twoweeks <- fct_relevel(outcome_df$freq_depressed_twoweeks, "not at all")
outcome_df$avg_total_household_income <- fct_relevel(outcome_df$avg_total_household_income, "31-50")
outcome_df$sick_disabled <- ifelse(outcome_df$employment == "sick or disabled", 1, 0)
outcome_df$retired <- ifelse(outcome_df$employment == "retired", 1, 0)
outcome_df$overall_health_rating <- as.factor(outcome_df$overall_health_rating)
levels(outcome_df$overall_health_rating) <- c("prefer not answer", "dont know", "excellent", "good", "fair", "poor")
outcome_df$shift <- ifelse(outcome_df$job_night_shift %in% c(3, 4) | outcome_df$job_shift_work %in% c(3, 4), 1, 0)
outcome_df <- outcome_df |> select(-job_night_shift, -job_shift_work)

# mark prefer not answer as missing

factor_vars <- Hmisc::Cs(
  smok_status, diagnosed_diabetes, highest_qual, ethnicity,
  chronotype, freq_depressed_twoweeks,
  avg_total_household_income, employment, overall_health_rating
)

outcome_df <- outcome_df |>
  mutate(across(all_of(factor_vars), as.character)) |>
  mutate(across(all_of(factor_vars), ~ na_if(., "prefer not answer"))) |>
  mutate(across(all_of(factor_vars), ~ na_if(., "dont know"))) |>
  mutate(across(all_of(factor_vars), as.factor)) |>
  mutate(across(where(is.factor), fct_drop))

###################################################

# Death from other causes censoring event
full_df <- left_join(outcome_df, source_df, by = "eid") |>
  mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |>
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, accel_date),
    competing == 1 ~ difftime(date_of_death, accel_date),
    TRUE ~ difftime("2022-01-01", accel_date)
  )) |>
  mutate(time_to_dem = as.integer(time_to_dem))

# TODO: Consult data dictionary to figure out why totals we are calculating don't match existing totals
full_df$awake_sleep <-
  full_df$dur_spt_wake_IN_min_pla +
  full_df$dur_spt_wake_LIG_min_pla +
  full_df$dur_spt_wake_MOD_min_pla +
  full_df$dur_spt_wake_VIG_min_pla

full_df$mins_worn <-
  full_df$dur_spt_sleep_min_pla +
  full_df$dur_day_total_IN_min_pla +
  full_df$dur_day_total_LIG_min_pla +
  full_df$dur_day_total_MOD_min_pla +
  full_df$dur_day_total_VIG_min_pla +
  full_df$awake_sleep

# Variables normalised to 1440 relative to their proportion of total wear time
full_df$sleep_n <-
  (full_df$dur_spt_sleep_min_pla / full_df$mins_worn) * mins_in_day
full_df$inactive_n <-
  ((full_df$dur_day_total_IN_min_pla + full_df$dur_spt_wake_IN_min_pla) / full_df$mins_worn) * mins_in_day
full_df$light_n <-
  ((full_df$dur_day_total_LIG_min_pla + full_df$dur_spt_wake_LIG_min_pla) / full_df$mins_worn) * mins_in_day
full_df$moderate_n <-
  ((full_df$dur_day_total_MOD_min_pla + full_df$dur_spt_wake_MOD_min_pla) / full_df$mins_worn) * mins_in_day
full_df$vigorous_n <-
  ((full_df$dur_day_total_VIG_min_pla + full_df$dur_spt_wake_VIG_min_pla) / full_df$mins_worn) * mins_in_day

full_df$mvpa_n <- full_df$moderate_n + full_df$vigorous_n

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)
v <- gsi.buildilrBase(t(sbp))

inc <- 0:sub_amount / 1440
comp <- acomp(data.frame(full_df$sleep_n, full_df$inactive_n, full_df$light_n, full_df$mvpa_n))
base_ilr <-
  ilr(comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

model_data <- select(
  full_df,
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
) |> cbind(base_ilr)

options(datadist = datadist(model_data))

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
  fitter = cph, xtrans = imp, data = model_data
)

avg_sleep_comp <- comp[comp$full_df.sleep_n >= 6 / 25 & comp$full_df.sleep_n <= 9 / 25]
short_sleep_comp <- comp[comp$full_df.sleep_n < 6 / 25]
long_sleep_comp <- comp[comp$full_df.sleep_n > 9 / 25]

avg_sleep_geo_mean <-
  apply(avg_sleep_comp, 2, function(x) exp(mean(log(x))))
short_sleep_geo_mean <-
  apply(short_sleep_comp, 2, function(x) exp(mean(log(x))))
long_sleep_geo_mean <-
  apply(long_sleep_comp, 2, function(x) exp(mean(log(x))))

avg_ilr <-
  ilr(acomp(avg_sleep_geo_mean), V = v)
short_sleep_mean_ilr <-
  ilr(acomp(short_sleep_geo_mean), V = v)
long_sleep_mean_ilr <-
  ilr(acomp(long_sleep_geo_mean), V = v)

avg_offset <- matrix(rep(avg_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
short_sleep_base_offset <- matrix(rep(short_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
long_sleep_base_offset <- matrix(rep(long_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)

# List of replacement prefixes and their respective pairs of indexes
# These are then used to generate and plot the relevant contrasts
sub_sleep_conditions <- list(
  inactive_sleep = c(inactive_idx, sleep_idx),
  light_sleep = c(light_idx, sleep_idx),
  vig_sleep = c(vig_idx, sleep_idx),
  sleep_inactive = c(sleep_idx, inactive_idx),
  sleep_light = c(sleep_idx, light_idx),
  sleep_vig = c(sleep_idx, vig_idx)
)

sub_inactive_conditions <- list(
  sleep_inactive = c(sleep_idx, inactive_idx),
  light_inactive = c(light_idx, inactive_idx),
  vig_inactive = c(vig_idx, inactive_idx),
  inactive_sleep = c(inactive_idx, sleep_idx),
  inactive_light = c(inactive_idx, light_idx),
  inactive_vig = c(inactive_idx, vig_idx)
)

sub_light_conditions <- list(
  inactive_light = c(inactive_idx, light_idx),
  sleep_light = c(sleep_idx, light_idx),
  vig_light = c(vig_idx, light_idx),
  light_inactive = c(light_idx, inactive_idx),
  light_sleep = c(light_idx, sleep_idx),
  light_vig = c(light_idx, vig_idx)
)

sub_vig_conditions <- list(
  inactive_vig = c(inactive_idx, vig_idx),
  light_vig = c(light_idx, vig_idx),
  sleep_vig = c(sleep_idx, vig_idx),
  vig_inactive = c(vig_idx, inactive_idx),
  vig_light = c(vig_idx, light_idx),
  vig_sleep = c(vig_idx, sleep_idx)
)

generate_plot <- function(conditions, base_offset, mean_ilr, comp) {
  # Initialize your contrast list
  contrasts <- list()

  # Process each condition
  for (cond in names(conditions)) {
    # Create offset matrix
    offset <- base_offset
    offset[, conditions[[cond]][1]] <- base_offset[, conditions[[cond]][1]] + inc
    offset[, conditions[[cond]][2]] <- base_offset[, conditions[[cond]][2]] - inc

    # Calculate ilrs
    ilrs <- t(apply(offset, 1, function(comp) ilr(acomp(comp), V = v)))

    # Get contrasts
    contrasts[[cond]] <- apply(ilrs, 1, function(new_ilr) {
      contrast_out <- contrast(
        model,
        list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
        list(R1 = mean_ilr[1], R2 = mean_ilr[2], R3 = mean_ilr[3])
      )
      return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
    })
  }

  # Create an empty data frame
  df <- data.frame()

  scaled_x <- inc * 1440
  # Iterate over conditions and add data to the dataframe
  for (cond in names(conditions)) {
    temp_df <- data.frame(
      X = scaled_x,
      Y = contrasts[[cond]][1, ],
      Y_lower = contrasts[[cond]][2, ],
      Y_upper = contrasts[[cond]][3, ],
      condition = cond
    )
    df <- rbind(df, temp_df)
  }

  return(ggplot(df, aes(x = X)) +
    geom_line(aes(y = Y, color = condition)) +
    geom_ribbon(aes(ymin = Y_lower, ymax = Y_upper, fill = condition), alpha = 0.2) +
    facet_wrap(~condition) +
    xlab("min/day realocated") +
    ylab("Delta log hazard") +
    geom_hline(yintercept = 0))
}

sub_sleep_avg_plot <-
  generate_plot(sub_sleep_conditions, avg_offset, avg_ilr, comp)
short_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, short_sleep_base_offset, short_sleep_mean_ilr, comp)
long_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, long_sleep_base_offset, long_sleep_mean_ilr, comp)

sub_inactive_avg_plot <-
  generate_plot(sub_inactive_conditions, avg_offset, avg_ilr, comp)
short_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, short_sleep_base_offset, short_sleep_mean_ilr, comp)
long_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, long_sleep_base_offset, long_sleep_mean_ilr, comp)

sub_light_avg_plot <-
  generate_plot(sub_light_conditions, avg_offset, avg_ilr, comp)
short_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, short_sleep_base_offset, short_sleep_mean_ilr, comp)
long_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, long_sleep_base_offset, long_sleep_mean_ilr, comp)

sub_vig_avg_plot <-
  generate_plot(sub_vig_conditions, avg_offset, avg_ilr, comp)
short_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, short_sleep_base_offset, short_sleep_mean_ilr, comp)
long_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, long_sleep_base_offset, long_sleep_mean_ilr, comp)

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
