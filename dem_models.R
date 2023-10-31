# load packages
# library(mice)
# library(tidyverse)
# library(survival)
# library(rms)
# library(compositions)
# library(data.table)
# library(boot)

list_of_packages <- c(
  "mice",
  "tidyverse",
  "survival",
  "rms",
  "compositions",
  "data.table",
  "parallel",
  "boot"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

if (length(new_packages)) install.packages(new_packages)

lapply(list_of_packages, library, character.only = TRUE)

# constants
short_sleep_hours <- 6
hrs_in_day <- 24

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

# source helper functions
source("model_helpers.R")

## Impute ##

# imp <- mice(boot_data, predictorMatrix = predmat,
#            method = imp_methods, m = 1)

# saveRDS(imp, file.path(data_dir, "full_imp.rds"))

imp <- read_rds(file.path(data_dir, "full_imp.rds"))

imp <- complete(imp)
imp$id <- seq_len(nrow(imp))

primary_form <- as.formula(dem ~ rcs(timegroup, 5) +
    poly(R1, 2) +
    poly(R2, 2) +
    poly(R3, 2) +
    sex +
    retired +
    shift +
    apoe_e4 +
    highest_qual +
    rcs(townsend_deprivation_index, 3) +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    smok_status
)

model <- fit_model(imp, primary_form)

## Get results for substitutions

plot_data <- sub_results(model, imp, timegroup = 55)
# write_rds(plot_data, file.path(data_dir,"primary_mod_res.rds"))
plot_data <- read_rds(file.path(data_dir, "primary_mod_res.rds"))

## Plot ##

plot1 <- plot_data |>
  ggplot(aes(x = offset, y = risk, colour = Reference)) +
  geom_line() +
  facet_wrap(~Substitution) +
  cowplot::theme_cowplot() +
  labs(x = "Time added to sleep", y = "Dementia cumulative incidence by age 75")

plot1

#### Primary model - bootstrap ####

## bootstrap

ncpus <- 6
boot_out <- boot(
  data = boot_data, statistic = boot_fn, reg_formula = primary_form,
  timegroup = timegroup, R = ncpus, parallel = "multicore", ncpus = ncpus
)


#### Figure 1 ####

# combine bootstrap and full sample estimates

#### Sensitivity 1 - full sample ####

s1_form <- as.formula(dem ~ rcs(timegroup, 5) +
    poly(R1, 2) +
    poly(R2, 2) +
    poly(R3, 2) +
    sex +
    retired +
    shift +
    apoe_e4 +
    highest_qual +
    rcs(townsend_deprivation_index, 3) +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    smok_status +
    rcs(avg_WASO, 3)
)

model_s1 <- fit_model(imp, s1_form)

## Get substitution results

plot_data_s1 <- sub_results(model_s1, imp)

plot_s1 <-
  plot_data_s1 |>
  ggplot(aes(x = offset, y = risk, colour = Reference)) +
  geom_line() +
  facet_wrap(~Substitution) +
  cowplot::theme_cowplot() +
  labs(x = "Time added to sleep", y = "Dementia cumulative incidence by age 75") +
  ylim(c(0, 0.05))

#### Sensitivity 1 - bootstrap ####

ncpus <- 6
boot_out_s1 <- boot(
  data = boot_data, statistic = boot_fn, reg_formula = s1_form,
  R = ncpus, parallel = "multicore", ncpus = ncpus
)


#### Sensitivity 2 - full sample ####

s2_form <- as.formula(dem ~ rcs(timegroup, 5) +
    poly(R1, 2) +
    poly(R2, 2) +
    poly(R3, 2) +
    sex +
    retired +
    shift +
    apoe_e4 +
    highest_qual +
    rcs(townsend_deprivation_index, 3) +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    smok_status +
    sick_disabled +
    prev_diabetes +
    prev_cancer +
    prev_mental_disorder +
    prev_nervous_system +
    prev_cvd +
    bp_med +
    rcs(BMI, c(22, 26, 32)) +
    rcs(bp_syst_avg, c(115, 135, 161))
)

model_s2 <- fit_model(imp, s2_form)

## Get substitution results

plot_data_s2 <- sub_results(model_s2, imp)

#### Sensitivity 2 - bootstrap ####

ncpus <- 6
boot_out_s2 <- boot(
  data = boot_data, statistic = boot_fn, reg_formula = s2_form,
  R = ncpus, parallel = "multicore", ncpus = ncpus
)


#### Sensitivity 3 - full sample ####

## Create dataset with time since accelerometry as the timescale ##

imp_long <- survSplit(Surv(time = time_to_dem, event = dem) ~ .,
  data = imp,
  cut = seq(
    from = min(imp$time_to_dem),
    to = max(imp$time_to_dem),
    length.out = 18
  ),
  episode = "timegroup", end = "time_start", event = "dem",
  start = "time_end"
)


model_s3 <- glm(dem ~ rcs(timegroup, 5) * (poly(R1, 2) + poly(R2, 2) + poly(R3, 2)) +
    sex +
    retired +
    shift +
    apoe_e4 +
    highest_qual +
    rcs(townsend_deprivation_index, 3) +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    smok_status, data = imp_long, family = binomial
)
