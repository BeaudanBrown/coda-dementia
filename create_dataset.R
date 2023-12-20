source("utils.R")

### Preparing UKB data ###
# packages
library(tidyverse)
library(lme4)
library(emmeans)
library(data.table)
library(lubridate)
library(here)

# read UKB data

d <- fread(file.path(data_dir, "../core_ukb_data/core_ukb_trimmed_Sep2023.csv"))

# Add in SNPs #

snps <- read_csv(file.path(data_dir, "../Raw UKB data/SNPs (basket 2)/snp_data.csv"))

d <- d |> left_join(snps, by = "eid")

#### Clean confounder variables ####

## continuous variables

# histograms
d |>
  select(
    BMI, num_non_cancer_illness, townsend_deprivation_index, sleep_duration_sr,
    neuroticism_score, bp_syst_avg, num_treatments_medications
  ) |>
  psych::multi.hist(global = F)

# remove negative values from sleep duration

d <-
  d |>
  mutate(sleep_duration_sr = ifelse(sleep_duration_sr < 0, NA_integer_, sleep_duration_sr))

### Categorical

d$apoe_e4 <- as.factor(d$apoe_e4)

d$bp_med <- as.factor(d$bp_med)

d$any_cvd <- as.factor(d$any_cvd)

# qualifications

quals <- d |>
  select(eid, contains("qualifications")) |>
  pivot_longer(-eid) |>
  mutate(highest_qual = case_when(
    value == 1 ~ 4,
    value == 6 ~ 4,
    value == 5 ~ 3,
    value == 2 ~ 2,
    value == 3 ~ 1,
    value == 4 ~ 1,
    TRUE ~ as.numeric(value)
  )) |>
  group_by(eid) |>
  summarise(highest_qual = max(highest_qual, na.rm = T)) |>
  mutate(highest_qual = ifelse(highest_qual == -Inf, NA_real_, highest_qual))

d <- left_join(d, quals, by = "eid")

d$highest_qual <- as.factor(d$highest_qual)

d$highest_qual <- relevel(d$highest_qual, "4") # grad qualification as reference

levels(d$highest_qual) <- c("Grad", "Other", "prefer not answer", "O", "A", "NVQ")

summary(d$highest_qual)

# remove redundant variables

d <- d |> select(-starts_with("qualification"))

## antidepressant medication

d$antidepressant_med <- as.factor(d$antidepressant_med)

summary(d$antidepressant_med)

## antipsychotic medication

d$antipsychotic_med <- as.factor(d$antipsychotic_med)

summary(d$antipsychotic_med)

## sleep medication

d$insomnia_med <- as.factor(d$insomnia_med)

summary(d$insomnia_med)

## Race

d$ethnicity <- as.factor(d$ethnicity)

# collapse categories

d$ethnicity <- fct_recode(d$ethnicity,
  "prefer not answer" = "-3",
  "other" = "-1",
  "white" = "1",
  "white" = "1001",
  "white" = "2001",
  "white" = "3001",
  "white" = "4001",
  "mixed" = "2",
  "mixed" = "1002",
  "mixed" = "2002",
  "mixed" = "3002",
  "mixed" = "4002",
  "asian" = "3",
  "asian" = "1003",
  "asian" = "2003",
  "asian" = "3003",
  "black" = "4",
  "black" = "2004",
  "black" = "3004",
  "black" = "4003",
  "asian" = "5",
  "other" = "6"
)

d$ethnicity <- fct_relevel(d$ethnicity, "white") # white as reference

summary(d$ethnicity)

## sex

d$sex <- as.factor(d$sex)

levels(d$sex) <- c("female", "male")

summary(d$sex)

## insomnia

d$insomnia_sr <- as.factor(d$insomnia_sr)

levels(d$insomnia_sr) <- c("prefer not answer", "never", "sometimes", "usually")

d$insomnia_sr <- fct_relevel(d$insomnia_sr, "never") # never as reference category

## chronotype

d$chronotype <- as.factor(d$chronotype)

levels(d$chronotype) <- c(
  "prefer not answer", "dont know", "morning",
  "more_morning", "more_evening", "evening"
)

d$chronotype <- fct_relevel(d$chronotype, "morning") # morning as reference

summary(d$chronotype)

## depression

d$freq_depressed_twoweeks <- as.factor(d$freq_depressed_twoweeks)

levels(d$freq_depressed_twoweeks) <- c(
  "prefer not answer", "dont know",
  "not at all", "several days", "more than half",
  "nearly every day"
)

d$freq_depressed_twoweeks <- fct_relevel(
  d$freq_depressed_twoweeks,
  "not at all"
) # not at all as ref category

summary(d$freq_depressed_twoweeks)

## Income

d$avg_total_household_income <- as.factor(d$avg_total_household_income)

levels(d$avg_total_household_income) <- c(
  "prefer not answer", "dont know",
  "<18", "18-30", "31-50", "52-100", ">100"
)

d$avg_total_household_income <- fct_relevel(
  d$avg_total_household_income,
  "31-50"
) # 31-50 as ref category
summary(d$avg_total_household_income)

## employment

d$employment <- as.factor(d$employment_1)

levels(d$employment) <- c(
  "none of above", "prefer not answer",
  "paid employment", "retired", "caring",
  "sick or disabled", "unemployed", "volunteer",
  "student"
)


d$employment <- fct_relevel(d$employment, "paid employment")

summary(d$employment)

# remove redundant variables

d <- d |> select(-c(employment_1:employment_7))

#### Trim dataset for selected sample ####

# Accelerometry data available

d2 <- d |> filter(!is.na(accel_data_available)) # 502359 -> 103661

# apply UKB data cleaning thresholds to accelerometry data

d2 <- d2 |> filter(is.na(d2$accel_data_problem)) # 98967

d2 <- d2 |> filter(accel_good_wear_time == 1) # 94497

d2 <- d2 |> filter(accel_good_calibration == 1) # 94494

d2 <- d2 |> filter(accel_calibrated_own_data == 1) # 94340

## must have GGIR data available

# Read in accelerometry data

a <- read_csv(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_daysumMM_output.csv"))

# filter to only those with GGIR data

d2 <- d2 |> filter(eid %in% a$eid) # 89626

### Merge in actigraphy data

# merge

d3 <- d2 |> left_join(a, by = "eid")

## At least 2 measurements of sleep duration

d3 <- d3 |>
  group_by(eid) |>
  mutate(n_valid = sum(!is.na(dur_spt_sleep_min))) |>
  ungroup() |>
  filter(n_valid >= 2)

### Create accelerometry time use variables

# total time awake during sleep window
d3$awake_sleep <-
  d3$dur_spt_wake_IN_min +
  d3$dur_spt_wake_LIG_min +
  d3$dur_spt_wake_MOD_min +
  d3$dur_spt_wake_VIG_min

# total wear time
d3$mins_worn <-
  d3$dur_spt_sleep_min +
  d3$dur_day_total_IN_min +
  d3$dur_day_total_LIG_min +
  d3$dur_day_total_MOD_min +
  d3$dur_day_total_VIG_min +
  d3$awake_sleep

# Variables normalised to 1440 relative to their proportion of total wear time

d3$sleep_n <-
  (d3$dur_spt_sleep_min / d3$mins_worn) * mins_in_day

d3$inactive_n <-
  ((d3$dur_day_total_IN_min + d3$dur_spt_wake_IN_min) / d3$mins_worn) * mins_in_day

d3$light_n <-
  ((d3$dur_day_total_LIG_min + d3$dur_spt_wake_LIG_min) / d3$mins_worn) * mins_in_day

d3$moderate_n <-
  ((d3$dur_day_total_MOD_min + d3$dur_spt_wake_MOD_min) / d3$mins_worn) * mins_in_day

d3$vigorous_n <-
  ((d3$dur_day_total_VIG_min + d3$dur_spt_wake_VIG_min) / d3$mins_worn) * mins_in_day

d3$mvpa_n <- d3$moderate_n + d3$vigorous_n

### Estimate average time spent in each of sleep, inactivity, moderate to vigorous activity,
### and light physical activity

##

# fit mixed model accounting for difference between days and also daylight savings
# crossover

d3$weekday <-
  factor(d3$weekday,
    levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
  )

m1 <-
  lmer(sleep_n ~ weekday + accel_daylight_savings_crossover + (1 | eid), data = d3)

summary(m1)

# Estimate average SRI for each participant, standardising for/marginalising over
# day of week.

nd <- expand_grid(
  weekday = levels(d3$weekday),
  accel_daylight_savings_crossover = 0,
  eid = unique(d3$eid)
)

nd$pred <- as_tibble(predict(m1, newdata = nd))$value

nd <- nd |>
  group_by(eid) |>
  summarise(avg_sleep = mean(pred))

## Repeat for inactivity

m2 <- lmer(inactive_n ~ weekday + accel_daylight_savings_crossover + (1 | eid), data = d3)

summary(m2)

# standardise

nd2 <- expand_grid(
  weekday = levels(d3$weekday),
  accel_daylight_savings_crossover = 0,
  eid = unique(d3$eid)
)

nd2$pred <- as_tibble(predict(m2, newdata = nd2))$value

nd2 <- nd2 |>
  group_by(eid) |>
  summarise(avg_inactivity = mean(pred))

## Repeat for light activity

m3 <- lmer(light_n ~ weekday + accel_daylight_savings_crossover + (1 | eid), data = d3)

summary(m3)

# standardise

nd3 <- expand_grid(
  weekday = levels(d3$weekday),
  accel_daylight_savings_crossover = 0,
  eid = unique(d3$eid)
)

nd3$pred <- as_tibble(predict(m3, newdata = nd3))$value

nd3 <- nd3 |>
  group_by(eid) |>
  summarise(avg_light = mean(pred))

## Repeat for MVPA

m4 <- lmer(mvpa_n ~ weekday + accel_daylight_savings_crossover + (1 | eid), data = d3)

summary(m4)

# standardise

nd4 <- expand_grid(
  weekday = levels(d3$weekday),
  accel_daylight_savings_crossover = 0,
  eid = unique(d3$eid)
)

nd4$pred <- as_tibble(predict(m4, newdata = nd4))$value

nd4 <- nd4 |>
  group_by(eid) |>
  summarise(avg_mvpa = mean(pred))

## Add standardised, averaged time use variables back into main data

nd <- left_join(nd, nd2, by = "eid") |>
  left_join(nd3, by = "eid") |>
  left_join(nd4, by = "eid")

d2 <- left_join(d2, nd, by = "eid")

# remove those with fewer than 2 sleep values

d2 <- d2 |> filter(!is.na(avg_sleep)) # 89545

## remove exclusionary neurological conditions

d2 <- d2 |> filter(neurological_exclude_bl == 0) # 88660

### Calculate time to dementia

# baseline (accelerometry) date

date <- a |>
  select(eid, calendar_date) |>
  filter(!duplicated(eid))

d2 <- d2 |> left_join(date, by = "eid")

d2 <-
  d2 |>
  mutate(calendar_date = lubridate::parse_date_time(calendar_date, orders = "ymd"))

# Date of first dementia diagnosis

demtime <- d |>
  mutate(
    date_all_cause_dementia = as_date(date_all_cause_dementia),
    dem_date_pc = as_date(dem_date_pc),
    time_dif = date_all_cause_dementia - dem_date_pc
  ) |>
  select(eid, date_all_cause_dementia, dem_date_pc, time_dif) |>
  pivot_longer(c(date_all_cause_dementia, dem_date_pc),
    values_to = "date_acdem2"
  ) |>
  arrange(eid, date_acdem2) |>
  group_by(eid) |>
  mutate(num = 1:n()) |>
  ungroup() |>
  filter(num == 1) |>
  select(eid, date_acdem2)

d2 <- left_join(d2, demtime, by = "eid")

# Update source of dementia cases

d2 <- d2 |>
  mutate(
    source_all_cause_dementia_report =
      ifelse(!date_acdem2 == date_all_cause_dementia,
        "GP", source_all_cause_dementia_report
      )
  ) |>
  mutate(
    source_all_cause_dementia_report =
      ifelse(!is.na(date_acdem2) & is.na(date_all_cause_dementia),
        "GP", source_all_cause_dementia_report
      )
  )

# dementia diagnosis

d2$dem <- ifelse(!is.na(d2$date_acdem2), 1, 0)

# calculate time to dementia

d2 <- d2 |>
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, calendar_date),
    !is.na(date_ltfu) ~ difftime(date_ltfu, calendar_date),
    TRUE ~ difftime("2022-02-01", calendar_date)
  )) |>
  mutate(time_to_dem = as.integer(time_to_dem))

## remove those with prevalent dementia

d2 <- d2 |> filter(time_to_dem > 0) ## 88654

#### Save created dataset ####

write_csv(d2, file.path(data_dir, "24hr_behaviours.csv"))

