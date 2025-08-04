create_data <- function(
  core_file,
  demdeath_file,
  snp_file,
  diet_file,
  accel_file,
  sleep_file,
  sri_file
) {
  # read UKB data

  d <- fread(core_file)

  # up to date event variables

  latest <- fread(demdeath_file)

  # rename variables
  latest <- latest |>
    rename(
      date_of_death = `40000-0.0`,
      date_all_cause_dementia = `42018-0.0`
    )
  d <- d |>
    select(
      -date_of_death,
      -date_all_cause_dementia
    )
  d <- latest |> left_join(d, by = "eid")

  # Add in SNPs #

  snps <- fread(snp_file)

  d <- d |> left_join(snps, by = "eid")

  # Add in diet and alcohol variables #

  alc_diet <- fread(diet_file)

  alc_diet <-
    alc_diet |>
    set_names(
      c(
        "eid",
        "alc_freq",
        "veg_cooked",
        "veg_raw",
        "fruit_fresh",
        "fruit_dry"
      )
    )

  d <- left_join(d, alc_diet, by = "eid")

  #### Clean confounder variables ####

  # remove negative values from sleep duration

  d <-
    d |>
    mutate(
      sleep_duration_sr = ifelse(
        sleep_duration_sr < 0,
        NA_integer_,
        sleep_duration_sr
      )
    )

  ### Set categorical variables to factor

  d$apoe_e4 <- as.factor(d$apoe_e4)

  d$bp_med <- as.factor(d$bp_med)

  d$any_cvd <- as.factor(d$any_cvd)

  # Highest qualification

  quals <- d |>
    select(eid, contains("qualifications")) |>
    pivot_longer(-eid) |>
    mutate(
      highest_qual = case_when(
        value == 1 ~ 4,
        value == 6 ~ 4,
        value == 5 ~ 3,
        value == 2 ~ 2,
        value == 3 ~ 1,
        value == 4 ~ 1,
        TRUE ~ as.numeric(value)
      )
    ) |>
    group_by(eid) |>
    summarise(highest_qual = max(highest_qual, na.rm = T)) |>
    mutate(
      highest_qual = ifelse(
        highest_qual == -Inf,
        NA_real_,
        highest_qual
      )
    )

  d <- left_join(d, quals, by = "eid")

  d$highest_qual <- as.factor(d$highest_qual)

  d$highest_qual <- relevel(d$highest_qual, "4")

  levels(d$highest_qual) <- c(
    "Grad",
    "Other",
    "prefer not answer",
    "O",
    "A",
    "NVQ"
  )

  # remove redundant variables

  d <- d |> select(-starts_with("qualification"))

  d$antidepressant_med <- as.numeric(d$antidepressant_med)
  d$antipsychotic_med <- as.numeric(d$antipsychotic_med)
  d$insomnia_med <- as.numeric(d$insomnia_med)
  d$psych_meds <- with(
    d,
    as.factor(
      (antidepressant_med == 1 | antipsychotic_med == 1 | insomnia_med == 1) * 1
    )
  )

  ## Race

  d$ethnicity <- as.factor(d$ethnicity)

  # collapse categories

  d$ethnicity <- fct_recode(
    d$ethnicity,
    "prefer not answer" = "-3",
    "other" = "-1",
    "white" = "1",
    "white" = "1001",
    "white" = "2001",
    "white" = "3001",
    "white" = "4001",
    "other" = "2",
    "other" = "1002",
    "other" = "2002",
    "other" = "3002",
    "other" = "4002",
    "other" = "3",
    "other" = "1003",
    "other" = "2003",
    "other" = "3003",
    "other" = "4",
    "other" = "2004",
    "other" = "3004",
    "other" = "4003",
    "other" = "5",
    "other" = "6"
  )

  d$ethnicity <- fct_relevel(d$ethnicity, "white")

  ## sex

  d$sex <- as.factor(d$sex)

  levels(d$sex) <- c("female", "male")

  ## insomnia

  d$insomnia_scale_sr <- as.factor(d$insomnia_sr)

  levels(d$insomnia_scale_sr) <- c(
    "prefer not answer",
    "never",
    "sometimes",
    "usually"
  )

  d$insomnia_scale_sr <- fct_relevel(d$insomnia_scale_sr, "never")

  ## chronotype

  d$chronotype <- as.factor(d$chronotype)

  levels(d$chronotype) <- c(
    "prefer not answer",
    "dont know",
    "morning",
    "more_morning",
    "more_evening",
    "evening"
  )

  d$chronotype <- fct_relevel(d$chronotype, "morning")

  ## depression

  d$freq_depressed_twoweeks <- as.factor(d$freq_depressed_twoweeks)

  levels(d$freq_depressed_twoweeks) <- c(
    "prefer not answer",
    "dont know",
    "not at all",
    "several days",
    "more than half",
    "nearly every day"
  )

  d$freq_depressed_twoweeks <- fct_relevel(
    d$freq_depressed_twoweeks,
    "not at all"
  ) # not at all as ref category

  ## Income

  d$avg_total_household_income <- as.factor(
    d$avg_total_household_income
  )

  levels(d$avg_total_household_income) <- c(
    "prefer not answer",
    "dont know",
    "<18",
    "18-30",
    "31-50",
    "52-100",
    ">100"
  )

  d$avg_total_household_income <- fct_relevel(
    d$avg_total_household_income,
    "31-50"
  ) # 31-50 as ref category

  ## employment

  d$employment <- as.factor(d$employment_1)

  levels(d$employment) <- c(
    "none of above",
    "prefer not answer",
    "paid employment",
    "retired",
    "caring",
    "sick or disabled",
    "unemployed",
    "volunteer",
    "student"
  )

  d$employment <- fct_relevel(d$employment, "paid employment")

  # remove redundant variables

  d <- d |> select(-c(employment_1:employment_7))

  #### Trim dataset for selected sample ####

  # tracking sample size

  sample_size_info <- list()

  # Accelerometry data available

  sample_size_info$total_cohort <- nrow(d)

  d2 <- d |> filter(!is.na(accel_data_available))

  sample_size_info$accel_available <- nrow(d2)

  # apply UKB data cleaning thresholds to accelerometry data

  d2 <- d2 |> filter(is.na(d2$accel_data_problem))

  sample_size_info$no_accel_data_problem <- nrow(d2)

  d2 <- d2 |> filter(accel_good_wear_time == 1)

  sample_size_info$good_wear_time <- nrow(d2)

  d2 <- d2 |> filter(accel_good_calibration == 1)

  d2 <- d2 |> filter(accel_calibrated_own_data == 1)

  sample_size_info$good_calibration <- nrow(d2)

  ## must have GGIR data available

  # Read in accelerometry data

  a <- fread(accel_file)

  # filter to only those with GGIR data

  d2 <- d2 |> filter(eid %in% a$eid)

  ### Merge in actigraphy data

  # merge

  d3 <- d2 |> left_join(a, by = "eid")

  ## At least 3 measurements of sleep duration

  # Set sleep to missing if cleaningcode indicates problem

  d3 <- d3 |>
    mutate(dur_spt_sleep_min = ifelse(cleaningcode == 1, dur_spt_sleep_min, NA))

  d3 <- d3 |>
    group_by(eid) |>
    mutate(n_valid = sum(!is.na(dur_spt_sleep_min))) |>
    ungroup() |>
    filter(n_valid >= 3)

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
    ((d3$dur_day_total_IN_min + d3$dur_spt_wake_IN_min) /
      d3$mins_worn) *
    mins_in_day

  d3$light_n <-
    ((d3$dur_day_total_LIG_min + d3$dur_spt_wake_LIG_min) /
      d3$mins_worn) *
    mins_in_day

  d3$moderate_n <-
    ((d3$dur_day_total_MOD_min + d3$dur_spt_wake_MOD_min) /
      d3$mins_worn) *
    mins_in_day

  d3$vigorous_n <-
    ((d3$dur_day_total_VIG_min + d3$dur_spt_wake_VIG_min) /
      d3$mins_worn) *
    mins_in_day

  d3$mvpa_n <- d3$moderate_n + d3$vigorous_n

  ### Estimate average time spent in each of sleep, inactivity,
  ### moderate to vigorous activity and light physical activity

  # fit mixed model accounting for difference between days
  # and daylight savings crossover

  d3$weekday <-
    factor(
      d3$weekday,
      levels = c(
        "Sunday",
        "Monday",
        "Tuesday",
        "Wednesday",
        "Thursday",
        "Friday",
        "Saturday"
      )
    )

  m1 <-
    lmer(
      sleep_n ~ weekday + accel_daylight_savings_crossover + (1 | eid),
      data = d3
    )

  # Estimate average SRI for each participant, standardising over
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

  m2 <- lmer(
    inactive_n ~ weekday + accel_daylight_savings_crossover + (1 | eid),
    data = d3
  )

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

  m3 <- lmer(
    light_n ~ weekday + accel_daylight_savings_crossover + (1 | eid),
    data = d3
  )

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

  m4 <- lmer(
    mvpa_n ~ weekday + accel_daylight_savings_crossover + (1 | eid),
    data = d3
  )

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

  # remove those with insufficient sleep values/failed GGIR quality checks

  d2 <- d2 |> filter(!is.na(avg_sleep))

  sample_size_info$GGIR_checks <- nrow(d2)

  # Remove observations outside of 0.1st and 99.9th percentiles
  d2 <- d2 |>
    filter(
      avg_sleep > quantile(avg_sleep, 0.001) &
        avg_sleep < quantile(avg_sleep, 0.999),
      avg_mvpa > quantile(avg_mvpa, 0.001) &
        avg_mvpa < quantile(avg_mvpa, 0.999),
      avg_light > quantile(avg_light, 0.001) &
        avg_light < quantile(avg_light, 0.999),
      avg_inactivity > quantile(avg_inactivity, 0.001) &
        avg_inactivity < quantile(avg_inactivity, 0.999)
    )

  sample_size_info$quantile_exclusion <- nrow(d2)

  ## remove exclusionary neurological conditions

  d2 <- d2 |> filter(neurological_exclude_bl == 0)

  sample_size_info$Non_neurological_exclusion <- nrow(d2)

  ### Calculate time to dementia

  # baseline (accelerometry) date

  date <- a |>
    select(eid, calendar_date) |>
    filter(!duplicated(eid))

  d2 <- d2 |> left_join(date, by = "eid")

  d2 <-
    d2 |>
    mutate(
      calendar_date = lubridate::parse_date_time(
        calendar_date,
        orders = "ymd"
      )
    )

  # Date of first dementia diagnosis

  demtime <- d |>
    mutate(
      date_all_cause_dementia = as_date(date_all_cause_dementia),
      dem_date_pc = as_date(dem_date_pc),
      time_dif = date_all_cause_dementia - dem_date_pc
    ) |>
    select(eid, date_all_cause_dementia, dem_date_pc, time_dif) |>
    pivot_longer(
      c(date_all_cause_dementia, dem_date_pc),
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
      source_all_cause_dementia_report = ifelse(
        !date_acdem2 == date_all_cause_dementia,
        "GP",
        source_all_cause_dementia_report
      )
    ) |>
    mutate(
      source_all_cause_dementia_report = ifelse(
        !is.na(date_acdem2) & is.na(date_all_cause_dementia),
        "GP",
        source_all_cause_dementia_report
      )
    )

  # dementia diagnosis

  d2$dem <- ifelse(!is.na(d2$date_acdem2), 1, 0)

  ## Construct risk set
  # Death from other causes remain in risk set until end of FU

  d2 <- d2 |> rename("date_accel" = "calendar_date")

  d2 <- d2 |>
    mutate(
      competing = ifelse(
        !is.na(date_of_death) & is.na(date_acdem2),
        1,
        0
      )
    ) |>
    mutate(
      time_to_dem = case_when(
        dem == 1 ~ difftime(date_acdem2, date_accel),
        competing == 1 ~ difftime("2023-01-01", date_accel),
        TRUE ~ difftime("2023-01-01", date_accel)
      )
    ) |>
    mutate(time_to_dem = as.integer(time_to_dem))

  # create death and time to death variable

  d2$death <- ifelse(is.na(d2$date_of_death), 0, 1)

  d2$time_to_death <- ifelse(
    d2$death == 1,
    difftime(d2$date_of_death, d2$date_accel),
    difftime("2023-01-01", d2$date_accel)
  )

  ## remove those with prevalent dementia

  d2 <- d2 |> filter(time_to_dem > 0)

  sample_size_info$no_prevalent_dementia <- nrow(d2)

  # Sleep disorders data
  sleep_dis_df <- fread(sleep_file)

  # Merge dataframes
  d2 <- left_join(d2, sleep_dis_df, by = "eid")

  # rename and remove some variables
  d2 <- d2 |> select(-starts_with("dur_day_total_"))

  ## Add WASO and SRI to dataset ##

  waso_dat <- fread(sri_file)

  waso_dat <- select(waso_dat, eid, avg_WASO, avg_sri)

  d2 <- left_join(d2, waso_dat, by = "eid")

  #### Save created dataset ####

  return(list(df = d2, sample_size_info = sample_size_info))
}
