prepare_dataset <- function(df, disease_file) {
  # sleep disorder prior to accelerometry?
  df$OSA_dx <-
    ifelse(
      df$OSA_sr == 1 | df$date_osa_dx <= df$date_accel,
      1,
      0
    )
  df$insomnia_dx <- ifelse(
    df$insomnia_sr == 1 |
      df$date_insomnia_dx <= df$date_accel,
    1,
    0
  )
  df$sleep_disorder_dx <- ifelse(
    df$sleep_disorder_sr == 1 |
      df$date_any_sleep_dx <= df$date_accel,
    1,
    0
  )

  # replace missing with zero
  df <-
    df |>
    mutate(across(
      c(OSA_dx, insomnia_dx, sleep_disorder_dx),
      ~ if_else(is.na(.), 0, .)
    ))

  ### Prepare confounding variables ###
  # set pack years to zero for non-smokers
  df$smok_pckyrs <- ifelse(df$smok_status == 0, 0, df$smok_pckyrs)

  # set alcohol prefer not to answer to missing
  df$alc_freq <- na_if(df$alc_freq, -3)

  # set do not know and prefer not to answer as missing for diet variables
  df <-
    df |>
    mutate(across(
      c("veg_cooked", "veg_raw", "fruit_fresh", "fruit_dry"),
      ~ na_if(., -1)
    ))
  df <-
    df |>
    mutate(across(
      c("veg_cooked", "veg_raw", "fruit_fresh", "fruit_dry"),
      ~ na_if(., -3)
    ))

  # set "less than one" to 0 for diet variables
  df <-
    df |>
    mutate(across(
      c("veg_cooked", "veg_raw", "fruit_fresh", "fruit_dry"),
      ~ if_else(. == -10, 0, .)
    ))

  # create total fruit and veg variable
  df$fruit_veg <- rowSums(
    df[, c(c("veg_cooked", "veg_raw", "fruit_fresh", "fruit_dry"))]
  )

  # set reference categories for factors
  df$diagnosed_diabetes <- as.factor(df$diagnosed_diabetes)
  levels(df$diagnosed_diabetes) <-
    c("prefer not answer", "dont know", "no", "yes")
  df$smok_status <- as.factor(df$smok_status)
  levels(df$smok_status) <-
    c("prefer not answer", "never", "former", "current")
  df$apoe_e4 <- as.factor(df$apoe_e4)
  df$highest_qual <- fct_relevel(df$highest_qual, "Grad")
  df$ethnicity <- fct_relevel(df$ethnicity, "white")
  df$insomnia_scale_sr <- fct_relevel(df$insomnia_scale_sr, "never")
  df$chronotype <- fct_relevel(df$chronotype, "morning")
  df$freq_depressed_twoweeks <- fct_relevel(
    df$freq_depressed_twoweeks,
    "not at all"
  )
  df$avg_total_household_income <- fct_relevel(
    df$avg_total_household_income,
    "31-50"
  )
  df$sick_disabled <- ifelse(df$employment == "sick or disabled", 1, 0)
  df$retired <- ifelse(df$employment == "retired", 1, 0)
  df$overall_health_rating <- as.factor(df$overall_health_rating)
  levels(df$overall_health_rating) <-
    c("prefer not answer", "dont know", "excellent", "good", "fair", "poor")
  df$shift <- ifelse(
    df$job_night_shift %in% c(3, 4) | df$job_shift_work %in% c(3, 4),
    1,
    0
  )

  # mark prefer not answer as missing
  df <- df |>
    mutate(across(where(is.factor), as.character)) |>
    mutate(across(where(is.character), ~ na_if(., "prefer not answer"))) |>
    mutate(across(where(is.character), ~ na_if(., "dont know"))) |>
    mutate(across(where(is.character), ~ na_if(., ""))) |>
    mutate(across(where(is.character), as.factor)) |>
    mutate(across(where(is.factor), fct_drop))

  # Update age variable to age at accelerometry study
  df$age_accel <-
    df$age_assessment +
    ((as.Date(df$date_accel) - as.Date(df$date_baseline)) / 365)

  df$age_accel <- as.numeric(df$age_accel)

  ### Create isometric log ratio coordinates

  # close time use variables

  dem_comp <-
    acomp(data.frame(
      df$avg_sleep,
      df$avg_inactivity,
      df$avg_light,
      df$avg_mvpa
    ))

  dem_base_ilr <-
    ilr(dem_comp, V = v) |>
    setNames(c("R1", "R2", "R3"))

  ### Add prevalent disease variables

  prev <- fread(disease_file, stringsAsFactors = TRUE)

  # add in date of actigraphy and number of cancers

  adate <- df |> select(eid, date_accel, num_sr_cancers)

  prev <- left_join(prev, adate, by = "eid")

  prev <- as_tibble(prev)

  # create prevalent illness variables

  prev <- prev |>
    mutate(across(
      starts_with("sr_"),
      ~ ifelse(
        is.na(.),
        0,
        .
      )
    ))

  prev$date_accel <- as_date(prev$date_accel)

  prev <- prev |>
    mutate(
      prev_diabetes = case_when(
        sr_diabetes == 1 ~ 1,
        sr_diabetes == 0 & date_diabetes < date_accel ~ 1,
        TRUE ~ 0
      )
    ) |>
    mutate(
      prev_cancer = case_when(
        num_sr_cancers > 0 ~ 1,
        num_sr_cancers == 0 & date_neoplasm < date_accel ~ 1,
        TRUE ~ 0
      )
    ) |>
    mutate(
      prev_mental_disorder = case_when(
        sr_mental_disorder == 1 ~ 1,
        sr_mental_disorder == 0 & date_mental_disorder < date_accel ~ 1,
        TRUE ~ 0
      )
    ) |>
    mutate(
      prev_nervous_system = case_when(
        sr_nervous_system == 1 ~ 1,
        sr_nervous_system == 0 & date_nervous_system_disorder < date_accel ~ 1,
        TRUE ~ 0
      )
    ) |>
    mutate(
      prev_cvd = case_when(
        sr_cvd == 1 ~ 1,
        sr_cvd == 0 & date_CVD < date_accel ~ 1,
        TRUE ~ 0
      )
    )

  prev <- prev |> select(eid, starts_with("prev_"))

  df <- left_join(df, prev, by = "eid")

  ### Select model data

  dem_model_data <- select(
    df,
    eid,
    fruit_veg,
    alc_freq,
    avg_sleep,
    avg_inactivity,
    avg_light,
    avg_mvpa,
    dem,
    death,
    time_to_dem,
    time_to_death,
    avg_WASO,
    bp_syst_avg,
    age_accel,
    retired,
    shift,
    townsend_deprivation_index,
    sex,
    psych_meds,
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
    sick_disabled,
    prev_diabetes,
    prev_cancer,
    prev_mental_disorder,
    prev_nervous_system,
    prev_cvd,
    bp_med,
    BMI,
    bp_syst_avg
  ) |>
    cbind(dem_base_ilr)

  return(dem_model_data)
}

impute_data <- function(df, m, maxit) {
  df <- ordinal_to_numeric(df)
  predmat <- quickpred(
    df,
    mincor = 0,
    exclude = c(
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa",
      "eid",
      "id"
    )
  )
  imp <- mice(
    df,
    m = m,
    predictorMatrix = predmat,
    maxit = maxit
  )
  print(imp$loggedEvents)
  imp <- complete(imp)
  imp <- back_to_factor(imp)
  return(imp)
}

back_to_factor <- function(df) {
  df$avg_total_household_income <- as.factor(
    df$avg_total_household_income
  )
  levels(df$avg_total_household_income) <- c(
    "<18",
    "18-30",
    "31-50",
    "52-100",
    ">100"
  )
  df$chronotype <- as.factor(df$chronotype)
  df$apoe_e4 <- as.factor(df$apoe_e4)
  df$highest_qual <- as.factor(df$highest_qual)
  df$smok_status <- as.factor(df$smok_status)
  levels(df$smok_status) <- c("current", "former", "never")
  return(as.data.table(df))
}

widen_data <- function(df, timegroup_cuts) {
  num_cuts <- length(timegroup_cuts)

  # Wide format
  lapply(
    df,
    function(.x) {
      .x$censoring_1 <- 1
      .x$dem_1 <- ifelse(
        .x$dem == 1 & .x$time_to_dem <= timegroup_cuts[2],
        1,
        0
      )
      .x$death_1 <- case_when(
        .x$death == 1 & .x$time_to_death < timegroup_cuts[2] ~ 1,
        .x$dem_1 == 1 ~ NA,
        .default = 0
      )

      for (i in 2:(num_cuts - 1)) {
        # Censoring
        .x[[paste0("censoring_", i)]] <- case_when(
          .x$dem == 0 & .x$death == 0 & .x$time_to_dem <= timegroup_cuts[i] ~ 0,
          .x[[paste0("dem_", i - 1)]] == 1 ~ NA,
          .x[[paste0("death_", i - 1)]] == 1 ~ NA,
          .default = 1
        )
        # Dem
        .x[[paste0("dem_", i)]] <- case_when(
          .x$dem == 1 & .x$time_to_dem <= timegroup_cuts[i + 1] ~ 1,
          .x[[paste0("death_", i - 1)]] == 1 ~ 0,
          .x[[paste0("censoring_", i)]] == 0 ~ NA,
          .default = 0
        )
        # Death
        .x[[paste0("death_", i)]] <- case_when(
          .x$death == 1 &
            .x$time_to_death <= timegroup_cuts[i + 1] &
            .x[[paste0("dem_", i)]] == 0 ~
            1,
          .x[[paste0("dem_", i)]] == 1 ~ NA,
          .x[[paste0("censoring_", i)]] == 0 ~ NA,
          .default = 0
        )
      }

      return(.x)
    }
  )
}
