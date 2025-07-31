get_primary_formula <- function(data) {
  knots_timegroup <- quantile(
    data[["timegroup"]],
    c(0.01, 0.5, 0.99)
  )
  knots_age <- quantile(
    data[["age_accel"]],
    c(0.1, 0.5, 0.9)
  )
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_fruit_veg <- quantile(
    data[["fruit_veg"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ poly(R1, 2) *
      rcs(timegroup, knots_timegroup) +
      poly(R2, 2) * rcs(timegroup, knots_timegroup) +
      poly(R3, 2) * rcs(timegroup, knots_timegroup) +
      poly(R1, 2) * as.numeric(avg_total_household_income) +
      poly(R2, 2) * as.numeric(avg_total_household_income) +
      poly(R3, 2) * as.numeric(avg_total_household_income) +
      poly(R1, 2) * sex +
      poly(R2, 2) * sex +
      poly(R3, 2) * sex +
      poly(R1, 2) * retired +
      poly(R2, 2) * retired +
      poly(R3, 2) * retired +
      poly(R1, 2) * smok_status +
      poly(R2, 2) * smok_status +
      poly(R3, 2) * smok_status +
      poly(R1, 2) * rcs(age_accel, knots_age) +
      poly(R2, 2) * rcs(age_accel, knots_age) +
      poly(R3, 2) * rcs(age_accel, knots_age) +
      rcs(timegroup, knots_timegroup) * rcs(age_accel, knots_age) +
      rcs(timegroup, knots_timegroup) * sex +
      rcs(timegroup, knots_timegroup) * as.numeric(apoe_e4) +
      highest_qual +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      shift +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity +
      rcs(age_accel, knots_age) * sex +
      rcs(age_accel, knots_age) * as.numeric(apoe_e4)
  )

  return(primary_formula)
}

get_s1_formula <- function(data) {
  knots_timegroup <- quantile(
    data[["timegroup"]],
    c(0.05, 0.275, 0.5, 0.725, 0.95)
  )
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_waso <- quantile(data[["avg_WASO"]], c(0.1, 0.5, 0.9))
  knots_fruit_veg <- quantile(data[["fruit_veg"]], c(0.1, 0.5, 0.9))

  s1_formula <- as.formula(
    ~ -1 +
      rcs(timegroup, knots_timegroup) +
      R1 +
      I(R1^2) +
      R2 +
      I(R2^2) +
      R3 +
      I(R3^2) +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      sex +
      retired +
      shift +
      apoe_e4 +
      highest_qual +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity +
      avg_total_household_income +
      smok_status +
      rcs(avg_WASO, knots_waso)
  )
  return(s1_formula)
}

get_s2_formula <- function(data) {
  knots_timegroup <- quantile(
    data[["timegroup"]],
    c(0.05, 0.275, 0.5, 0.725, 0.95)
  )
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_bmi <- quantile(data[["BMI"]], c(0.1, 0.5, 0.9))
  knots_bp <- quantile(data[["bp_syst_avg"]], c(0.1, 0.5, 0.9))
  knots_fruit_veg <- quantile(data[["fruit_veg"]], c(0.1, 0.5, 0.9))

  s2_formula <- as.formula(
    ~ -1 +
      rcs(timegroup, knots_timegroup) +
      R1 +
      I(R1^2) +
      R2 +
      I(R2^2) +
      R3 +
      I(R3^2) +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      sex +
      retired +
      shift +
      apoe_e4 +
      highest_qual +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
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
      rcs(BMI, knots_bmi) +
      rcs(bp_syst_avg, knots_bp)
  )
  return(s2_formula)
}

get_s3_formula <- function(data) {
  knots_timegroup <- quantile(
    data[["timegroup"]],
    c(0.05, 0.275, 0.5, 0.725, 0.95)
  )
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_fruit_veg <- quantile(
    data[["fruit_veg"]],
    c(0.1, 0.5, 0.9)
  )

  s3_formula <- as.formula(
    ~ -1 +
      rcs(timegroup, knots_timegroup) +
      (sex + retired + avg_total_household_income + smok_status) *
        (R1 + I(R1^2)) +
      (sex + retired + avg_total_household_income + smok_status) *
        (R2 + I(R2^2)) +
      (sex + retired + avg_total_household_income + smok_status) *
        (R3 + I(R3^2)) +
      shift +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      apoe_e4 +
      highest_qual +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity
  )
  return(s3_formula)
}

fit_models <- function(imp, timegroup_cuts, create_formula_fn) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  imp_survival <- survSplit(
    Surv(time = time_to_dem, event = dem) ~ .,
    data = imp,
    cut = timegroup_cuts,
    episode = "timegroup",
    end = "end",
    event = "dem",
    start = "start"
  )

  setDT(imp_survival)

  # start timegroup at 1
  imp_survival[, timegroup := timegroup - 1]

  # add indicators for death
  imp_survival[,
    death := fcase(
      death == 1 & end > time_to_death,
      1,
      death == 0 & end > time_to_death,
      NA,
      default = 0
    )
  ]

  imp_survival[, death := ifelse(dem == 1, NA, death)]

  # remove rows after death
  imp_survival[, sumdeath := cumsum(death), by = "id"]
  imp_survival <- imp_survival[sumdeath < 2 | is.na(sumdeath), ]
  imp_survival[, sumdeath := NULL]

  # model formula
  model_formula <- create_formula_fn(imp_survival)
  dem_model_formula <- update(model_formula, dem ~ .)
  death_model_formula <- update(model_formula, death ~ .)

  # fit models
  model_dem <- glm(
    dem_model_formula,
    data = imp_survival[death == 0 | is.na(death), ],
    family = binomial()
  )

  model_dem <- strip_glm(model_dem)

  model_death <- glm(
    death_model_formula,
    data = imp_survival[dem == 0, ],
    family = binomial()
  )

  model_death <- strip_glm(model_death)

  return(list(
    model_dem = model_dem,
    model_death = model_death
  ))
}

calc_substitution <- function(
  base_comp,
  imp_stacked_dt,
  model_dem,
  model_death,
  model_formula,
  substitution,
  timegroup
) {
  # The list of substitutions to be calculated in minutes
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  # Initialize a list to hold generated data.tables
  sub_comps_list <- vector("list", length(inc))

  # Loop over inc and create a sub_comps data table for each element
  for (i in seq_along(inc)) {
    # The list of compositions to be fed into the model after applying the substitutions
    sub_comps <- as.data.table(t(base_comp))
    setnames(
      sub_comps,
      c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
    )

    sub_comps[, (substitution[1]) := .SD[[substitution[1]]] + inc[i]]
    sub_comps[, (substitution[2]) := .SD[[substitution[2]]] - inc[i]]

    # Store the data.table into list
    sub_comps_list[[i]] <- sub_comps
  }

  # Combine all data.tables in the list
  sub_comps <- rbindlist(sub_comps_list)

  # The risk predicted by the model for each of these composition
  sub_risks <-
    rbindlist(lapply(
      seq_len(nrow(sub_comps)),
      function(i) {
        predict_composition_risk(
          acomp(sub_comps[i]),
          imp_stacked_dt,
          model_dem,
          model_death,
          model_formula,
          timegroup
        )
      }
    ))

  result <-
    setnames(
      data.table(offset = inc * mins_in_day, risks = sub_risks),
      c(
        "offset",
        paste0(
          substitution[2],
          "_",
          deparse(
            substitute(base_comp)
          )
        )
      )
    )

  return(result)
}
