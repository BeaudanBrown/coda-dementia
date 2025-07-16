get_primary_outcome_formula <- function(data) {
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_fruit_veg <- quantile(
    data[["fruit_veg"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ (R1 + R2 + R3) *
      (sex +
        retired +
        avg_total_household_income +
        highest_qual +
        smok_status) +
      I(R1^2) +
      I(R2^2) +
      I(R3^2) +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      shift +
      apoe_e4 +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity
  )

  return(primary_formula)
}

get_primary_treatment_formula <- function(data) {
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_fruit_veg <- quantile(
    data[["fruit_veg"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ (sex +
      retired +
      avg_total_household_income +
      highest_qual +
      smok_status)^2 +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      shift +
      apoe_e4 +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity
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

fit_model <- function(imp, timegroup_cuts, create_formula_fn) {
  imp_long <- survSplit(
    Surv(time = age_accel, event = dem, time2 = age_dem) ~ .,
    data = imp,
    cut = timegroup_cuts,
    episode = "timegroup",
    end = "age_end",
    event = "dem",
    start = "age_start"
  )

  # add indicators for death
  setDT(imp_long)

  imp_long$death <-
    ifelse(imp_long$death == 1 & imp_long$age_at_death < imp_long$age_end, 1, 0)

  # remove rows after death
  imp_long[, sumdeath := cumsum(death), by = "eid"]
  imp_long <- imp_long[sumdeath < 2, ]
  model_formula <- create_formula_fn(imp_long)

  # predictor and outcome matrices for fastglm
  Xdem <- model.matrix(
    model_formula,
    imp_long[death == 0, ]
  )
  Ydem <- as.matrix(imp_long[death == 0, ]$dem)

  Xdeath <- model.matrix(
    model_formula,
    imp_long[dem == 0, ]
  )
  Ydeath <- as.matrix(imp_long[dem == 0, ]$death)

  # fit models
  model_dem <- fastglm(
    x = Xdem,
    y = Ydem,
    family = binomial(),
    method = 3
  )
  model_death <- fastglm(
    x = Xdeath,
    y = Ydeath,
    family = binomial,
    method = 3
  )

  return(list(
    model_dem = model_dem,
    model_death = model_death,
    model_formula = model_formula
  ))
}

predict_composition_risk <- function(
  composition,
  stacked_data_table,
  model_dem,
  model_death,
  model_formula,
  timegroup
) {
  ilr <- ilr(composition, V = v)

  newx <- model.matrix(model_formula, stacked_data_table)

  newx[, "R1"] <- ilr[1]
  newx[, "I(R1^2)"] <- ilr[1]^2
  newx[, "R2"] <- ilr[2]
  newx[, "I(R2^2)"] <- ilr[2]^2
  newx[, "R3"] <- ilr[3]
  newx[, "I(R3^2)"] <- ilr[3]^2

  stacked_data_table[,
    haz_dem := predict(
      model_dem,
      newdata = newx,
      type = "response"
    )
  ]
  stacked_data_table[,
    haz_death := predict(
      model_death,
      newdata = newx,
      type = "response"
    )
  ]

  setkey(stacked_data_table, id, timegroup) # sort and set keys
  stacked_data_table[,
    risk := cumsum(
      haz_dem * cumprod((1 - lag(haz_dem, default = 0)) * (1 - haz_death))
    ),
    by = id
  ]

  risk <- stacked_data_table[timegroup == timegroup, .(mean = mean(risk))]

  return(risk)
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
