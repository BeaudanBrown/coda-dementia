get_mri_formula <- function(data) {
  knots_age_mri <- quantile(
    data[["age_mri"]],
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
  knots_tfmri <- quantile(
    data[["mean_tfmri_headmot"]],
    c(0.1, 0.5, 0.9)
  )
  knots_head_scale <- quantile(
    data[["head_scale"]],
    c(0.1, 0.5, 0.9)
  )
  knots_lat_bpos <- quantile(
    data[["scan_lat_bpos"]],
    c(0.1, 0.5, 0.9)
  )
  knots_trans_bpos <- quantile(
    data[["scan_trans_bpos"]],
    c(0.1, 0.5, 0.9)
  )
  knots_long_bpos <- quantile(
    data[["scan_long_bpos"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ poly(R1, 2) *
      as.numeric(avg_total_household_income) +
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
      poly(R1, 2) * rcs(age_mri, knots_age_mri) +
      poly(R2, 2) * rcs(age_mri, knots_age_mri) +
      poly(R3, 2) * rcs(age_mri, knots_age_mri) +
      as.numeric(apoe_e4) +
      highest_qual +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      shift +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity +
      rcs(age_mri, knots_age_mri) * sex +
      rcs(age_mri, knots_age_mri) * as.numeric(apoe_e4) +
      rcs(mean_tfmri_headmot, knots_tfmri) +
      rcs(head_scale, knots_head_scale) +
      rcs(scan_lat_bpos, knots_lat_bpos) +
      rcs(scan_trans_bpos, knots_trans_bpos) +
      scan_tabpos +
      assessment_centre_mri1
  )

  return(primary_formula)
}

calc_mri_substitution <- function(base_comp, imp, model, substitution) {
  # The list of substitutions to be calculated in minutes
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  # Initialize a list to hold generated data frames
  sub_comps_list <- vector("list", length(inc))

  # Loop over inc and create a sub_comps data frame for each element
  for (i in seq_along(inc)) {
    # The list of compositions to be fed into the model after
    # applying the substitutions
    sub_comps <- as.data.frame(t(base_comp))
    colnames(sub_comps) <-
      c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")

    sub_comps[[substitution[1]]] <- sub_comps[[substitution[1]]] + inc[i]
    sub_comps[[substitution[2]]] <- sub_comps[[substitution[2]]] - inc[i]

    # Store the data frame into list
    sub_comps_list[[i]] <- sub_comps
  }

  # Combine all data frames in the list
  sub_comps <- do.call(rbind, sub_comps_list)

  # The risk predicted by the model for each of these composition

  sub_risks <-
    do.call(
      rbind,
      lapply(
        seq_len(nrow(sub_comps)),
        function(i) {
          ilr <- ilr(sub_comps[i, ], V = v)
          imp$R1 <- ilr[1]
          imp$R2 <- ilr[2]
          imp$R3 <- ilr[3]
          prediction <- predict(model, newdata = imp)
          mean(prediction)
        }
      )
    )

  result <- data.frame(offset = inc * mins_in_day, risks = sub_risks)
  colnames(result) <- c(
    "offset",
    paste0(
      substitution[2],
      "_",
      deparse(substitute(base_comp))
    )
  )

  return(result)
}

predict_mri_substitutions <- function(
  outcome_var,
  model_data,
  short_sleep_geo_mean,
  avg_sleep_geo_mean
) {
  model_formula <- get_mri_formula(outcome_var, model_data)

  model <- ols(model_formula, data = model_data)
  print(summary(model))

  sleep_inactive <- c("avg_sleep", "avg_inactivity")
  sleep_light <- c("avg_sleep", "avg_light")
  sleep_mvpa <- c("avg_sleep", "avg_mvpa")

  avg_sleep_inactive <- calc_mri_substitution(
    avg_sleep_geo_mean,
    model_data,
    model,
    sleep_inactive
  )
  avg_sleep_light <- calc_mri_substitution(
    avg_sleep_geo_mean,
    model_data,
    model,
    sleep_light
  )
  avg_sleep_mvpa <- calc_mri_substitution(
    avg_sleep_geo_mean,
    model_data,
    model,
    sleep_mvpa
  )
  short_sleep_inactive <- calc_mri_substitution(
    short_sleep_geo_mean,
    model_data,
    model,
    sleep_inactive
  )
  short_sleep_light <- calc_mri_substitution(
    short_sleep_geo_mean,
    model_data,
    model,
    sleep_light
  )
  short_sleep_mvpa <- calc_mri_substitution(
    short_sleep_geo_mean,
    model_data,
    model,
    sleep_mvpa
  )

  full_df <- full_join(
    short_sleep_inactive,
    short_sleep_light,
    by = "offset"
  ) |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset")

  return(full_df)
}

predict_mri_outcome <- function(outcome_var, model_data, best_and_worst) {
  model_formula <- get_mri_formula(outcome_var, model_data)

  model <- ols(model_formula, data = model_data)

  do_prediction <- function(data, comp) {
    data[, "R1"] <- comp[[1]]
    data[, "R2"] <- comp[[2]]
    data[, "R3"] <- comp[[3]]
    predictions <- predict(model, newdata = data)
    mean_pred <- mean(predictions)
    return(mean_pred)
  }

  best <- do_prediction(model_data, best_and_worst$best)
  worst <- do_prediction(model_data, best_and_worst$worst)
  typical <- do_prediction(model_data, best_and_worst$typical)

  return(data.frame(
    value = c(worst, typical, best)
  ))
}

predict_mri_results <- function(model_data, best_and_worst) {
  results <- lapply(
    c("tbv", "wmv", "gmv", "hip", "log_wmh"),
    function(outcome) {
      df <- predict_mri_outcome(outcome, model_data, best_and_worst)
      return(df)
    }
  )

  final_df <- do.call(rbind, results)
  return(final_df)
}

predict_all_substitutions <- function(
  model_data,
  short_sleep_geo_mean,
  avg_sleep_geo_mean
) {
  results <- lapply(
    c("tbv", "wmv", "gmv", "hip", "log_wmh"),
    function(outcome) {
      df <- predict_mri_substitutions(
        outcome,
        model_data,
        short_sleep_geo_mean,
        avg_sleep_geo_mean
      )
      return(df)
    }
  )

  final_df <- do.call(rbind, results)
  return(final_df)
}
