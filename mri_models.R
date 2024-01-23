source("ideal_comp.R")
source("utils.R")

get_mri_formula <- function(outcome_var, model_data) {
  knots_age_mri_str <-
    paste(quantile(model_data[["age_mri"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_deprivation_str <-
    paste(quantile(model_data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_tfmri_str <-
    paste(quantile(model_data[["mean_tfmri_headmot"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_lat_bpos_str <-
    paste(quantile(model_data[["scan_lat_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_trans_bpos_str <-
    paste(quantile(model_data[["scan_trans_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")
  knots_long_bpos_str <-
    paste(quantile(model_data[["scan_long_bpos"]], c(0.1, 0.5, 0.9)), collapse = ", ")

  model_formula <- as.formula(paste(outcome_var, " ~
      pol(R1, 2) + pol(R2, 2) + pol(R3, 2) +
      rcs(age_mri, c(", knots_age_mri_str, ")) +
      retired +
      rcs(townsend_deprivation_index, c(", knots_deprivation_str, ")) +
      sex +
      antidepressant_med +
      antipsychotic_med +
      insomnia_med +
      ethnicity +
      avg_total_household_income +
      highest_qual +
      apoe_e4 +
      smok_status +
      rcs(mean_tfmri_headmot, c(", knots_tfmri_str, ")) +
      rcs(scan_lat_bpos, c(", knots_lat_bpos_str, ")) +
      rcs(scan_trans_bpos, c(", knots_trans_bpos_str, ")) +
      rcs(scan_long_bpos, c(", knots_long_bpos_str, ")) +
      scan_tabpos +
      assessment_centre_mri1"))

  return(model_formula)
}

calc_mri_substitution <- function(base_comp, imp, model, substitution) {
  # The list of substitutions to be calculated in minutes
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  # Initialize a list to hold generated data frames
  sub_comps_list <- vector("list", length(inc))

  # Loop over inc and create a sub_comps data frame for each element
  for (i in seq_along(inc)) {
    # The list of compositions to be fed into the model after applying the substitutions
    sub_comps <- as.data.frame(t(base_comp))
    colnames(sub_comps) <- c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")

    sub_comps[[substitution[1]]] <- sub_comps[[substitution[1]]] + inc[i]
    sub_comps[[substitution[2]]] <- sub_comps[[substitution[2]]] - inc[i]

    # Store the data frame into list
    sub_comps_list[[i]] <- sub_comps
  }

  # Combine all data frames in the list
  sub_comps <- do.call(rbind, sub_comps_list)

  # The risk predicted by the model for each of these composition

  sub_risks <-
    do.call(rbind,
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
  colnames(result) <- c("offset", paste0(substitution[2], "_", deparse(substitute(base_comp))))

  return(result)
}

predict_mri_substitutions <- function(outcome_var, model_data, short_sleep_geo_mean, avg_sleep_geo_mean) {
  model_formula <- get_mri_formula(outcome_var, model_data)

  model <- ols(model_formula, data = model_data)
  print(summary(model))

  sleep_inactive <- c("avg_sleep", "avg_inactivity")
  sleep_light <- c("avg_sleep", "avg_light")
  sleep_mvpa <- c("avg_sleep", "avg_mvpa")

  avg_sleep_inactive <- calc_mri_substitution(avg_sleep_geo_mean, model_data, model, sleep_inactive)
  avg_sleep_light <- calc_mri_substitution(avg_sleep_geo_mean, model_data, model, sleep_light)
  avg_sleep_mvpa <- calc_mri_substitution(avg_sleep_geo_mean, model_data, model, sleep_mvpa)
  short_sleep_inactive <- calc_mri_substitution(short_sleep_geo_mean, model_data, model, sleep_inactive)
  short_sleep_light <- calc_mri_substitution(short_sleep_geo_mean, model_data, model, sleep_light)
  short_sleep_mvpa <- calc_mri_substitution(short_sleep_geo_mean, model_data, model, sleep_mvpa)

  full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset")

  return(full_df)
}

predict_mri_outcome <- function(outcome_var, model_data, best_and_worst) {
  model_formula <- get_mri_formula(outcome_var, model_data)

  model <- ols(model_formula, data = model_data)
  print(summary(model))

  # Update the local copy with 'best' values and predict outcomes
  model_data[, "R1"] <- best_and_worst$best[[1]]
  model_data[, "R2"] <- best_and_worst$best[[2]]
  model_data[, "R3"] <- best_and_worst$best[[3]]
  predictions_best <- predict(model, newdata = model_data)
  best <- mean(predictions_best)

  # Update the local copy with 'worst' values and predict outcomes
  model_data[, "R1"] <- best_and_worst$worst[[1]]
  model_data[, "R2"] <- best_and_worst$worst[[2]]
  model_data[, "R3"] <- best_and_worst$worst[[3]]
  predictions_worst <- predict(model, newdata = model_data)
  worst <- mean(predictions_worst)

  # Update the local copy with 'most_common' values and predict outcomes
  model_data[, "R1"] <- best_and_worst$most_common[[1]]
  model_data[, "R2"] <- best_and_worst$most_common[[2]]
  model_data[, "R3"] <- best_and_worst$most_common[[3]]
  predictions_common <- predict(model, newdata = model_data)
  common <- mean(predictions_common)

  return(data.frame(
    value = c(worst, common, best)
  ))
}

predict_mri_results <- function(model_data, best_and_worst) {
  results <- lapply(c("tbv", "wmv", "gmv", "hip", "log_wmh"), function(outcome) {
    df <- predict_mri_outcome(outcome, model_data, best_and_worst)
    return(df)
  })

  final_df <- do.call(rbind, results)
  return(final_df)
}

predict_all_substitutions <- function(model_data, short_sleep_geo_mean, avg_sleep_geo_mean) {
  results <- lapply(c("tbv", "wmv", "gmv", "hip", "log_wmh"), function(outcome) {
    df <- predict_mri_substitutions(outcome, model_data, short_sleep_geo_mean, avg_sleep_geo_mean)
    return(df)
  })

  final_df <- do.call(rbind, results)
  return(final_df)
}
