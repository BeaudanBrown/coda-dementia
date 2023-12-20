source("ideal_comp.R")

# sbp <- matrix(
#   c(
#     1, 1, -1, -1,
#     1, -1, 0, 0,
#     0, 0, 1, -1
#   ),
#   ncol = 4, byrow = TRUE
# )
# v <- gsi.buildilrBase(t(sbp))

# make_ilr <- function(comp) {
#   return(ilr(acomp(comp), V = v))
# }
# best_and_worst <- get_best_and_worst_comp()
# best_and_worst <- as.data.frame(apply(best_and_worst, 2, make_ilr))

# mri_df <-
#   fread(file.path(data_dir, "mri_data.csv"), stringsAsFactors = TRUE) |>
#   as_tibble()

# mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)

# ## merge MRI and main dataset

# boot_df <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

# mri_df <- mri_df |> left_join(boot_df, by = "eid")

# ## Remove participants without accel data

# mri_df <- mri_df |> filter(!is.na(age_accel))

# ## Remove participants with accelerometry after MRI

# mri_df <- mri_df |>
#   mutate(date_mri1 = as_date(date_mri1),
#          date_accel = as_date(date_accel)) |>
#   mutate(mri_accel_time_dif = as.integer(date_mri1 - date_accel)) |>
#   mutate(mri_before_accel = ifelse(mri_accel_time_dif < 0, 1, 0))

# mri_df <- mri_df |> filter(mri_before_accel == 0)

# # Rename age at MRI scan
# mri_df <- rename(mri_df, "age_mri" = "age_assessment_mri1")

get_mri_formula <- function(outcome_var, model_data) {
  knots_avg_sri_str <-
    paste(quantile(model_data[["avg_sri"]], c(0.1, 0.5, 0.9)), collapse = ", ")
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
      rcs(R1, 3) + rcs(R2, 3) + rcs(R3, 3) +
      rcs(avg_sri, c(", knots_avg_sri_str, ")) +
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

predict_mri_outcome <- function(outcome_var, model_data, best_and_worst) {
  model_formula <- get_mri_formula(outcome_var, model_data)

  model <- ols(model_formula,
    data = model_data
  )
  print(summary(model))

  ref_row <- as.data.frame(model_data[1, ])
  ref_row[c("R1", "R2", "R3")] <- best_and_worst$best
  best <- predict(model, newdata = ref_row)

  ref_row[c("R1", "R2", "R3")] <- best_and_worst$worst
  worst <- predict(model, newdata = ref_row)

  ref_row[c("R1", "R2", "R3")] <- best_and_worst$most_common
  common <- predict(model, newdata = ref_row)

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

# ggplot(final_df, aes(x = reference, y = value, color = outcome)) +
#   geom_point() +
#   facet_wrap(~outcome, scales = "free") +
#   xlab("") +
#   ylab("Value") +
#   theme_bw()
