imps <- list(seq_len(m))
names(imps) <- c("imp")

# comp_vars <- c("avg_sleep", "avg_mvpa", "avg_light", "avg_inactivity")
# sub_pairs <-
#   expand.grid(
#     from_var = comp_vars,
#     to_var = comp_vars,
#     stringsAsFactors = FALSE
#   ) |>
#   dplyr::filter(
#     !from_var == to_var &
#       (from_var == "avg_sleep" | to_var == "avg_sleep")
#   )
sub_pairs <- tidyr::tibble(
  from_var = c("avg_sleep"),
  to_var = c("avg_mvpa")
)

# sub_durations <- seq(from = -60, to = 60, by = 15)
sub_durations <- seq(from = 60, to = 60, by = 15)
sub_durations <- list(sub_durations[sub_durations != 0])
names(sub_durations) <- c("duration")

sub_targets <- list(
  tar_target(
    references,
    estimate_lmtp_reference(
      imp_wide,
      baseline_covars = c(
        "fruit_veg",
        "alc_freq",
        "sex",
        "retired",
        "shift",
        "apoe_e4",
        "highest_qual",
        "townsend_deprivation_index",
        "psych_meds",
        "ethnicity",
        "avg_total_household_income",
        "smok_status",
        "age_accel"
      )
    ),
    pattern = map(imp_wide)
  ),
  tar_map(
    values = sub_pairs,
    tar_map(
      values = sub_durations,
      tar_target(
        sub_df,
        apply_substitution(
          imp_wide,
          from_var,
          to_var,
          duration
        ),
        pattern = map(imp_wide),
        iteration = "list"
      ),
      tar_target(
        lmtp,
        estimate_lmtp_subs(
          imp_wide,
          sub_df,
          baseline_covars = c(
            "fruit_veg",
            "alc_freq",
            "sex",
            "retired",
            "shift",
            "apoe_e4",
            "highest_qual",
            "townsend_deprivation_index",
            "psych_meds",
            "ethnicity",
            "avg_total_household_income",
            "smok_status",
            "age_accel"
          )
        ),
        pattern = map(imp_wide, sub_df)
      ),
      tar_target(
        lmtp_combined,
        process_lmtp_imps(lmtp, references)
      )
    )
  )
)

process_lmtp_imps <- function(sub_imps, ref_imps) {
  subs_df <- bind_rows(lapply(sub_imps, function(sub_imp) {
    lmtp_result <- lmtp::tidy(sub_imp)
    lmtp_result[nrow(lmtp_result), ] |>
      select(estimate, std.error) |>
      mutate(estimate = 1 - estimate) |>
      rename(
        sub_risk = estimate,
        sub_err = std.error
      )
  }))

  refs_df <- bind_rows(lapply(ref_imps, function(ref_imp) {
    lmtp_result <- lmtp::tidy(ref_imp)
    lmtp_result[nrow(lmtp_result), ] |>
      select(estimate, std.error) |>
      mutate(estimate = 1 - estimate) |>
      rename(
        ref_risk = estimate,
        ref_err = std.error
      )
  }))

  M <- nrow(subs_df) # number of imputations

  combined <- bind_cols(subs_df, refs_df) |>
    mutate(
      log_rr = log(sub_risk / ref_risk),
      sub_var = (sub_err / sub_risk)^2, # delta–method
      ref_var = (ref_err / ref_risk)^2,
      se_rr = sqrt(sub_var + ref_var) # SE of log-RR in each imp.
    )

  pooled_log_rr <- mean(combined$log_rr)

  within_var <- mean(combined$se_rr^2) # Ŵ
  between_var <- sum((combined$log_rr - pooled_log_rr)^2) / (M - 1) # B
  total_var <- within_var + (1 + 1 / M) * between_var # T

  list(
    pooled_log_rr = pooled_log_rr,
    within_var = within_var,
    between_var = between_var,
    total_var = total_var,
    risk_ratio = exp(pooled_log_rr),
    refs_df = refs_df,
    subs_df = subs_df
  )
}
