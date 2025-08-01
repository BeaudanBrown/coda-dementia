train_model <- function(
  train_data,
  timegroup_cuts,
  create_formula_fn
) {
  # fit model
  train_data[, id := .I]
  model_formula <- create_formula_fn(train_data)
  models <- fit_models(train_data, timegroup_cuts, create_formula_fn)
  models
}

get_covar_data <- function(train_data) {
  # minimal data for predicting
  pred_data <- expand_grid(
    avg_total_household_income = unique(train_data$avg_total_household_income),
    sex = unique(train_data$sex),
    retired = unique(train_data$retired),
    smok_status = unique(train_data$smok_status),
    age_accel = unique(round(train_data$age_accel))
  )
  setDT(pred_data)
  # Add proportion in data for each row
  train_data[, age_round := round(age_accel)]
  prop <- vector("numeric", nrow(pred_data))
  for (i in seq_len(nrow(pred_data))) {
    prop[i] <- nrow(train_data[
      avg_total_household_income == pred_data[i, ]$avg_total_household_income &
        sex == pred_data[i, ]$sex &
        retired == pred_data[i, ]$retired &
        smok_status == pred_data[i, ]$smok_status &
        age_round == pred_data[i, ]$age_accel,
    ]) /
      nrow(train_data)
  }
  pred_data$prop <- prop
  pred_data <- pred_data[prop > 0, ]

  # add mean for other covars
  pred_data[, `:=`(
    apoe_e4 = Mode(train_data$apoe_e4),
    shift = Mode(train_data$shift),
    highest_qual = Mode(train_data$highest_qual),
    fruit_veg = mean(train_data$fruit_veg),
    alc_freq = Mode(train_data$alc_freq),
    townsend_deprivation_index = mean(train_data$townsend_deprivation_index),
    psych_meds = Mode(train_data$psych_meds),
    ethnicity = Mode(train_data$ethnicity)
  )]

  pred_data[, id := .I]
  pred_data
}

get_synth_risk <- function(pred_data, synth_comps, models, final_time) {
  setDT(pred_data)
  setDT(synth_comps)
  future::plan("multicore", workers = RhpcBLASctl::blas_get_num_procs())
  rbindlist(future.apply::future_lapply(
    seq_len(nrow(synth_comps)),
    function(i) {
      this_comp <- synth_comps[i]
      pred_data <- cbind(pred_data, this_comp[, .(R1, R2, R3)])
      pred_data_long_cuts <- pred_data[rep(
        seq_len(nrow(pred_data)),
        each = final_time
      )]
      pred_data_long_cuts[, timegroup := rep(1:final_time, nrow(pred_data))]

      pred_data_long_cuts[, `:=`(
        haz_dem = predict(models$model_dem, newdata = .SD, type = "response"),
        haz_death = predict(
          models$model_death,
          newdata = .SD,
          type = "response"
        )
      )]

      setkey(pred_data_long_cuts, id, timegroup)
      pred_data_long_cuts[,
        risk := cumsum(
          haz_dem *
            cumprod((1 - shift(haz_dem, 1, 0)) * (1 - haz_death))
        ),
        by = id
      ]
      final_risk <- pred_data_long_cuts[
        timegroup == final_time,
        weighted.mean(risk, prop)
      ]

      this_comp[, risk := final_risk]
    }
  ))
}
