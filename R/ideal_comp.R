generate_compositions <- function() {
  step_size <- 15

  comp_limits <- list(
    avg_sleep = list(
      lower = ceiling(181 / step_size) * step_size,
      upper = floor(544 / step_size) * step_size
    ),
    avg_inactivity = list(
      lower = ceiling(348 / step_size) * step_size,
      upper = floor(1059 / step_size) * step_size
    ),
    avg_light = list(
      lower = ceiling(76 / step_size) * step_size,
      upper = floor(511 / step_size) * step_size
    ),
    avg_mvpa = list(
      lower = ceiling(20 / step_size) * step_size,
      upper = floor(384 / step_size) * step_size
    )
  )

  # Create an empty data frame
  df <- expand_grid(
    avg_sleep = seq(
      comp_limits[["avg_sleep"]]$lower,
      comp_limits[["avg_sleep"]]$upper,
      by = step_size
    ),
    avg_inactivity = seq(
      comp_limits[["avg_inactivity"]]$lower,
      comp_limits[["avg_inactivity"]]$upper,
      by = step_size
    ),
    avg_light = seq(
      comp_limits[["avg_light"]]$lower,
      comp_limits[["avg_light"]]$upper,
      by = step_size
    ),
    avg_mvpa = seq(
      comp_limits[["avg_mvpa"]]$lower,
      comp_limits[["avg_mvpa"]]$upper,
      by = step_size
    )
  )

  setDT(df)

  # Calculate total times
  df[, total := rowSums(.SD)]

  # Filter out rows where total time is not 24 * 60
  df <- df[total == mins_in_day, ]
  df[, total := NULL]

  df
}

add_density <- function(df, synth_comps, dens_threshold) {
  # only need ILR variables for df
  df <- df[, list(R1, R2, R3)]
  # add ILRs for the synth comps
  synth_comps <- cbind(
    synth_comps,
    as.data.table(ilr(acomp(synth_comps), V = v)) |>
      setnames(c("R1", "R2", "R3"))
  )
  # estimate params of multivariate dist using kernel dens estimator
  df <- as.matrix(scale(df))
  fit <- kde(df, gridsize = rep(100, 3))
  # density for observed compositions
  dens_obs <- log(predict(fit, x = df))
  dens_threshold <- quantile(dens_obs, dens_threshold)
  # predict density for each composition
  synth_comps$dens <- log(
    predict(
      fit,
      x = as.matrix(scale(synth_comps[, list(R1, R2, R3)]))
    )
  )

  synth_comps$dens_threshold <- dens_threshold
  synth_comps
}

generate_hazards <- function(comps, model, model_formula, ref_row) {
  ilrs <- t(apply(
    comps[, c(1, 2, 3, 4)],
    1,
    function(comp) ilr(acomp(comp), V = v)
  ))

  # Convert contrasts into a data frame
  risks <- apply(ilrs, 1, function(new_ilr) {
    ref_row[c("R1", "R2", "R3")] <- new_ilr
    newx <- model.matrix(model_formula, ref_row)
    predict(model, newdata = newx, type = "response")
  })
  return(risks)
}

get_best_and_worst_comp <- function(df, timegroup_cuts = NULL) {
  if (is.null(timegroup_cuts)) {
    min_age_of_dem <- min(df$age_dem)
    max_age_of_dem <- max(df$age_dem)
    age_range <- max_age_of_dem - min_age_of_dem
    timegroup_steps <- ceiling(age_range * 2)
    median_age_of_dem <- median(df[df$dem == 1, ]$age_dem)

    timegroup_cuts <-
      seq(
        from = min_age_of_dem,
        to = max_age_of_dem,
        length.out = timegroup_steps
      )
  }
  ## Load data
  predmat <- quickpred(
    df,
    mincor = 0,
    exclude = c(
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa",
      "eid",
      "time_to_dem"
    )
  )
  dem_df <- mice(df, m = 1, predictorMatrix = predmat, maxit = maxit)
  dem_df <- complete(dem_df)

  dem_comp_df <- as.data.frame(
    dem_df[, c("avg_sleep", "avg_inactivity", "avg_mvpa", "avg_light")]
  )

  include_cols <- c("avg_sleep", "avg_inactivity", "avg_mvpa")
  mean_vec <- as.numeric(apply(dem_comp_df[, include_cols], 2, mean))
  vcv_mat <- cov(dem_comp_df[, include_cols])
  dem_comp_df$dens <- apply(
    dem_comp_df[, include_cols],
    1,
    dmvnorm,
    mean = mean_vec,
    sigma = vcv_mat,
    log = TRUE
  )

  threshold <- quantile(dem_comp_df$dens, probs = c(0.025), na.rm = TRUE)

  quantiles <- apply(dem_comp_df, 2, function(column) {
    quantile(column, probs = c(0.025, 0.975), na.rm = TRUE)
  })

  # Store in "lower" and "upper" vectors
  lower <- quantiles[1, ]
  upper <- quantiles[2, ]

  generated_comps <- generate_compositions(lower, upper)
  generated_comps$dens <-
    apply(
      generated_comps[, c("sleep", "inactive", "modvig")],
      1,
      dmvnorm,
      mean = mean_vec,
      sigma = vcv_mat,
      log = TRUE
    )
  generated_comps <- generated_comps[generated_comps$dens > threshold, ]

  # fit model
  models <- fit_model(
    dem_df,
    timegroup_cuts,
    get_primary_formula
  )
  dem_model <- models[["model_dem"]]
  model_formula <- models[["model_formula"]]
  ref_row <- as.data.frame(dem_df[1, ])
  # Arbitrary, since hazards will be equivalently ranked across timegroups
  ref_row$timegroup <- 55

  generated_comps$haz <- generate_hazards(
    generated_comps,
    dem_model,
    model_formula,
    ref_row
  )
  generated_comps <- generated_comps[order(generated_comps$haz), ]
  best_comp <- acomp(generated_comps[1, c(1, 2, 3, 4)])
  worst_comp <- acomp(generated_comps[nrow(generated_comps), c(1, 2, 3, 4)])
  typical_comp <- acomp(generated_comps[
    generated_comps$dens == max(generated_comps$dens),
    c(1, 2, 3, 4)
  ])

  result <- as.data.frame(t(rbind(best_comp, worst_comp, typical_comp)))
  colnames(result) <- c("best", "worst", "typical")

  result
}
