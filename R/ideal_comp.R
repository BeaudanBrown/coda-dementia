generate_compositions <- function(df) {
  step_size <- 15

  sleep_quants <- quantile(
    df$avg_sleep,
    probs = c(uni_threshold, 1 - uni_threshold)
  )
  inactivity_quants <- quantile(
    df$avg_inactivity,
    probs = c(uni_threshold, 1 - uni_threshold)
  )
  light_quants <- quantile(
    df$avg_light,
    probs = c(uni_threshold, 1 - uni_threshold)
  )
  mvpa_quants <- quantile(
    df$avg_mvpa,
    probs = c(uni_threshold, 1 - uni_threshold)
  )

  comp_limits <- list(
    avg_sleep = list(
      lower = ceiling(sleep_quants[1] / step_size) * step_size,
      upper = floor(sleep_quants[2] / step_size) * step_size
    ),
    avg_inactivity = list(
      lower = ceiling(inactivity_quants[1] / step_size) * step_size,
      upper = floor(inactivity_quants[2] / step_size) * step_size
    ),
    avg_light = list(
      lower = ceiling(light_quants[1] / step_size) * step_size,
      upper = floor(light_quants[2] / step_size) * step_size
    ),
    avg_mvpa = list(
      lower = ceiling(mvpa_quants[1] / step_size) * step_size,
      upper = floor(mvpa_quants[2] / step_size) * step_size
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
