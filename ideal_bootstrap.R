source("ideal_comp.R")
source("utils.R")

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

best_and_worst <- get_best_and_worst_comp(boot_data)

run_cum_bootstrap <- function(data, output_name) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(data,
    mincor = 0,
    exclude = c(
      "avg_sleep", "avg_inactivity", "avg_light",
      "avg_mvpa", "eid", "time_to_dem"
    )
  )

  min_age_of_dem <- min(data$age_dem)
  max_age_of_dem <- max(data$age_dem)
  age_range <- max_age_of_dem - min_age_of_dem
  num_timegroups <- ceiling(age_range * 2)

  result <- boot(
    data = boot_data,
    statistic = bootstrap_ideal_fn,
    create_formula_fn = get_primary_formula,
    best_and_worst = best_and_worst,
    predmat = predmat,
    num_timegroups = num_timegroups,
    R = bootstrap_iterations,
    parallel = "multicore",
    ncpus = ncpus
  )

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  output_name_with_timestamp <- paste0(output_name, "_", timestamp, ".rds")
  saveRDS(result, file.path(output_dir, output_name_with_timestamp))
}

bootstrap_ideal_fn <- function(
    data,
    indices,
    create_formula_fn,
    predmat,
    best_and_worst,
    num_timegroups) {
  this_sample <- data[indices, ]

  print("Imputing")
  print(format(Sys.time(), "%H:%M:%S"))
  imp <- mice(this_sample, m = 1, predictorMatrix = predmat, maxit = maxit)
  imp <- complete(imp)
  setDT(imp)
  imp[, id := .I]
  imp_len <- nrow(imp)
  # fit model
  print("Fitting model")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  models <- fit_model(imp, create_formula_fn)
  dem_model <- models[["model_dem"]]
  death_model <- models[["model_death"]]

  imp <- imp[rep(seq_len(imp_len), each = num_timegroups)]
  imp[, timegroup := rep(1:num_timegroups, imp_len)]

  best_ilr <- ilr(acomp(best_and_worst$best), V = v)
  worst_ilr <- ilr(acomp(best_and_worst$worst), V = v)
  common_ilr <- ilr(acomp(best_and_worst$most_common), V = v)

  calculate_risk <- function(risk_data, ilr_values, risk_name) {
    newx <- model.matrix(create_formula_fn(risk_data), risk_data)

    newx[, "R1"] <- ilr_values[[1]]
    newx[, "I(R1^2)"] <- ilr_values[[1]]^2
    newx[, "R2"] <- ilr_values[[2]]
    newx[, "I(R2^2)"] <- ilr_values[[2]]^2
    newx[, "R3"] <- ilr_values[[3]]
    newx[, "I(R3^2)"] <- ilr_values[[3]]^2

    risk_data[, haz_dem := predict(
      dem_model,
      newdata = newx, type = "response"
    )]
    risk_data[, haz_death := predict(
      death_model,
      newdata = newx, type = "response"
    )]
    risk_data[, risk := cumsum(
      haz_dem * cumprod(1 - lag(haz_dem, default = 0) * (1 - haz_death))
    ), by = id]
    summarised_risk <- risk_data[
      , .(avg_risk = mean(risk, na.rm = TRUE)),
      by = .(timegroup)
    ]
    setnames(summarised_risk, c("avg_risk"), c(risk_name))
    return(summarised_risk)
  }

  best_risk <- calculate_risk(
    imp, list(best_ilr[1], best_ilr[2], best_ilr[3]), "best"
  )
  worst_risk <- calculate_risk(
    imp, list(worst_ilr[1], worst_ilr[2], worst_ilr[3]), "worst"
  )
  common_risk <- calculate_risk(
    imp, list(common_ilr[1], common_ilr[2], common_ilr[3]), "common"
  )

  print("Returning results")
  print(format(Sys.time(), "%H:%M:%S"))
  full_df <- full_join(best_risk, worst_risk, by = "timegroup") |>
    full_join(common_risk, by = "timegroup")
  return(as.matrix(full_df))
}

process_ideal_output <- function(rds_path) {
  data <- readRDS(file.path(output_dir, rds_path))
  num_timegroups <- 78

  plot_data <- as.data.frame(data$t0) |>
    pivot_longer(
      cols = -timegroup,
      names_to = "Reference",
      values_to = "Risk"
    )

  get_quantiles <- function(start, reference) {
    slice <- data$t[, start:(start + num_timegroups - 1)]

    quantiles <-
      as.data.frame(t(as.data.frame(apply(slice, 2, function(column) quantile(column, probs = c(0.025, 0.975))))))
    colnames(quantiles) <- c("lower", "upper")
    quantiles$Reference <- reference
    quantiles$timegroup <- 1:num_timegroups
    return(quantiles)
  }

  reference_start <- num_timegroups + 1

  best_quantiles <- get_quantiles(reference_start, "best")
  reference_start <- reference_start + num_timegroups

  worst_quantiles <- get_quantiles(reference_start, "worst")
  reference_start <- reference_start + num_timegroups

  common_quantiles <- get_quantiles(reference_start, "common")
  reference_start <- reference_start + num_timegroups

  all_quantiles <- rbind(
    best_quantiles,
    worst_quantiles,
    common_quantiles
  )


  plot_data <- full_join(plot_data, all_quantiles, by = c("timegroup", "Reference"))

  plot_data |>
    mutate(age = 47 + (timegroup / 2)) |>
    mutate(Composition = fct_recode(Reference,
      "Ideal" = "best",
      "Typical" = "common",
      "Worst" = "worst"
    )) |>
    ggplot(aes(x = age, y = Risk)) +
    geom_line(aes(colour = Composition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Composition),
      alpha = 0.25
    ) +
    labs(x = "Age", y = "Cumulative all-cause dementia incidence") +
    cowplot::theme_cowplot() +
    scale_color_manual(
      labels = c("Ideal", "Typical", "Worst"),
      values = c("#7AC36A", "#56B4E9", "#DC3912")
    ) +
    scale_fill_manual(
      labels = c("Ideal", "Typical", "Worst"),
      values = c("#7AC36A", "#56B4E9", "#DC3912")
    ) +
    theme(text = element_text(family = "serif"))
}
