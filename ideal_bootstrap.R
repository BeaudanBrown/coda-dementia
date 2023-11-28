source("ideal_comp.R")

ncpus <- 1
bootstrap_iterations <- ncpus

best_and_worst <- get_best_and_worst_comp()

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))
boot_copy <- boot_data

# set date variables to strings to avoid errors
boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)

boot_data <- boot_data[1:5000, ]

run_cum_bootstrap <- function(timegroup, output_name) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(boot_data,
    mincor = 0,
    exclude = c(
      "date_acdem2", "date_accel", "date_of_death",
      "avg_sleep", "avg_inactivity", "avg_light",
      "avg_mvpa"
    )
  )
  predmat["date_acdem2", ] <- 0
  predmat["date_of_death", ] <- 0
  predmat["date_accel", ] <- 0

  # method for each imputed variable
  imp_methods <- make.method(boot_data)
  # exclude dates from being imputed
  imp_methods["date_acdem2"] <- ""
  imp_methods["date_of_death"] <- ""
  imp_methods["date_accel"] <- ""

  result <- boot(
    data = boot_data,
    statistic = bootstrap_ideal_fn,
    create_formula_fn = get_primary_formula,
    timegroup = timegroup,
    best_and_worst = best_and_worst,
    predmat = predmat,
    imp_methods = imp_methods,
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
  timegroup,
  predmat,
  imp_methods,
  best_and_worst
) {
  this_sample <- data[indices, ]

  print("Imputing")
  print(format(Sys.time(), "%H:%M:%S"))
  imp <- mice(this_sample, m = 1, maxit = 1, predictorMatrix = predmat, methods = imp_methods)
  imp <- complete(imp)
  setDT(imp)
  imp[, id := .I]
  imp_len <- nrow(imp)
  # fit model
  print("Fitting model")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  model <- fit_model(imp, create_formula_fn)

  # data for g-computation/standardisation
  imp <- imp[rep(seq_len(imp_len), each = timegroup)]
  imp[, timegroup := rep(1:timegroup, imp_len)]
  setkey(imp, id, timegroup) # sort and set keys for efficient grouping and joining

  best_ilr = ilr(acomp(best_and_worst$best), V = v)
  worst_ilr = ilr(acomp(best_and_worst$worst), V = v)
  common_ilr = ilr(acomp(best_and_worst$most_common), V = v)

  ## Best composition risk
  imp[, c("R1", "R2", "R3") := list(best_ilr[1], best_ilr[2], best_ilr[3])]
  imp[, haz := predict(model, newdata = .SD, type = "response")]
  imp[, risk := 1 - cumprod(1 - haz), by = id]
  best_risk <- imp %>%
    group_by(timegroup) %>%
    summarise(best_risk = mean(risk, na.rm = TRUE))

  ## Worst composition risk
  imp[, c("R1", "R2", "R3") := list(worst_ilr[1], worst_ilr[2], worst_ilr[3])]
  imp[, haz := predict(model, newdata = .SD, type = "response")]
  imp[, risk := 1 - cumprod(1 - haz), by = id]
  worst_risk <- imp %>%
    group_by(timegroup) %>%
    summarise(worst_risk = mean(risk, na.rm = TRUE))

  ## Most common composition risk
  imp[, c("R1", "R2", "R3") := list(common_ilr[1], common_ilr[2], common_ilr[3])]
  imp[, haz := predict(model, newdata = .SD, type = "response")]
  setkey(imp, id, timegroup) # sort and set keys for efficient grouping and joining
  imp[, risk := 1 - cumprod(1 - haz), by = id]
  common_risk <- imp %>%
    group_by(timegroup) %>%
    summarise(common_risk = mean(risk, na.rm = TRUE))

  print("Returning results")
  print(format(Sys.time(), "%H:%M:%S"))
  full_df <- full_join(best_risk, worst_risk, by = "timegroup") |>
    full_join(common_risk, by = "timegroup")
  return(as.matrix(full_df))
}

process_ideal_output <- function(rds_path) {
  data <- readRDS(file.path(data_dir, rds_path))
  num_timegroups <- 76

  plot_data <- as.data.frame(data$t0) |>
    pivot_longer(
      cols = -timegroup,
      names_to = "Reference",
      values_to = "Risk",
      names_pattern = "(.*)_risk"
    ) |>
    mutate(
      Reference = case_when(
        Reference == "best"   ~ "best",
        Reference == "worst"  ~ "worst",
        Reference == "common" ~ "common"
      )
    )

  slice <- data$t[, start:(start + num_timegroups - 1)]

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

  all_quantiles <- rbind(best_quantiles,
                         worst_quantiles,
                         common_quantiles)


  plot_data <- full_join(plot_data, all_quantiles, by = c("timegroup", "Reference"))

  ggplot(plot_data, aes(x = timegroup, y = Risk)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.25) +
    facet_wrap(~ Reference) +
    cowplot::theme_cowplot()
}

process_ideal_output("ideal.rds")
