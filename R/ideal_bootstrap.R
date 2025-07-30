produce_ideal_plot <- function() {
  ideal_rds <- file.path(data_dir, "boot_ideal.rds")
  ideal_plot <- process_ideal_output(ideal_rds)
  ggsave(
    file.path(
      data_dir,
      "../Manuscript/Main_figures/Cumulative.svg"
    ),
    ideal_plot,
    device = "svg",
    width = 10,
    height = 10.5
  )
}

run_cum_bootstrap <- function(output_name, intervals = TRUE) {
  data <- read_rds(boot_data_file)
  # Matrix of variables to include in imputation model
  predmat <- quickpred(
    data,
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

  min_age_of_dem <- min(data$age_dem)
  max_age_of_dem <- max(data$age_dem)
  age_range <- max_age_of_dem - min_age_of_dem
  num_timegroups <- ceiling(age_range * 2)

  timegroup_cuts <-
    seq(
      from = min_age_of_dem,
      to = max_age_of_dem,
      length.out = num_timegroups
    )

  ## Split data into two folds - one for finding best and worst comp
  ## and another for estimating dementia risk for identified comps

  set.seed(seed_val)
  split <- sample(1:nrow(data), size = floor(nrow(data) / 2), replace = FALSE)

  fold1 <- data[split, ]
  fold2 <- data[-split, ]

  filename <- paste("best_and_worst", seed_val, sep = "-")
  best_and_worst <-
    get_best_and_worst_comp(fold1, timegroup_cuts, filename = filename)

  if (isTRUE(intervals)) {
    library(parallel)
    outfile <- paste("ideal_outfile", seed_val, sep = "-")
    cl <- makeCluster(ncpus, outfile = paste0(outfile, ".txt"))
    # clusterEvalQ(cl, {
    #   source("utils.R")
    #   source("ideal_comp.R")
    # })

    result <- boot(
      data = fold2,
      statistic = bootstrap_ideal_fn,
      create_formula_fn = get_primary_formula,
      predmat = predmat,
      best_and_worst = best_and_worst,
      timegroup_cuts = timegroup_cuts,
      num_timegroups = num_timegroups,
      R = bootstrap_iterations,
      parallel = "snow",
      cl = cl,
      ncpus = ncpus
    )
    stopCluster(cl)
  } else {
    result <- bootstrap_ideal_fn(
      data = fold2,
      indices = seq(1, nrow(data)),
      create_formula_fn = get_primary_formula,
      predmat = predmat,
      best_and_worst = best_and_worst,
      timegroup_cuts = timegroup_cuts,
      num_timegroups = num_timegroups
    )
  }

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
  timegroup_cuts,
  num_timegroups
) {
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
  models <- fit_model(imp, timegroup_cuts, create_formula_fn)
  dem_model <- models[["model_dem"]]
  death_model <- models[["model_death"]]

  imp <- imp[rep(seq_len(imp_len), each = num_timegroups)]
  imp[, timegroup := rep(1:num_timegroups, imp_len)]

  best_ilr <- ilr(acomp(best_and_worst$best), V = v)
  worst_ilr <- ilr(acomp(best_and_worst$worst), V = v)
  common_ilr <- ilr(acomp(best_and_worst$typical), V = v)

  calculate_risk <- function(model_formula, risk_data, ilr_values, risk_name) {
    newx <- model.matrix(model_formula, risk_data)

    newx[, "R1"] <- ilr_values[[1]]
    newx[, "I(R1^2)"] <- ilr_values[[1]]^2
    newx[, "R2"] <- ilr_values[[2]]
    newx[, "I(R2^2)"] <- ilr_values[[2]]^2
    newx[, "R3"] <- ilr_values[[3]]
    newx[, "I(R3^2)"] <- ilr_values[[3]]^2

    risk_data[,
      haz_dem := predict(
        dem_model,
        newdata = newx,
        type = "response"
      )
    ]
    risk_data[,
      haz_death := predict(
        death_model,
        newdata = newx,
        type = "response"
      )
    ]
    risk_data[,
      risk := cumsum(
        haz_dem * cumprod(1 - lag(haz_dem, default = 0) * (1 - haz_death))
      ),
      by = id
    ]
    summarised_risk <- risk_data[,
      .(avg_risk = mean(risk, na.rm = TRUE)),
      by = .(timegroup)
    ]
    setnames(summarised_risk, c("avg_risk"), c(risk_name))
    return(summarised_risk)
  }

  print(paste(format(Sys.time(), "%H:%M:%S"), "Calculating best risk"))
  print(gc())
  best_risk <- calculate_risk(
    models[["model_formula"]],
    imp,
    list(best_ilr[1], best_ilr[2], best_ilr[3]),
    "best"
  )
  print(paste(format(Sys.time(), "%H:%M:%S"), "Calculating worst risk"))
  print(gc())
  worst_risk <- calculate_risk(
    models[["model_formula"]],
    imp,
    list(worst_ilr[1], worst_ilr[2], worst_ilr[3]),
    "worst"
  )
  print(paste(format(Sys.time(), "%H:%M:%S"), "Calculating common risk"))
  print(gc())
  common_risk <- calculate_risk(
    models[["model_formula"]],
    imp,
    list(common_ilr[1], common_ilr[2], common_ilr[3]),
    "common"
  )

  print(paste(format(Sys.time(), "%H:%M:%S"), "Returning results"))
  full_df <- full_join(best_risk, worst_risk, by = "timegroup") |>
    full_join(common_risk, by = "timegroup")
  return(as.matrix(full_df))
}

process_ideal_output <- function(rds_path, intervals = TRUE) {
  data <- readRDS(rds_path)
  num_timegroups <- 78
  if (isTRUE(intervals)) {
    full_sample_results <- data$t0
  } else {
    full_sample_results <- data
  }

  plot_data <- as.data.frame(full_sample_results) |>
    pivot_longer(
      cols = -timegroup,
      names_to = "Reference",
      values_to = "Risk"
    )

  if (isTRUE(intervals)) {
    get_quantiles <- function(start, reference) {
      slice <- data$t[, start:(start + num_timegroups - 1)]

      quantiles <-
        as.data.frame(
          t(
            as.data.frame(
              apply(slice, 2, function(column) {
                quantile(
                  column,
                  probs = c(0.025, 0.975)
                )
              })
            )
          )
        )
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

    plot_data <-
      full_join(
        plot_data,
        all_quantiles,
        by = c("timegroup", "Reference")
      )
  }

  p <- plot_data |>
    mutate(age = 47 + (timegroup / 2)) |>
    mutate(
      Composition = fct_recode(
        Reference,
        "Ideal" = "best",
        "Typical" = "common",
        "Worst" = "worst"
      )
    ) |>
    mutate(Composition = fct_rev(Composition)) |>
    ggplot(aes(x = age, y = Risk)) +
    geom_line(aes(colour = Composition)) +
    labs(x = "Age (years)", y = "Cumulative all-cause dementia incidence") +
    scale_color_manual(
      labels = c("Worst", "Typical", "Ideal"),
      values = c("#DC3912", "#56B4E9", "#7AC36A")
    ) +
    scale_fill_manual(
      labels = c("Worst", "Typical", "Ideal"),
      values = c("#DC3912", "#56B4E9", "#7AC36A")
    ) +
    cowplot::theme_cowplot(
      font_size = 16,
      font_family = "serif",
      line_size = 0.25
    ) +
    theme(
      panel.border = element_rect(fill = NA, colour = "#585656"),
      panel.grid = element_line(colour = "grey92"),
      panel.grid.minor = element_line(linewidth = rel(0.5)),
      axis.ticks.y = element_blank(),
      axis.line = element_line(color = "#585656"),
      axis.title.x = element_text(family = "serif", size = 20),
      axis.title.y = element_text(family = "serif", size = 20)
    )
  if (isTRUE(intervals)) {
    p <- p +
      geom_ribbon(
        aes(ymin = lower, ymax = upper, fill = Composition),
        alpha = 0.25
      )
  }
  return(p)
}

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

get_synth_risk <- function(pred_data, synth_comp, models, final_time) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)

  # append synthetic composition and predict risk
  pred_data <- cbind(pred_data, synth_comp[, list(R1, R2, R3)])
  pred_data_len <- nrow(pred_data)
  pred_data_long_cuts <- pred_data[rep(
    seq_len(pred_data_len),
    each = final_time
  )]
  pred_data_long_cuts[, timegroup := rep(1:final_time, pred_data_len)]

  pred_data_long_cuts[,
    haz_dem := predict(
      models$model_dem,
      newdata = .SD,
      type = "response"
    )
  ]

  pred_data_long_cuts[,
    haz_death := predict(
      models$model_death,
      newdata = .SD,
      type = "response"
    )
  ]
  setkey(pred_data_long_cuts, id, timegroup) # sort and set keys
  pred_data_long_cuts[,
    risk := cumsum(
      haz_dem * cumprod((1 - lag(haz_dem, default = 0)) * (1 - haz_death))
    ),
    by = id
  ]
  risk <- weighted.mean(
    pred_data_long_cuts[timegroup == final_time, ]$risk,
    w = pred_data_long_cuts[timegroup == final_time, ]$prop
  )

  synth_comp$risk <- risk
  synth_comp
}
