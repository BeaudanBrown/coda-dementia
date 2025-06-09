process_mri_output <- function(data) {
  # tidy bootstrap output
  ## prepare data for plotting
  # Define the phenos and comps
  phenos <- rep(c("tbv", "wmv", "gmv", "hip", "log_wmh"), each = 3)
  comps <- rep(c("worst", "typical", "ideal"), times = 5)

  col_names <- paste(comps, phenos, sep = "_")
  boot_reps <- as_tibble(data$t)
  colnames(boot_reps) <- col_names

  # Combine the phenos, comps, and values into a new data frame
  plot_data <- data.frame(comp = comps, pheno = phenos, value = data$t0)

  num_comps <- 3

  get_quantiles <- function(start, pheno) {
    slice <- data$t[, start:(start + num_comps - 1)]

    quantiles <-
      as.data.frame(t(apply(
        slice,
        2,
        function(column) quantile(column, probs = c(0.025, 0.975))
      )))
    colnames(quantiles) <- c("lower", "upper")
    comps <- c("worst", "typical", "ideal")
    quantiles$comp <- comps
    quantiles$pheno <- pheno
    return(quantiles)
  }

  pheno_start <- 1

  tbv_quants <- get_quantiles(pheno_start, "tbv")
  pheno_start <- pheno_start + num_comps

  wmv_quants <- get_quantiles(pheno_start, "wmv")
  pheno_start <- pheno_start + num_comps

  gmv_quants <- get_quantiles(pheno_start, "gmv")
  pheno_start <- pheno_start + num_comps

  hip_quants <- get_quantiles(pheno_start, "hip")
  pheno_start <- pheno_start + num_comps

  log_wmh_quants <- get_quantiles(pheno_start, "log_wmh")

  all_quantiles <- rbind(
    tbv_quants,
    wmv_quants,
    gmv_quants,
    hip_quants,
    log_wmh_quants
  )

  plot_data <- full_join(plot_data, all_quantiles, by = c("comp", "pheno"))

  return(list(plot_data, boot_reps))
}


### Plot

plot_mri <- function(result_list, mri_model_data) {
  plot_data <- result_list[[1]]

  plot_data$comp <- str_to_title(plot_data$comp)
  plot_data$comp <- ifelse(
    plot_data$comp == "Common",
    "Typical",
    plot_data$comp
  )
  plot_data$comp <- factor(
    plot_data$comp,
    levels = c("Worst", "Typical", "Ideal")
  )

  single_plot <- function(pheno) {
    plot_data[plot_data$pheno == pheno, ] |>
      ggplot(aes(x = comp, y = value)) +
      geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.1) +
      labs(x = "Composition") +
      theme_cowplot(
        font_size = 12,
        font_family = "serif",
        line_size = 0.25,
      ) +
      theme(
        panel.border = element_rect(fill = NA, colour = "#585656"),
        panel.grid = element_line(colour = "grey92"),
        panel.grid.minor = element_line(linewidth = rel(0.5)),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color = "#585656")
      )
  }

  tbv_plot <- single_plot("tbv")
  gmv_plot <- single_plot("gmv")
  wmv_plot <- single_plot("wmv")
  hip_plot <- single_plot("hip")
  wmh_plot <- single_plot("log_wmh")

  plot_grid(
    tbv_plot +
      labs(
        y = expression(paste(
          "Total brain volume ",
          (cm^{
            3
          })
        ))
      ) +
      labs(x = "") +
      ylim(c(
        mean(mri_model_data$tbv, na.rm = TRUE) -
          0.5 *
            sd(mri_model_data$tbv, na.rm = TRUE),
        mean(mri_model_data$tbv, na.rm = TRUE) +
          0.5 *
            sd(mri_model_data$tbv, na.rm = TRUE)
      )),
    gmv_plot +
      labs(
        y = expression(paste(
          "Grey matter volume ",
          (cm^{
            3
          })
        ))
      ) +
      labs(x = "") +
      ylim(c(
        mean(mri_model_data$gmv, na.rm = TRUE) -
          0.5 *
            sd(mri_model_data$gmv, na.rm = TRUE),
        mean(mri_model_data$gmv, na.rm = TRUE) +
          0.5 *
            sd(mri_model_data$gmv, na.rm = TRUE)
      )),
    wmv_plot +
      labs(
        y = expression(paste(
          "White matter volume ",
          (cm^{
            3
          })
        )),
        x = ""
      ) +
      ylim(c(
        mean(mri_model_data$wmv, na.rm = TRUE) -
          0.5 *
            sd(mri_model_data$wmv, na.rm = TRUE),
        mean(mri_model_data$wmv, na.rm = TRUE) +
          0.5 *
            sd(mri_model_data$wmv, na.rm = TRUE)
      )),
    hip_plot +
      labs(
        y = expression(paste(
          "Hippocampal volume ",
          (cm^{
            3
          })
        ))
      ) +
      ylim(c(
        mean(mri_model_data$hip, na.rm = TRUE) -
          0.5 *
            sd(mri_model_data$hip, na.rm = TRUE),
        mean(mri_model_data$hip, na.rm = TRUE) +
          0.5 *
            sd(mri_model_data$hip, na.rm = TRUE)
      )),
    wmh_plot +
      labs(y = "Log white matter hyperintensities") +
      ylim(c(
        mean(mri_model_data$log_wmh, na.rm = TRUE) -
          0.5 *
            sd(mri_model_data$log_wmh, na.rm = TRUE),
        mean(mri_model_data$log_wmh, na.rm = TRUE) +
          0.5 *
            sd(mri_model_data$log_wmh, na.rm = TRUE)
      )),
    nrow = 3
  )
}

run_mri_bootstrap <- function(mri_df, mri_best_and_worst, create_formula_fn) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(
    mri_df,
    mincor = 0,
    exclude = c(
      "eid",
      "calendar_date",
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa"
    )
  )

  make_ilr <- function(comp) {
    return(ilr(acomp(comp), V = v))
  }

  mri_best_and_worst <- as.data.frame(apply(mri_best_and_worst, 2, make_ilr))

  result <- boot(
    data = mri_df,
    statistic = bootstrap_mri_fn,
    create_formula_fn = create_formula_fn,
    predmat = predmat,
    best_and_worst = mri_best_and_worst,
    R = bootstrap_iterations,
    parallel = "multicore",
    ncpus = ncpus
  )

  return(result)
}

bootstrap_mri_fn <- function(
  data,
  indices,
  create_formula_fn,
  predmat,
  best_and_worst
) {
  this_sample <- data[indices, ]

  imp <- mice(
    this_sample,
    m = 1,
    maxit = maxit,
    predictorMatrix = predmat
  )
  imp <- complete(imp)

  result_df <- predict_mri_results(imp, best_and_worst)

  return(as.matrix(result_df))
}
