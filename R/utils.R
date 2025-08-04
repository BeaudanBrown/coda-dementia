Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

bootstrap_sample <- function(data) {
  out <- data[sample(nrow(data), nrow(data), TRUE)]
  out[, id := .I]
  return(out)
}

## Define SBP
sbp <- matrix(
  c(
    1,
    1,
    -1,
    -1,
    1,
    -1,
    0,
    0,
    0,
    0,
    1,
    -1
  ),
  ncol = 4,
  byrow = TRUE
)

v <- compositions::gsi.buildilrBase(t(sbp))


## strip unneccessary model components

strip_lm <- function(cm) {
  cm$y <- c()
  cm$model <- c()

  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()

  cm
}

strip_glm <- function(cm) {
  cm$y <- c()
  cm$model <- c()

  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()

  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()

  return(cm)
}

save_plot <- function(plot, file_path) {
  ggsave(
    file_path,
    plot = plot,
    device = "svg",
    width = 10,
    height = 12
  )
}

ordinal_to_numeric <- function(data) {
  ## set binary/ordinal factor variables to numeric to allow use of pmm
  data$ethnicity <- as.numeric(data$ethnicity)
  data$avg_total_household_income <- ifelse(
    data$avg_total_household_income == "<18",
    0,
    ifelse(
      data$avg_total_household_income == "18-30",
      1,
      ifelse(data$avg_total_household_income == "31-50", 2, 3)
    )
  )
  data$highest_qual <- ifelse(
    data$highest_qual == "O",
    0,
    ifelse(
      data$highest_qual == "NVQ",
      1,
      ifelse(data$highest_qual == "A", 2, 3)
    )
  )
  data$smok_status <- as.numeric(data$smok_status) - 1
  data$bp_med <- as.numeric(data$bp_med) - 1
  data$any_cvd <- as.numeric(data$any_cvd) - 1
  data$chronotype <- as.numeric(data$chronotype) - 1
  data$apoe_e4 <- as.numeric(data$apoe_e4) - 1

  data$freq_depressed_twoweeks <- ifelse(
    data$freq_depressed_twoweeks == "not at all",
    0,
    ifelse(
      data$freq_depressed_twoweeks == "several days",
      1,
      ifelse(data$freq_depressed_twoweeks == "more than half", 2, 3)
    )
  )
  data$diagnosed_diabetes <- as.numeric(data$diagnosed_diabetes) - 1

  data
}

apply_substitution <- function(imp, from_var, to_var, duration, comp_limits) {
  lower_from <- comp_limits[[from_var]]$lower
  upper_to <- comp_limits[[to_var]]$upper

  max_from_change <- imp[[from_var]] - lower_from
  max_to_change <- upper_to - imp[[to_var]]
  can_substitute <- (max_from_change >= duration) & (max_to_change >= duration)
  prop_substituted <- sum(can_substitute) / nrow(imp)

  sub <- imp
  sub[[from_var]] <- sub[[from_var]] - (can_substitute * duration)
  sub[[to_var]] <- sub[[to_var]] + (can_substitute * duration)

  comp <- acomp(sub[, c(
    "avg_sleep",
    "avg_inactivity",
    "avg_light",
    "avg_mvpa"
  )])
  ilr_vars <- ilr(comp, V = v) |>
    setNames(c("R1", "R2", "R3"))

  sub[, c("R1", "R2", "R3")] <- as.data.table(ilr_vars)
  data.table(
    results = list(sub),
    B = unique(imp$tar_batch),
    from_var = from_var,
    to_var = to_var,
    duration = duration,
    prop_substituted = prop_substituted
  )
}

make_plot_grid <- function(
  plot_data,
  grid_layout
) {
  plot_data[, cohort := factor(cohort, levels = grid_layout$cohort_order)]
  plot_data[, sub_type := factor(sub_type, levels = grid_layout$subtype_order)]

  title_row <- lapply(grid_layout$cohort_order, function(cohort) {
    cohort_name <- switch(
      cohort,
      short_sleeper = "Short Sleepers",
      avg_sleeper = "Normal Sleepers",
      long_sleeper = "Long Sleepers",
      full_cohort = "Full Cohort",
      s1 = "WASO Adjusted",
      s2 = "Comorbidity Adjusted",
      s3 = "Interactions Adjusted",
      representative = "Representative Covar Adjusted",
      reverse_causation = "Reverse Causation Adjusted",
      retired = "Retired",
      not_retired = "Not Retired",
      cohort
    )
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = cohort_name,
        size = 5,
        vjust = 0,
        hjust = 0.5
      ) +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      theme_void()
  })
  grid_list <- lapply(grid_layout$subtype_order, function(this_subtype) {
    plots_row <- plot_data[sub_type == this_subtype][order(cohort), plot]
    wrap_plots(c(plots_row), nrow = 1)
  })
  patchwork::wrap_plots(
    c(list(wrap_plots(c(title_row), nrow = 1)), grid_list),
    heights = c(0.2, rep(1, length(grid_list))),
    ncol = 1
  )
}

save_plots <- function(plot_pattern = "plot_grid", sc = TRUE) {
  source("primary_targets.R")
  source("mri_targets.R")
  source("cum_targets.R")
  source("sensitivity_reverse_causation_targets.R")
  source("sensitivity_covar_targets.R")
  source("sensitivity_representative_targets.R")
  quint_targets <- tarchetypes::tar_select_names(
    list(
      covar_sensitivity_targets
    ),
    tidyr::contains(plot_pattern)
  )
  triple_targets <- tarchetypes::tar_select_names(
    list(
      primary_targets,
      mri_targets
    ),
    tidyr::contains(plot_pattern)
  )
  # double_targets <- tarchetypes::tar_select_names(
  #   list(
  #     primary_targets,
  #     mri_targets
  #   ),
  #   tidyr::contains(plot_pattern)
  # )
  single_targets <- tarchetypes::tar_select_names(
    list(
      cum_targets
    ),
    tidyr::contains(plot_pattern)
  )
  save_plots_cust(
    plot_width = 25,
    plot_targets = quint_targets,
    sc = sc,
    plot_pattern = plot_pattern
  )
  save_plots_cust(
    plot_width = 15,
    plot_targets = triple_targets,
    sc = sc,
    plot_pattern = plot_pattern
  )
  # save_plots_cust(
  #   plot_width = 20,
  #   plot_targets = quad_targets,
  #   sc = sc,
  #   plot_pattern = plot_pattern
  # )
  # save_plots_cust(
  #   plot_width = 10,
  #   plot_targets = double_targets,
  #   sc = sc,
  #   plot_pattern = plot_pattern
  # )
  save_plots_cust(
    plot_width = 5,
    plot_targets = single_targets,
    sc = sc,
    plot_pattern = plot_pattern
  )
}

save_plots_cust <- function(
  plot_dir = "plots",
  plot_width = 15,
  plot_height = 12,
  plot_targets,
  plot_pattern = "plot_grid",
  sc = TRUE
) {
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  tar_make(names = plot_targets, shortcut = sc)
  for (target in plot_targets) {
    name_part <- sub(paste0(plot_pattern, "_"), "", target)
    if (name_part == "") {
      name_part <- target
    }
    file_path <- file.path(plot_dir, paste0(name_part, ".png"))
    plot_obj <- tar_read_raw(target)
    ggplot2::ggsave(
      file_path,
      plot_obj,
      width = plot_width,
      height = plot_height
    )
  }
}
