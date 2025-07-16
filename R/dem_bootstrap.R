# Generate bootstrap target definitions
generate_bootstrap_targets <- function(
  data_target_name,
  create_formula_fn_name,
  output_target_name,
  empirical = TRUE,
  intervals = TRUE,
  iterations = bootstrap_iterations,
  seed_val = 5678
) {
  if (FALSE) {
    data_target_name <- "df_raw"
    create_formula_fn_name <- "get_primary_formula"
    output_target_name <- "test_output"
    empirical <- FALSE
    intervals <- TRUE
    iterations <- 10
    seed_val <- 5678
  }
  # Create a list to hold all the bootstrap iteration targets
  bootstrap_targets <- list()

  # Create a setup target that computes everything needed for bootstrapping
  setup_target_name <- paste0(output_target_name, "_setup")
  bootstrap_targets[[1]] <- tar_target_raw(
    name = setup_target_name,
    command = substitute(
      {
        # Matrix of variables to include in imputation model
        data <- data_target_name
        predmat <- quickpred(
          data,
          mincor = 0,
          exclude = c(
            "avg_sleep",
            "avg_inactivity",
            "avg_light",
            "avg_mvpa",
            "eid"
          )
        )

        ## reference compositions
        all_comp <-
          acomp(data[,
            c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
          ])

        short_sleep_comp <-
          all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
        short_sleep_geo_mean <-
          acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

        avg_sleep_comp <-
          all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
        avg_sleep_geo_mean <-
          acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

        ## Constants for dementia risk
        min_age_of_dem <- min(data$age_dem)
        max_age_of_dem <- max(data$age_dem)
        age_range <- max_age_of_dem - min_age_of_dem
        timegroup_steps <- ceiling(age_range * 2)
        median_age_of_dem <- median(data[data$dem == 1, ]$age_dem)

        timegroup_cuts <-
          seq(
            from = min_age_of_dem,
            to = max_age_of_dem,
            length.out = timegroup_steps
          )

        median_age_of_dem_timegroup <-
          which(timegroup_cuts > median_age_of_dem)[1] - 1

        # Return everything needed for bootstrapping
        list(
          predmat = predmat,
          timegroup_cuts = timegroup_cuts,
          median_age_of_dem_timegroup = median_age_of_dem_timegroup,
          short_sleep_geo_mean = short_sleep_geo_mean,
          avg_sleep_geo_mean = avg_sleep_geo_mean,
          nrows = nrow(data),
          seed_val = seed_val
        )
      },
      list(
        data_target_name = as.name(data_target_name)
      )
    )
  )

  # Create a target for each bootstrap iteration
  for (i in 1:iterations) {
    iteration_name <- paste0(output_target_name, "_iter_", i)

    # Create an expression for the bootstrap iteration
    bootstrap_targets[[i + 1]] <- tar_target_raw(
      name = iteration_name,
      command = substitute(
        {
          setup <- setup_target
          data <- data_target
          create_formula_fn <- get(create_formula_fn_name)
          set.seed(seed_val + iter_num) # Unique seed for each iteration

          data <- ordinal_to_numeric(data)
          predmat <- quickpred(
            data,
            mincor = 0,
            exclude = c(
              "avg_sleep",
              "avg_inactivity",
              "avg_light",
              "avg_mvpa",
              "eid"
            )
          )

          # Generate random indices for bootstrap sample (with replacement)
          indices <- sample(seq_len(setup$nrows), setup$nrows, replace = TRUE)

          # Run the bootstrap function
          bootstrap_substitutions_fn(
            data = data,
            indices = indices,
            create_formula_fn = create_formula_fn,
            predmat = predmat,
            timegroup_cuts = setup$timegroup_cuts,
            median_age_of_dem_timegroup = setup$median_age_of_dem_timegroup,
            short_sleep_geo_mean = setup$short_sleep_geo_mean,
            avg_sleep_geo_mean = setup$avg_sleep_geo_mean,
            empirical = empirical_val
          )
        },
        list(
          setup_target = as.name(setup_target_name),
          data_target = as.name(data_target_name),
          create_formula_fn_name = create_formula_fn_name,
          iter_num = i,
          seed_val = seed_val,
          empirical_val = empirical
        )
      )
    )
  }

  bootstrap_targets
}

# Run bootstrap for a particular model formula and timegroup target
run_bootstrap <- function(
  data,
  create_formula_fn,
  empirical = TRUE,
  intervals = TRUE
) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(
    data,
    mincor = 0,
    exclude = c(
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa",
      "eid"
    )
  )

  ## reference compositions
  all_comp <-
    acomp(data[,
      c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
    ])

  short_sleep_comp <-
    all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
  short_sleep_geo_mean <-
    acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

  avg_sleep_comp <-
    all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
  avg_sleep_geo_mean <-
    acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

  ## Ordinal factor variables to numeric for faster imputation
  data <- ordinal_to_numeric(data)

  ## Constants for dementia risk

  min_age_of_dem <- min(data$age_dem)
  max_age_of_dem <- max(data$age_dem)
  age_range <- max_age_of_dem - min_age_of_dem
  timegroup_steps <- ceiling(age_range * 2)
  median_age_of_dem <- median(data[data$dem == 1, ]$age_dem)

  timegroup_cuts <-
    seq(
      from = min_age_of_dem,
      to = max_age_of_dem,
      length.out = timegroup_steps
    )

  median_age_of_dem_timegroup <-
    which(timegroup_cuts > median_age_of_dem)[1] - 1

  if (isTRUE(intervals)) {
    result <- boot(
      data = data,
      statistic = bootstrap_substitutions_fn,
      create_formula_fn = create_formula_fn,
      predmat = predmat,
      empirical = empirical,
      timegroup_cuts = timegroup_cuts,
      median_age_of_dem_timegroup = median_age_of_dem_timegroup,
      short_sleep_geo_mean = short_sleep_geo_mean,
      avg_sleep_geo_mean = avg_sleep_geo_mean,
      R = bootstrap_iterations,
      parallel = "multicore",
      ncpus = ncpus
    )
  } else {
    result <- bootstrap_substitutions_fn(
      data = data,
      indices = seq(1, nrow(data)),
      create_formula_fn = create_formula_fn,
      predmat = predmat,
      timegroup_cuts = timegroup_cuts,
      median_age_of_dem_timegroup = median_age_of_dem_timegroup,
      short_sleep_geo_mean = short_sleep_geo_mean,
      avg_sleep_geo_mean = avg_sleep_geo_mean,
      empirical = empirical
    )
  }

  return(result)
}

# Ran each bootstrap iteration
# Imputes the data, fits the model for the given formula and
# predicts the risks for a variety of substitutions and cohorts
bootstrap_substitutions_fn <- function(
  data,
  indices,
  create_formula_fn,
  predmat,
  timegroup_cuts,
  median_age_of_dem_timegroup,
  short_sleep_geo_mean,
  avg_sleep_geo_mean,
  empirical = TRUE
) {
  this_sample <- data[indices, ]

  print("Imputing")
  print(format(Sys.time(), "%H:%M:%S"))

  imp <- mice(this_sample, m = 1, predictorMatrix = predmat, maxit = maxit)

  imp <- complete(imp)
  setDT(imp)
  imp[, id := .I]
  imp_len <- nrow(imp)

  # set ordinal variables back to factors

  imp$avg_total_household_income <- as.factor(imp$avg_total_household_income)
  levels(imp$avg_total_household_income) <- c(
    "<18",
    "18-30",
    "31-50",
    "52-100",
    ">100"
  )
  imp$chronotype <- as.factor(imp$chronotype)
  imp$apoe_e4 <- as.factor(imp$apoe_e4)
  imp$highest_qual <- as.factor(imp$highest_qual)
  imp$smok_status <- as.factor(imp$smok_status)
  levels(imp$smok_status) <- c("current", "former", "never")

  # fit model
  print("Fitting model")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  models <- fit_model(imp, timegroup_cuts, create_formula_fn)

  if (isFALSE(empirical)) {
    # shift covariates to match mean (continuous vars) or
    # probability (categorical vars)
    # of Schoeler et al pseudo-pop (see paper)
    imp <-
      imp %>%
      mutate(
        sex = sample(
          c("female", "male"),
          n(),
          replace = TRUE,
          prob = c(0.504, 0.496)
        ),
        retired = rbinom(n(), 1, prob = 0.193),
        avg_total_household_income = sample(
          c("<18", "18-30", "31-50", "52-100", ">100"),
          n(),
          replace = TRUE,
          prob = c(0.264, 0.141, 0.205, 0.145, 0.245)
        ),
        smok_status = sample(
          c("current", "former", "never"),
          n(),
          replace = TRUE,
          prob = c(0.208, 0.359, 0.433)
        )
      )
  }

  # data for g-computation/standardisation
  imp <- imp[rep(seq_len(imp_len), each = median_age_of_dem_timegroup)]
  imp[, timegroup := rep(1:median_age_of_dem_timegroup, imp_len)]

  print("short_sleep_inactive")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  short_sleep_inactive <-
    calc_substitution(
      short_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_inactivity"),
      timegroup = median_age_of_dem_timegroup
    )

  print("short_sleep_light")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  short_sleep_light <-
    calc_substitution(
      short_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_light"),
      timegroup = median_age_of_dem_timegroup
    )

  print("short_sleep_mvpa")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  short_sleep_mvpa <-
    calc_substitution(
      short_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_mvpa"),
      timegroup = median_age_of_dem_timegroup
    )

  print("avg_sleep_inactive")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  avg_sleep_inactive <-
    calc_substitution(
      avg_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_inactivity"),
      timegroup = median_age_of_dem_timegroup
    )

  print("avg_sleep_light")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  avg_sleep_light <-
    calc_substitution(
      avg_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_light"),
      timegroup = median_age_of_dem_timegroup
    )

  print("avg_sleep_mvpa")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  avg_sleep_mvpa <-
    calc_substitution(
      avg_sleep_geo_mean,
      imp,
      models[["model_dem"]],
      models[["model_death"]],
      models[["model_formula"]],
      c("avg_sleep", "avg_mvpa"),
      timegroup = median_age_of_dem_timegroup
    )

  full_df <- full_join(
    short_sleep_inactive,
    short_sleep_light,
    by = "offset"
  ) |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset")

  print("Returning results")
  print(format(Sys.time(), "%H:%M:%S"))
  return(as.matrix(full_df))
}

produce_plots <- function(primary, s1, s2, s3) {
  primary_plot <- process_dem_output(primary)
  s1_plot <- process_dem_output(s1)
  s2_plot <- process_dem_output(s2)
  s3_plot <- process_dem_output(s3, intervals = FALSE)

  save_plot(
    primary_plot[[1]],
    file.path(
      data_dir,
      "../Manuscript/Main_figures/Figure 1.svg"
    )
  )
  save_plot(
    s1_plot[[1]],
    file.path(
      data_dir,
      "../Manuscript/Appendix_figures/Substitutions_s1.svg"
    )
  )
  save_plot(
    s2_plot[[1]],
    file.path(
      data_dir,
      "../Manuscript/Appendix_figures/Substitutions_s2.svg"
    )
  )
  save_plot(
    s3_plot[[1]],
    file.path(
      data_dir,
      "../Manuscript/Appendix_figures/Substitutions_s3.svg"
    )
  )
}

process_dem_output <- function(result, intervals = TRUE) {
  data <- result
  if (isTRUE(intervals)) {
    full_sample_results <- data$t0
  } else {
    full_sample_results <- data
  }

  plot_data <-
    pivot_longer(
      as.data.frame(full_sample_results),
      -offset,
      values_to = "risk",
      names_to = "Substitution"
    ) |>
    group_by(Substitution) |>
    mutate(ref_risk = ifelse(offset == 0, risk, NA_real_)) |>
    fill(ref_risk, .direction = "downup") |>
    mutate(
      risk_dif = risk - ref_risk,
      risk_ratio = risk / ref_risk
    ) |>
    mutate(
      Reference = ifelse(
        str_detect(Substitution, "short_sleep"),
        "Short sleepers",
        "Normal sleepers"
      )
    ) |>
    mutate(Substitution = str_remove(Substitution, "_short_sleep_geo_mean")) |>
    mutate(Substitution = str_remove(Substitution, "_avg_sleep_geo_mean")) |>
    mutate(
      Substitution = ifelse(
        Substitution == "avg_inactivity",
        "Inactivity",
        ifelse(
          Substitution == "avg_light",
          "Light activity",
          "MVPA"
        )
      )
    )

  if (isTRUE(intervals)) {
    num_subs <- sub_steps * 2 + 1
    sub_col_names <- data$t[1, 1:num_subs]
    get_quantiles <- function(start, substitution, reference) {
      slice <- data$t[, start:(start + num_subs - 1)]
      middle_col <- ceiling(ncol(slice) / 2)

      slice <- t(apply(
        slice,
        1,
        function(row, idx) {
          zero_offset <- row[middle_col]
          row / zero_offset
        },
        idx = middle_col
      ))

      quantiles <-
        as.data.frame(t(apply(
          slice,
          2,
          function(column) quantile(column, probs = c(0.025, 0.975))
        )))
      colnames(quantiles) <- c("lower", "upper")
      quantiles$Substitution <- substitution
      quantiles$Reference <- reference
      quantiles$offset <- sub_col_names
      quantiles
    }

    sub_start <- num_subs + 1

    inactivity_short_sleep <-
      get_quantiles(sub_start, "Inactivity", "Short sleepers")
    sub_start <- sub_start + num_subs

    light_short_sleep <-
      get_quantiles(sub_start, "Light activity", "Short sleepers")
    sub_start <- sub_start + num_subs

    mvpa_short_sleep <-
      get_quantiles(sub_start, "MVPA", "Short sleepers")
    sub_start <- sub_start + num_subs

    inactivity_avg_sleep <-
      get_quantiles(sub_start, "Inactivity", "Normal sleepers")
    sub_start <- sub_start + num_subs

    light_avg_sleep <-
      get_quantiles(sub_start, "Light activity", "Normal sleepers")
    sub_start <- sub_start + num_subs

    mvpa_avg_sleep <-
      get_quantiles(sub_start, "MVPA", "Normal sleepers")

    all_quantiles <- rbind(
      inactivity_short_sleep,
      light_short_sleep,
      mvpa_short_sleep,
      inactivity_avg_sleep,
      light_avg_sleep,
      mvpa_avg_sleep
    )

    plot_data <-
      full_join(
        plot_data,
        all_quantiles,
        by = c("offset", "Substitution", "Reference")
      )
  }

  rr_plot <- function(sub, refcomp, colour) {
    plot_data$Substitution <-
      ifelse(
        plot_data$Substitution == "Inactivity",
        "inactivity",
        ifelse(
          plot_data$Substitution == "Light activity",
          "light activity",
          "MVPA"
        )
      )

    plot_data2 <-
      plot_data[
        plot_data$Substitution == sub & plot_data$Reference == refcomp,
      ]

    p <- plot_data2 |>
      ggplot(aes(x = offset, y = risk_ratio)) +
      geom_line(colour = colour) +
      facet_wrap(~Substitution, nrow = 2) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      xlab("") +
      ylab("Risk ratio") +
      annotate(
        geom = "text",
        x = 0,
        y = -0.15,
        hjust = 0.5,
        fontface = 1,
        size = 14 / .pt,
        label = "Minutes",
        family = "serif"
      ) +
      annotate(
        geom = "text",
        x = -20,
        y = -0.5,
        hjust = 1,
        fontface = 1,
        size = 12 / .pt,
        label = "Less sleep",
        family = "serif"
      ) +
      annotate(
        geom = "text",
        x = 20,
        y = -0.5,
        hjust = 0,
        fontface = 1,
        size = 12 / .pt,
        label = "More sleep",
        family = "serif"
      ) +
      geom_segment(
        aes(
          x = 1,
          y = -0.625,
          xend = 15,
          yend = -0.625
        ),
        arrow = arrow(length = unit(0.15, "cm"))
      ) +
      geom_segment(
        aes(
          x = -1,
          y = -0.625,
          xend = -15,
          yend = -0.625
        ),
        arrow = arrow(length = unit(0.15, "cm"))
      ) +
      annotate(
        geom = "text",
        x = -20,
        y = -0.75,
        hjust = 1,
        size = 12 / .pt,
        label = paste("More", plot_data2$Substitution),
        family = "serif",
        fontface = 1,
        size = 12 / .pt
      ) +
      annotate(
        geom = "text",
        x = 20,
        y = -0.75,
        hjust = 0,
        label = paste("Less", plot_data2$Substitution),
        family = "serif",
        fontface = 1,
        size = 12 / .pt
      ) +
      coord_cartesian(ylim = c(0.33, 3), expand = FALSE, clip = "off") +
      cowplot::theme_cowplot(
        font_size = 12,
        font_family = "serif",
        line_size = 0.25
      ) +
      theme(
        plot.margin = unit(c(1, 1, 4, 1), "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#585656"),
        panel.grid = element_line(colour = "grey92"),
        panel.grid.minor = element_line(linewidth = rel(0.5)),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color = "#585656")
      )
    if (isTRUE(intervals)) {
      p <- p +
        geom_ribbon(
          aes(ymin = lower, ymax = upper),
          alpha = 0.25,
          fill = colour
        )
    }
    p
  }

  # normal sleepers
  p1 <- rr_plot("inactivity", "Normal sleepers", "#fc020f")
  p2 <- rr_plot("light activity", "Normal sleepers", "#145e01")
  p3 <- rr_plot("MVPA", "Normal sleepers", "#011869")

  pnorm <-
    plot_grid(
      NULL,
      p1 +
        labs(x = "", title = "A") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      p2 +
        labs(x = "", title = "B") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      p3 +
        labs(x = "", title = "C") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      align = "vh",
      rel_heights = c(0.05, 1, 1, 1),
      nrow = 4,
      labels = "Normal sleepers",
      label_fontfamily = "serif",
      hjust = -1.1
    )

  # short sleepers
  p4 <- rr_plot("inactivity", "Short sleepers", "#ff747b")
  p5 <- rr_plot("light activity", "Short sleepers", "#6ed853")
  p6 <- rr_plot("MVPA", "Short sleepers", "#708ff9")

  pshort <-
    plot_grid(
      NULL,
      p4 +
        labs(title = "D", y = "") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      p5 +
        labs(y = "", title = "E") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      p6 +
        labs(y = "", title = "F") +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
      align = "vh",
      rel_heights = c(0.05, 1, 1, 1),
      nrow = 4,
      labels = "Short sleepers",
      label_fontfamily = "serif",
      hjust = -1.3
    )

  plot <- plot_grid(pnorm, pshort, nrow = 1)

  return(list(plot, plot_data))
}

apply_substitution <- function(
  df,
  from_var,
  to_var,
  duration,
  timegroup_cuts
) {
  lapply(
    df,
    function(.x) {
      # 1) compute 1st and 99th quantiles of the 'from' and 'to' cols
      from_q <- quantile(.x[[from_var]], probs = c(0.01, 0.99), na.rm = TRUE)
      to_q <- quantile(.x[[to_var]], probs = c(0.01, 0.99), na.rm = TRUE)

      min_from <- from_q[1]
      max_from <- from_q[2]
      min_to <- to_q[1]
      max_to <- to_q[2]

      sub_df <- .x |>
        mutate(
          new_from = .data[[from_var]] - duration,
          # clamp it
          "{from_var}" := pmin(pmax(new_from, min_from), max_from),
          new_to = .data[[to_var]] + duration,
          "{to_var}" := pmin(pmax(new_to, min_to), max_to),
          sub_name = paste0(from_var, "_", to_var, "_", duration),
          censoring = 1L # no (meaningful) censoring
        ) |>
        select(-new_from, -new_to)

      # 2) composition -> ILR
      comp <- acomp(sub_df[, c(
        "avg_sleep",
        "avg_inactivity",
        "avg_light",
        "avg_mvpa"
      )])
      ilr_vars <- ilr(comp, V = v) |>
        setNames(c("R1", "R2", "R3"))

      sub_df[, c("R1", "R2", "R3")] <- as.data.table(ilr_vars)
      sub_df
    }
  )
}

make_cuts <- function(df) {
  max_follow_up <- max(df$time_to_dem)
  timegroup_steps <- ceiling(max_follow_up / 365)

  timegroup_cuts <-
    seq(
      from = 0,
      to = max_follow_up,
      length.out = timegroup_steps + 1
    )
  timegroup_cuts
}

estimate_lmtp_reference <- function(df, baseline_covars) {
  lapply(df, function(.x) {
    cens <- grep("^censoring_", names(.x), value = TRUE)
    trt <- c("R1", "R2", "R3")
    outcomes <- grep("^dem_", names(.x), value = TRUE)
    compete <- grep("^death_", names(.x), value = TRUE)
    .x <- select(.x, all_of(c(baseline_covars, trt, cens, outcomes, compete)))
    trt <- list(c("R1", "R2", "R3"))
    lmtp::lmtp_survival(
      data = .x,
      trt = trt,
      outcomes = outcomes,
      cens = cens,
      baseline = baseline_covars,
      compete = compete,
      folds = 1,
      learners_outcome = "SL.glm.Q",
      learners_trt = "SL.glm.g",
      control = lmtp::lmtp_control(
        .learners_outcome_folds = 2,
        .learners_trt_folds = 2
      ),
      mtp = TRUE
    )
  })
}

estimate_lmtp_subs <- function(df, sub_df, baseline_covars) {
  pmap(
    c(df, unlist(sub_df, recursive = FALSE)),
    function(.x, .y) {
      cens <- grep("^censoring_", names(.x), value = TRUE)
      trt <- c("R1", "R2", "R3")
      outcomes <- grep("^dem_", names(.x), value = TRUE)
      compete <- grep("^death_", names(.x), value = TRUE)
      .x <- select(.x, all_of(c(baseline_covars, cens, trt, outcomes, compete)))
      .y <- select(.y, all_of(c(baseline_covars, cens, trt, outcomes, compete)))
      trt <- list(c("R1", "R2", "R3"))
      lmtp::lmtp_survival(
        .x,
        trt = trt,
        outcomes = outcomes,
        cens = cens,
        baseline = baseline_covars,
        compete = compete,
        shifted = .y,
        folds = 1,
        control = lmtp::lmtp_control(
          .learners_outcome_folds = 2,
          .learners_trt_folds = 2
        ),
        mtp = TRUE
      )
    }
  )
}
