cum_targets <- list(
  #### SYNTH COMPS ####
  tar_target(synth_comps, generate_compositions(df)),
  tar_target(synth_comps_dens, add_density(df, synth_comps, synth_threshold)),
  tar_target(
    synth_comps_filtered,
    synth_comps_dens[dens > dens_threshold & !is.na(dens), ]
  ),
  tar_target(
    train_indices,
    sample(seq_len(nrow(df)), size = floor(0.5 * nrow(df)))
  ),
  tar_target(df_train, df[train_indices]),
  tar_target(df_test, df[-train_indices]),
  tar_target(imp_train, impute_data(df_train, m, maxit)),
  tar_rep(boots_test, bootstrap_sample(df_test), batches = n_boots),
  tar_target(
    imp_test,
    impute_data(boots_test, m, maxit),
    pattern = map(boots_test)
  ),
  tar_target(
    train_models,
    train_model(imp_train, timegroup_cuts, get_primary_formula)
  ),
  tar_target(
    test_models,
    train_model(imp_test, timegroup_cuts, get_primary_formula),
    pattern = map(imp_test)
  ),
  tar_target(covar_data, get_covar_data(imp_train)),
  tar_target(
    synth_comp_risk,
    get_synth_risk(
      covar_data,
      synth_comps_filtered,
      train_models,
      final_time
    )
  ),
  tar_target(
    reference_comps,
    {
      list(
        Best = synth_comp_risk[order(risk), ][1],
        Worst = synth_comp_risk[order(-risk), ][1],
        Typical = synth_comp_risk[order(-dens), ][1]
      )
    }
  ),
  tar_target(
    old_reference_comps,
    {
      best_vars <- data.table(
        avg_sleep = 7.25,
        avg_inactivity = 11.5,
        avg_light = 2,
        avg_mvpa = 3.25
      ) *
        60
      old_best <- acomp(best_vars)

      typical_vars <- data.table(
        avg_sleep = 6.5,
        avg_inactivity = 12,
        avg_light = 3.5,
        avg_mvpa = 2
      ) *
        60
      old_typical <- acomp(typical_vars)
      worst_vars <- data.table(
        avg_sleep = 4.5,
        avg_inactivity = 15.25,
        avg_light = 3.75,
        avg_mvpa = 0.5
      ) *
        60
      old_worst <- acomp(worst_vars)

      best_ilr <-
        ilr(old_best, V = v) |>
        setNames(c("R1", "R2", "R3")) |>
        t() |>
        as.data.table()

      typical_ilr <-
        ilr(old_typical, V = v) |>
        setNames(c("R1", "R2", "R3")) |>
        t() |>
        as.data.table()

      worst_ilr <-
        ilr(old_worst, V = v) |>
        setNames(c("R1", "R2", "R3")) |>
        t() |>
        as.data.table()

      best_comp <- cbind(best_vars, best_ilr)
      worst_comp <- cbind(worst_vars, worst_ilr)
      typical_comp <- cbind(typical_vars, typical_ilr)

      list(
        Best = best_comp,
        Worst = worst_comp,
        Typical = typical_comp
      )
    }
  ),

  tar_map(
    values = list(
      comp = c(
        "Best",
        "Typical",
        "Worst"
      )
    ),
    names = comp,
    tar_target(
      old_cum_risks,
      {
        imp_test[, c("R1", "R2", "R3")] <-
          old_reference_comps[[comp]][, c("R1", "R2", "R3")]
        risks <- get_risk(imp_test, test_models, final_time)
        risks <- risks[, .(eid, risk, timegroup)]
        risks[,
          .(
            risk = mean(risk, na.rm = TRUE),
            B = unique(imp_test$tar_batch)
          ),
          by = .(timegroup)
        ]
      },
      pattern = map(imp_test, test_models)
    ),
    tar_target(
      old_cum_avg_risks,
      {
        old_cum_risks[,
          .(
            risk = mean(risk),
            lower_risk = quantile(risk, 0.025),
            upper_risk = quantile(risk, 0.975),
            Composition = comp
          ),
          by = timegroup
        ]
      }
    )
  ),

  tar_target(
    old_cum_plot_data,
    {
      bind_rows(
        old_cum_avg_risks_Worst,
        old_cum_avg_risks_Typical,
        old_cum_avg_risks_Best
      )
    }
  ),
  tar_target(
    old_cum_plot,
    {
      composition_colors <- c(
        "Worst" = "#ff747b",
        "Typical" = "#708ff9",
        "Best" = "#6ed853"
      )
      old_cum_plot_data |>
        ggplot(aes(x = timegroup, y = risk)) +
        geom_ribbon(
          aes(ymin = lower_risk, ymax = upper_risk, fill = Composition),
          alpha = 0.25
        ) +
        geom_line(aes(colour = Composition)) +
        labs(
          x = "Time since baseline (years)",
          y = "Cumulative all-cause dementia incidence"
        ) +
        scale_color_manual(
          values = composition_colors
        ) +
        scale_fill_manual(
          values = composition_colors
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
    }
  ),

  tar_map(
    values = list(
      comp = c(
        "Best",
        "Typical",
        "Worst"
      )
    ),
    names = comp,
    tar_target(
      cum_risks,
      {
        imp_test[, c("R1", "R2", "R3")] <-
          reference_comps[[comp]][, c("R1", "R2", "R3")]
        risks <- get_risk(imp_test, test_models, final_time)
        risks <- risks[, .(eid, risk, timegroup)]
        risks[,
          .(
            risk = mean(risk, na.rm = TRUE),
            B = unique(imp_test$tar_batch)
          ),
          by = .(timegroup)
        ]
      },
      pattern = map(imp_test, test_models)
    ),
    tar_target(
      cum_avg_risks,
      {
        cum_risks[,
          .(
            risk = mean(risk),
            lower_risk = quantile(risk, 0.025),
            upper_risk = quantile(risk, 0.975),
            Composition = comp
          ),
          by = timegroup
        ]
      }
    )
  ),
  tar_target(
    cum_plot_data,
    {
      bind_rows(
        cum_avg_risks_Worst,
        cum_avg_risks_Typical,
        cum_avg_risks_Best
      )
    }
  ),
  tar_target(
    cum_plot,
    {
      composition_colors <- c(
        "Worst" = "#ff747b",
        "Typical" = "#708ff9",
        "Best" = "#6ed853"
      )
      cum_plot_data |>
        ggplot(aes(x = timegroup, y = risk)) +
        geom_ribbon(
          aes(ymin = lower_risk, ymax = upper_risk, fill = Composition),
          alpha = 0.25
        ) +
        geom_line(aes(colour = Composition)) +
        labs(
          x = "Time since baseline (years)",
          y = "Cumulative all-cause dementia incidence"
        ) +
        scale_color_manual(
          values = composition_colors
        ) +
        scale_fill_manual(
          values = composition_colors
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
    }
  )
)
