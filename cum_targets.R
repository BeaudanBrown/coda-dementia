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
    ),
    tar_target(
      cum_pie,
      {
        mins <- as.numeric(reference_comps[[comp]][, c(
          avg_sleep,
          avg_inactivity,
          avg_light,
          avg_mvpa
        )])
        data <- data.table(
          activity = c("Sleep", "Inactivity", "Light Activity", "MVPA"),
          duration_min = mins
        ) |>
          arrange(desc(activity)) |>
          mutate(
            prop = duration_min / sum(duration_min) * 100,
            ypos = cumsum(prop) - 0.5 * prop,
            hours_label = sprintf("%.2f h", duration_min / 60)
          )

        ggplot(data, aes(y = prop, fill = activity)) +
          geom_bar(aes(x = 2), stat = "identity", width = 1, color = "white") +
          coord_polar(theta = "y", start = 0, clip = "off") +
          theme_void() +
          # theme(legend.position = "none", plot.margin = margin(-5, 0, -5, 0)) +
          theme(
            legend.position = "none",
            plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
          ) +
          geom_text(
            aes(x = 2.2, y = ypos, label = hours_label),
            color = "white",
            size = 4
          ) +
          ggrepel::geom_text_repel(
            aes(x = 2.5, y = ypos, label = activity),
            direction = "x",
            nudge_x = 0.5,
            hjust = 0,
            segment.size = 0.3,
            size = 4,
            box.padding = 0.3,
            point.padding = 0,
            min.segment.length = 0
          ) +
          scale_fill_brewer(palette = "Set1")
      }
    )
  ),
  tar_target(
    cum_pie_combined,
    {
      tar_make(c(cum_pie_Worst, cum_pie_Typical, cum_pie_Best))
      tar_load(c(cum_pie_Worst, cum_pie_Typical, cum_pie_Best))
      tar_read(cum_pie_Worst)
      pies <- list(cum_pie_Worst, cum_pie_Typical, cum_pie_Best)
      wrap_plots(pies, nrow = 3, axis_titles = "collect")
    }
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
      plot <- cum_plot_data |>
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
          font_size = 12,
          font_family = "serif",
          line_size = 0.25
        ) +
        theme(
          panel.border = element_rect(fill = NA, colour = "#585656"),
          panel.grid = element_line(colour = "grey92"),
          panel.grid.minor = element_line(linewidth = rel(0.5)),
          plot.background = element_rect(fill = "white", colour = NA),
          axis.ticks.y = element_blank(),
          axis.line = element_line(color = "#585656")
        )
      list(
        plot_grid = plot,
        n_cohorts = 1,
        n_rows = 1
      )
    }
  )
)
