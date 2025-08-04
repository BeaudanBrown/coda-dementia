cum_test_targets <- list(
  #### SYNTH COMPS ####
  tar_rep(
    cum_boot_indices,
    sample(seq_len(nrow(df)), size = floor(0.5 * nrow(df))),
    batches = n_boots / 2
  ),
  tar_target(
    cum_boot_train,
    df[cum_boot_indices],
    pattern = map(cum_boot_indices)
  ),
  tar_target(
    cum_boot_df_test,
    df[-cum_boot_indices],
    pattern = map(cum_boot_indices)
  ),
  tar_target(
    cum_boot_imp_train,
    impute_data(cum_boot_train, m, maxit),
    pattern = map(cum_boot_train)
  ),
  tar_rep(
    cum_boot_sample_test,
    bootstrap_sample(cum_boot_df_test),
    batches = n_boots
  ),
  tar_target(
    cum_boot_imp_test,
    impute_data(cum_boot_sample_test, m, maxit),
    pattern = map(cum_boot_sample_test)
  ),
  tar_target(
    cum_boot_train_models,
    train_model(cum_boot_imp_train, timegroup_cuts, get_primary_formula),
    pattern = map(cum_boot_imp_train)
  ),
  tar_target(
    cum_boot_test_models,
    train_model(cum_boot_imp_test, timegroup_cuts, get_primary_formula),
    pattern = map(cum_boot_imp_test)
  ),
  tar_target(
    cum_boot_covar_data,
    get_covar_data(cum_boot_imp_train),
    pattern = map(cum_boot_imp_train)
  ),
  tar_target(
    cum_boot_synth_comp_risk,
    get_synth_risk(
      cum_boot_covar_data,
      synth_comps_filtered,
      cum_boot_train_models,
      final_time
    ),
    pattern = map(cum_boot_covar_data, cum_boot_train_models)
  ),
  tar_target(
    cum_boot_ref_comps,
    {
      list(
        Best = cum_boot_synth_comp_risk[order(risk), ][1],
        Worst = cum_boot_synth_comp_risk[order(-risk), ][1],
        Typical = cum_boot_synth_comp_risk[order(-dens), ][1]
      )
    },
    pattern = map(cum_boot_synth_comp_risk)
  )
  # tar_map(
  #   values = list(
  #     comp = c(
  #       "Best",
  #       "Typical",
  #       "Worst"
  #     )
  #   ),
  #   names = comp,
  #   tar_target(
  #     cum_risks,
  #     {
  #       cum_boot_imp_test[, c("R1", "R2", "R3")] <-
  #         cum_boot_ref_comps[[comp]][, c("R1", "R2", "R3")]
  #       risks <- get_risk(cum_boot_imp_test, cum_boot_test_models, final_time)
  #       risks <- risks[, .(eid, risk, timegroup)]
  #       risks[,
  #         .(
  #           risk = mean(risk, na.rm = TRUE),
  #           B = unique(cum_boot_imp_test$tar_batch)
  #         ),
  #         by = .(timegroup)
  #       ]
  #     },
  #     pattern = map(cum_boot_imp_test, cum_boot_test_models)
  #   ),
  #   tar_target(
  #     cum_avg_risks,
  #     {
  #       cum_risks[,
  #         .(
  #           risk = mean(risk),
  #           lower_risk = quantile(risk, 0.025),
  #           upper_risk = quantile(risk, 0.975),
  #           Composition = comp
  #         ),
  #         by = timegroup
  #       ]
  #     }
  #   )
  # ),
  # tar_target(
  #   cum_plot_data,
  #   {
  #     bind_rows(
  #       cum_avg_risks_Worst,
  #       cum_avg_risks_Typical,
  #       cum_avg_risks_Best
  #     )
  #   }
  # ),
  # tar_target(
  #   cum_plot,
  #   {
  #     composition_colors <- c(
  #       "Worst" = "#ff747b",
  #       "Typical" = "#708ff9",
  #       "Best" = "#6ed853"
  #     )
  #     cum_plot_data |>
  #       ggplot(aes(x = timegroup, y = risk)) +
  #       geom_ribbon(
  #         aes(ymin = lower_risk, ymax = upper_risk, fill = Composition),
  #         alpha = 0.25
  #       ) +
  #       geom_line(aes(colour = Composition)) +
  #       labs(
  #         x = "Time since baseline (years)",
  #         y = "Cumulative all-cause dementia incidence"
  #       ) +
  #       scale_color_manual(
  #         values = composition_colors
  #       ) +
  #       scale_fill_manual(
  #         values = composition_colors
  #       ) +
  #       cowplot::theme_cowplot(
  #         font_size = 16,
  #         font_family = "serif",
  #         line_size = 0.25
  #       ) +
  #       theme(
  #         panel.border = element_rect(fill = NA, colour = "#585656"),
  #         panel.grid = element_line(colour = "grey92"),
  #         panel.grid.minor = element_line(linewidth = rel(0.5)),
  #         axis.ticks.y = element_blank(),
  #         axis.line = element_line(color = "#585656"),
  #         axis.title.x = element_text(family = "serif", size = 20),
  #         axis.title.y = element_text(family = "serif", size = 20)
  #       )
  #   }
  # )
)
