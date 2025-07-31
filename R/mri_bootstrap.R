get_mri_model <- function(imp, outcome) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)

  # Fit linear regression
  model_formula <- get_mri_formula(imp)
  model_formula <- update(model_formula, as.formula(paste(outcome, "~ .")))
  model <- lm(model_formula, imp)
  strip_lm(model)
}

get_mri_ref <- function(imp, outcome, model) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  imp$estimate <- predict(model, newdata = imp)
  data.table(
    outcome = outcome,
    results = list(result = imp |> select(eid, estimate)),
    B = unique(imp$tar_batch)
  )
}

get_mri_subs <- function(
  imp,
  outcome,
  model,
  from_var,
  to_var,
  duration,
  comp_limits
) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  subbed <- apply_substitution(imp, from_var, to_var, duration, comp_limits)
  subbed_results <- subbed$results[[1]]
  subbed_results$estimate <- predict(model, newdata = subbed_results)
  subbed$results <- list(subbed_results[, .(eid, estimate)])
  subbed[, outcome := outcome]
  subbed
}

merge_estimates <- function(sub_estimates, ref_estimates, outcome, cohort) {
  sub_estimates |>
    rename("sub_estimate" = "results") |>
    left_join(
      ref_estimates |> rename("ref_estimate" = "results"),
      by = "B"
    ) |>
    mutate(
      md = sub_estimate - ref_estimate
    ) |>
    dplyr::group_by(from_var, to_var, duration) |>
    dplyr::summarize(
      cohort = cohort,
      outcome = outcome,
      prop_substituted = mean(prop_substituted, na.rm = TRUE),
      mean_sub_estimate = mean(sub_estimate, na.rm = TRUE),
      mean_ref_estimate = mean(ref_estimate, na.rm = TRUE),
      lower_md = quantile(md, 0.025, names = FALSE),
      upper_md = quantile(md, 0.975, names = FALSE),
      md = mean(md, na.rm = TRUE),
      .groups = "drop"
    )
}

get_mri_labels <- function(mri_results) {
  outcome <- unique(mri_results$outcome)
  cohort <- unique(mri_results$cohort)
  list(
    ylabel = if (outcome == "tbv") {
      ylabel <- expression(paste(
        "Total brain volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "wmv") {
      ylabel <- expression(paste(
        "White matter volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "gmv") {
      ylabel <- expression(paste(
        "Grey matter volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "hip") {
      ylabel <- expression(paste(
        "Hippocampal volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "log_wmh") {
      ylabel <- "Log WMH"
    } else {
      ylabel <- expression(as.character(outcome))
    },
    sub_name = unique(case_when(
      mri_results$from_var == "avg_inactivity" ~ "Inactivity",
      mri_results$from_var == "avg_light" ~ "Light Activity",
      mri_results$from_var == "avg_mvpa" ~ "MVPA",
      .default = as.character(mri_results$from_var)
    )),
    outcome = outcome,
    cohort = cohort
  )
}

make_mri_plots <- function(mri_results, colour) {
  sub_types <- c("avg_mvpa", "avg_light", "avg_inactivity")

  y_mean <- mean(mri_results$md, na.rm = TRUE)
  y_sd <- sd(mri_results$md, na.rm = TRUE)
  limit <- max(abs(min(mri_results$lower_md)), abs(max(mri_results$upper_md)))
  lower_lim <- -limit
  upper_lim <- limit

  rbindlist(lapply(sub_types, function(sub_type) {
    sub_results <- mri_results |>
      filter(from_var == sub_type) |>
      filter(prop_substituted > 0.8)
    labels <- get_mri_labels(sub_results)

    left_centre <- unit(0.4, "npc")
    right_centre <- unit(0.6, "npc")

    minutes_offset <- unit(-0.8, "cm")
    sub_offset <- unit(-1.2, "cm")
    sleep_offset <- unit(-1.8, "cm")
    arrow_offset <- unit(-1.5, "cm")

    left_arrow_xs <- unit(c(0.4, 0.2), "npc")
    right_arrow_xs <- unit(c(0.6, 0.8), "npc")
    arrow_ys <- unit(c(-1.5, -1.5), "cm")

    p <- sub_results |>
      ggplot(aes(x = duration, y = md)) +
      geom_line(colour = colour) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_ribbon(
        aes(ymin = lower_md, ymax = upper_md),
        alpha = 0.2,
        fill = colour
      ) +
      labs(x = "", y = paste("Change in", labels$ylabel)) +
      annotation_custom(
        textGrob(
          "Minutes",
          gp = gpar(fontsize = 14, fontfamily = "serif", fontface = 1),
          x = unit(0.5, "npc"),
          y = minutes_offset,
          hjust = 0.5
        )
      ) +
      annotation_custom(
        textGrob(
          "Less sleep",
          gp = gpar(fontsize = 12, fontfamily = "serif", fontface = 1),
          x = left_centre,
          y = sleep_offset,
          hjust = 1
        )
      ) +
      annotation_custom(
        textGrob(
          "More sleep",
          gp = gpar(fontsize = 12, fontfamily = "serif", fontface = 1),
          x = right_centre,
          y = sleep_offset,
          hjust = 0
        )
      ) +
      annotation_custom(
        textGrob(
          paste("More", labels$sub_name),
          gp = gpar(fontsize = 12, fontfamily = "serif", fontface = 1),
          x = left_centre,
          y = sub_offset,
          hjust = 1
        )
      ) +
      annotation_custom(
        textGrob(
          paste("Less", labels$sub_name),
          gp = gpar(fontsize = 12, fontfamily = "serif", fontface = 1),
          x = right_centre,
          y = sub_offset,
          hjust = 0
        )
      ) +

      # Arrows
      annotation_custom(
        grid::linesGrob(
          x = right_arrow_xs,
          y = arrow_ys,
          arrow = grid::arrow(
            angle = 30,
            length = unit(0.15, "cm"),
            ends = "last"
          )
        )
      ) +
      annotation_custom(
        grid::linesGrob(
          x = left_arrow_xs,
          y = arrow_ys,
          arrow = grid::arrow(
            angle = 30,
            length = unit(0.15, "cm"),
            ends = "last"
          )
        )
      ) +

      # Coordinate system and theme
      coord_cartesian(
        xlim = c(-60, 60),
        ylim = c(lower_lim, upper_lim),
        expand = FALSE,
        clip = "off"
      ) +
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
        axis.line = element_line(color = "#585656"),
        axis.title.y = element_text(margin = margin(r = 10))
      )
    data.table(
      sub_type = sub_type,
      plot = list(p),
      cohort = labels$cohort,
      outcome = labels$outcome
    )
  }))
}

make_plot_grid <- function(plot_data, cohort_order, subtype_order) {
  plot_data[, cohort := factor(cohort, levels = cohort_order)]
  plot_data[, sub_type := factor(sub_type, levels = subtype_order)]

  # Reshape for grid
  grid_list <- lapply(cohort_order, function(this_cohort) {
    plots_row <- plot_data[cohort == this_cohort][order(sub_type), plot]
    wrap_plots(plots_row, nrow = 1)
  })

  # Assemble rows into a grid
  final_plot <- wrap_plots(grid_list, ncol = 1)
  final_plot
}
