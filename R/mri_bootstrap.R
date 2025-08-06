get_mri_model <- function(imp, outcome) {
  # Fit linear regression
  model_formula <- get_mri_formula(imp)
  model_formula <- update(model_formula, as.formula(paste(outcome, "~ .")))
  model <- lm(model_formula, imp)
  strip_lm(model)
}

get_mri_ref <- function(imp, outcome, model) {
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
  subbed <- apply_substitution(imp, from_var, to_var, duration, comp_limits)
  subbed_results <- subbed$results[[1]]
  subbed_results$estimate <- predict(model, newdata = subbed_results)
  subbed$results <- list(subbed_results[, .(eid, estimate, substituted)])
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
        "Change in Total brain volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "wmv") {
      ylabel <- expression(paste(
        "Change in White matter volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "gmv") {
      ylabel <- expression(paste(
        "Change in Grey matter volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "hip") {
      ylabel <- expression(paste(
        "Change in Hippocampal volume ",
        (cm^{
          "3"
        })
      ))
    } else if (outcome == "log_wmh") {
      ylabel <- "Change in Log WMH"
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

make_mri_plots <- function(mri_results, cohort) {
  sub_types <- c("avg_mvpa", "avg_light", "avg_inactivity")
  colours <- c("#708ff9", "#6ed853", "#ff747b")

  left_centre <- unit(0.4, "npc")
  right_centre <- unit(0.6, "npc")

  minutes_offset <- unit(-0.8, "cm")
  sub_offset <- unit(-1.2, "cm")
  sleep_offset <- unit(-1.8, "cm")
  arrow_offset <- unit(-1.5, "cm")

  left_arrow_xs <- unit(c(0.4, 0.2), "npc")
  right_arrow_xs <- unit(c(0.6, 0.8), "npc")
  arrow_ys <- unit(c(-1.5, -1.5), "cm")

  limit <- max(abs(min(mri_results$lower_md)), abs(max(mri_results$upper_md)))
  lower_lim <- -limit
  upper_lim <- limit

  unique_cohorts <- unique(mri_results$cohort)

  rbindlist(Map(
    function(sub_type, colour) {
      sub_results <- mri_results |>
        filter(from_var == sub_type) |>
        filter(prop_substituted > intervention_threshold)

      lapply(unique_cohorts, function(cohort_name) {
        labels <- get_mri_labels(sub_results)
        cohort_data <- sub_results |>
          filter(cohort == cohort_name) |>
          # Add in the pretend zero
          bind_rows(data.frame(
            from_var = sub_type,
            to_var = "avg_sleep",
            duration = 0,
            prop_substituted = 1,
            mean_sub_estimate = 1,
            mean_ref_estimate = 1,
            md = 0,
            lower_md = 0,
            upper_md = 0,
            cohort = labels$cohort,
            outcome = labels$outcome
          ))
        dark_factor <- switch(
          cohort_name,
          short_sleeper = 1.2,
          long_sleeper = 1.2,
          0.7
        )
        colour <- adjust_colour(colour, dark_factor)

        p <- cohort_data |>
          ggplot(aes(x = duration, y = md)) +
          geom_line(colour = colour) +
          geom_hline(yintercept = 0, linetype = "dotted") +
          geom_ribbon(
            aes(ymin = lower_md, ymax = upper_md),
            alpha = 0.2,
            fill = colour
          ) +
          labs(x = "", y = labels$ylabel) +
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
          cohort = cohort_name,
          outcome = labels$outcome
        )
      }) |>
        rbindlist()
    },
    sub_types,
    colours
  ))
}

get_mri_synth_labels <- function(mri_results) {
  outcome <- unique(mri_results$outcome)
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
    outcome = outcome
  )
}

make_mri_synth_plot <- function(synth_results) {
  labels <- get_mri_synth_labels(synth_results)
  synth_results <- synth_results |>
    mutate(comp = factor(comp, levels = c("Worst", "Typical", "Best")))
  synth_results$comp <- recode(synth_results$comp, "Best" = "Ideal")
  ggplot(synth_results, aes(x = comp, y = estimate)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    labs(
      x = "Composition",
      y = labels$ylabel,
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
}
