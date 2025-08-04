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

get_risk <- function(imp, models, final_time) {
  imp_len <- nrow(imp)
  imp_long_cuts <- imp[rep(seq_len(imp_len), each = final_time)]
  imp_long_cuts[, timegroup := rep(1:final_time, imp_len)]

  imp_long_cuts[,
    haz_dem := predict(
      models[["model_dem"]],
      newdata = .SD,
      type = "response"
    )
  ]
  imp_long_cuts[,
    haz_death := predict(
      models[["model_death"]],
      newdata = .SD,
      type = "response"
    )
  ]
  setkey(imp_long_cuts, id, timegroup) # sort and set keys
  imp_long_cuts[,
    risk := cumsum(
      haz_dem * cumprod((1 - lag(haz_dem, default = 0)) * (1 - haz_death))
    ),
    by = id
  ]
  imp_long_cuts[,
    .(eid, risk, timegroup)
  ]
  imp_long_cuts
}

get_ref_risk <- function(imp, models, final_time) {
  risks <- get_risk(imp, models, final_time)
  data.table(
    results = list(
      result = risks[timegroup == final_time, .(eid, risk)]
    ),
    B = unique(imp$tar_batch)
  )
}

get_sub_risk <- function(
  imp,
  from_var,
  to_var,
  duration,
  models,
  final_time,
  comp_limits
) {
  subbed <- apply_substitution(imp, from_var, to_var, duration, comp_limits)
  risks <- get_risk(subbed$results[[1]], models, final_time)
  subbed$results <- list(risks[timegroup == final_time, .(eid, risk)])
  subbed
}

average_sub_results <- function(results, df, filter_fn, result_name = "risk") {
  eids <- filter_fn(df)[, .(eid)]

  # replace the data table of results with the filtered mean
  results[,
    "results" := sapply(.SD[[1]], function(inner) {
      mean(inner[eids, on = "eid", nomatch = 0L][[result_name]], na.rm = TRUE)
    }),
    .SDcols = "results"
  ]
  results
}

merge_risks <- function(sub_risks, ref_risks, cohort) {
  sub_risks |>
    rename("sub_risk" = "results") |>
    left_join(
      ref_risks |> rename("ref_risk" = "results"),
      by = "B"
    ) |>
    mutate(
      rr = sub_risk / ref_risk
    ) |>
    dplyr::group_by(from_var, to_var, duration) |>
    dplyr::summarize(
      prop_substituted = mean(prop_substituted, na.rm = TRUE),
      mean_sub_risk = mean(sub_risk, na.rm = TRUE),
      mean_ref_risk = mean(ref_risk, na.rm = TRUE),
      mean_rr = mean(rr, na.rm = TRUE),
      lower_rr = quantile(rr, 0.025),
      upper_rr = quantile(rr, 0.975),
      cohort = cohort,
      .groups = "drop"
    )
}

adjust_colour <- function(col, factor = 1.2) {
  col_rgb <- col2rgb(col)
  col_new <- pmin(255, col_rgb * factor)
  rgb(col_new[1], col_new[2], col_new[3], max = 255)
}

make_plot <- function(df, cohort) {
  dark_factor <- switch(
    cohort,
    short_sleeper = 1.2,
    long_sleeper = 1.2,
    0.7
  )

  left_centre <- unit(0.4, "npc")
  right_centre <- unit(0.6, "npc")

  minutes_offset <- unit(-0.8, "cm")
  sub_offset <- unit(-1.2, "cm")
  sleep_offset <- unit(-1.8, "cm")
  arrow_offset <- unit(-1.5, "cm")

  left_arrow_xs <- unit(c(0.4, 0.2), "npc")
  right_arrow_xs <- unit(c(0.6, 0.8), "npc")
  arrow_ys <- unit(c(-1.5, -1.5), "cm")

  sub_types <- c("avg_mvpa", "avg_light", "avg_inactivity")
  colours <- c("#708ff9", "#6ed853", "#ff747b")
  rbindlist(Map(
    function(sub_type, colour) {
      colour <- adjust_colour(colour, dark_factor)
      sub_results <- df |>
        filter(from_var == sub_type) |>
        filter(prop_substituted > intervention_threshold) |>
        # Add in the pretend zero
        bind_rows(data.frame(
          from_var = sub_type,
          to_var = "avg_sleep",
          duration = 0,
          prop_substituted = 1,
          mean_sub_risk = 1,
          mean_ref_risk = 1,
          mean_rr = 1.0,
          lower_rr = 1.0,
          upper_rr = 1.0,
          cohort = cohort
        ))
      sub_name <- unique(case_when(
        sub_type == "avg_inactivity" ~ "Inactivity",
        sub_type == "avg_light" ~ "Light Activity",
        sub_type == "avg_mvpa" ~ "MVPA",
        .default = as.character(sub_type)
      ))
      p <- sub_results |>
        ggplot(aes(x = duration, y = mean_rr)) +
        geom_line(colour = colour) +
        geom_hline(yintercept = 1, linetype = "dotted") +
        xlab("") +
        ylab("Risk ratio") +
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
            paste("More ", sub_name),
            gp = gpar(fontsize = 12, fontfamily = "serif", fontface = 1),
            x = left_centre,
            y = sub_offset,
            hjust = 1
          )
        ) +
        annotation_custom(
          textGrob(
            paste("Less ", sub_name),
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
        scale_y_log10(
          breaks = c(0.5, 0.75, 1, 1.5, 2),
          labels = scales::label_number()
        ) +
        coord_cartesian(
          xlim = c(-60, 60),
          ylim = c(0.5, 2),
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
          axis.line = element_line(color = "#585656")
        ) +
        geom_ribbon(
          aes(ymin = lower_rr, ymax = upper_rr),
          alpha = 0.25,
          fill = colour
        )
      data.table(
        sub_type = sub_type,
        cohort = cohort,
        plot = list(p)
      )
    },
    sub_types,
    colours
  ))
}
