source("utils.R")
source("dem_models.R")

# load packages
list_of_packages <- c(
  "mice",
  "tidyverse",
  "survival",
  "rms",
  "compositions",
  "data.table",
  "parallel",
  "boot",
  "renv",
  "rlang",
  "cowplot",
  "extrafont"
)

lapply(list_of_packages, library, character.only = TRUE)

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data_26_04_24.rds"))

## Single imputation

predmat <- quickpred(boot_data,
  mincor = 0,
  exclude = c(
    "avg_sleep", "avg_inactivity", "avg_light",
    "avg_mvpa", "eid"
  )
)

imp <- mice(
  boot_data,
  m = 1,
  predictorMatrix = predmat, maxit = maxit
)

imp <- complete(imp)
setDT(imp)
imp[, id := .I]
imp_len <- nrow(imp)

## reference compositions

all_comp <- acomp(
  imp[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")]
)

short_sleep_comp <- all_comp[
  all_comp$avg_sleep < short_sleep_hours / hrs_in_day,
]
short_sleep_geo_mean <- acomp(apply(
  short_sleep_comp, 2, function(x) exp(mean(log(x)))
))

avg_sleep_comp <- all_comp[
  all_comp$avg_sleep >= short_sleep_hours / hrs_in_day,
]
avg_sleep_geo_mean <- acomp(
  apply(avg_sleep_comp, 2, function(x) exp(mean(log(x))))
)

## Remove first three years of follow-up

imp2 <- imp[imp$time_to_dem > (3 * 365), ]
imp_len <- nrow(imp2)

## constants
# data for g-computation/standardisation

min_age_of_dem <- min(imp2$age_dem)
max_age_of_dem <- max(imp2$age_dem)
age_range <- max_age_of_dem - min_age_of_dem
timegroup_steps <- ceiling(age_range * 2)
median_age_of_dem <- median(imp2[imp2$dem == 1, ]$age_dem)

timegroup_cuts <-
  seq(
    from = min_age_of_dem,
    to = max_age_of_dem,
    length.out = timegroup_steps
  )

median_age_of_dem_timegroup <-
  which(timegroup_cuts > median_age_of_dem)[1] - 1

## Fit models

models <- fit_model(imp2, timegroup_cuts, get_primary_formula)

# data for g-computation/standardisation
imp2 <- imp2[rep(seq_len(imp_len), each = median_age_of_dem_timegroup)]
imp2[, timegroup := rep(1:median_age_of_dem_timegroup, imp_len)]

## Estimate substitution effects

short_sleep_inactive <-
  calc_substitution(short_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_inactivity"),
    timegroup = timegroup
  )

short_sleep_light <-
  calc_substitution(short_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_light"),
    timegroup = timegroup
  )

short_sleep_mvpa <-
  calc_substitution(short_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_mvpa"),
    timegroup = timegroup
  )

avg_sleep_inactive <-
  calc_substitution(avg_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_inactivity"),
    timegroup = timegroup
  )

avg_sleep_light <-
  calc_substitution(avg_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_light"),
    timegroup = timegroup
  )

avg_sleep_mvpa <-
  calc_substitution(avg_sleep_geo_mean,
    imp2,
    models[["model_dem"]],
    models[["model_death"]],
    models[["model_formula"]],
    c("avg_sleep", "avg_mvpa"),
    timegroup = timegroup
  )

full_df <-
  full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
  full_join(short_sleep_mvpa, by = "offset") |>
  full_join(avg_sleep_inactive, by = "offset") |>
  full_join(avg_sleep_light, by = "offset") |>
  full_join(avg_sleep_mvpa, by = "offset")

## Plot

plot_data <-
  pivot_longer(full_df,
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
  mutate(Reference = ifelse(
    str_detect(Substitution, "short_sleep"),
    "Short sleepers",
    "Normal sleepers"
  )) |>
  mutate(Substitution = str_remove(Substitution, "_short_sleep_geo_mean")) |>
  mutate(Substitution = str_remove(Substitution, "_avg_sleep_geo_mean")) |>
  mutate(Substitution = ifelse(
    Substitution == "avg_inactivity",
    "Inactivity",
    ifelse(Substitution == "avg_light", "Light activity", "MVPA")
  ))

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

  plot_data2 <- plot_data[
    plot_data$Substitution == sub & plot_data$Reference == refcomp,
  ]


  plot_data2 |>
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
    coord_cartesian(
      ylim = c(0.33, 3),
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
    )
}

# normal sleepers
p1 <- rr_plot("inactivity", "Normal sleepers", "#fc020f")
p2 <- rr_plot("light activity", "Normal sleepers", "#145e01")
p3 <- rr_plot("MVPA", "Normal sleepers", "#011869")

pnorm <-
  plot_grid(
    NULL,
    p1 + labs(x = "", title = "A") + theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 16)
    ),
    p2 + labs(x = "", title = "B") + theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 16)
    ),
    p3 + labs(x = "", title = "C") + theme(
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
    p4 + labs(x = "", title = "D") + theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 16)
    ),
    p5 + labs(x = "", title = "E") + theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 16)
    ),
    p6 + labs(x = "", title = "F") + theme(
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

plot <- plot_grid(pnorm,
  pshort,
  nrow = 1
)

ggsave(
  file.path(
    data_dir,
    "../Manuscript/Appendix_figures/truncate_3years.svg"
  ),
  plot = plot,
  device = "svg",
  width = 10,
  height = 12
)
