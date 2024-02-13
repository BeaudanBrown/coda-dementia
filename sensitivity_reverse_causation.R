#### censoring competing events at death ####
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

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

if (length(new_packages)) install.packages(new_packages)

lapply(list_of_packages, library, character.only = TRUE)

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

## Construct risk set 

# death from non-dementia causes is a censoring event

boot_data <- boot_data |>
  mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |>
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, date_accel),
    competing == 1 ~ difftime(date_of_death, date_accel),
    TRUE ~ difftime("2022-01-01", date_accel)
  )) |>
  mutate(time_to_dem = as.integer(time_to_dem))

## Single imputation

# set date variables to strings to avoid errors
boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)

# run imputation (or read from disc if already complete)

if (file.exists(file.path(data_dir,"rc_imp.rds"))) {
  imp <- read_rds(file.path(data_dir,"rc_imp.rds"))
} else {
  predmat <- quickpred(boot_data,
                       mincor = 0,
                       exclude = c(
                         "date_acdem2", "date_accel", "date_of_death",
                         "avg_sleep", "avg_inactivity", "avg_light",
                         "avg_mvpa", "eid"
                       )
  )
  
  predmat["date_acdem2", ] <- 0
  predmat["date_of_death", ] <- 0
  
  # method for each imputed variable
  imp_methods <- make.method(boot_data)
  # exclude dates from being imputed
  imp_methods["date_acdem2"] <- ""
  imp_methods["date_of_death"] <- ""
  
  imp <-
    mice(boot_data,
         m = 1,
         predictorMatrix = predmat,
         methods = imp_methods)
  
  imp <- complete(imp)
  
  write_rds(imp, file.path(data_dir,"rc_imp.rds"))
}


## reference compositions

all_comp <- acomp(imp[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

short_sleep_comp <-
  all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
short_sleep_geo_mean <-
  acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

avg_sleep_comp <-
  all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
avg_sleep_geo_mean <-
  acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))


## fit model

# Person period format dataset

imp_long <- survSplit(Surv(time = time_to_dem, event = dem) ~ .,
                      data = imp,
                      cut = seq(
                        from = 0,
                        to = max(imp$time_to_dem),
                        length.out = 10
                      ),
                      episode = "timegroup", end = "dem_time", event = "dem"
)

## fit model 

knots_timegroup <- quantile(imp_long[["timegroup"]], c(0.1, 0.5, 0.9))
knots_deprivation <- quantile(imp_long[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))

# compare models with/without timegroup interaction

mod <- lrm(dem ~ 
             rcs(timegroup, knots_timegroup)*pol(R1,2) +
             rcs(timegroup, knots_timegroup)*pol(R2,2) +
             rcs(timegroup, knots_timegroup)*pol(R3,2) +
             sex +
             retired +
             shift +
             apoe_e4 +
             highest_qual +
             rcs(townsend_deprivation_index, knots_deprivation) +
             antidepressant_med +
             antipsychotic_med +
             insomnia_med +
             ethnicity +
             avg_total_household_income +
             smok_status,
           data = imp_long)

mod2 <- lrm(dem ~ 
              rcs(timegroup, knots_timegroup) +
              pol(R1,2) +
              pol(R2,2) +
              pol(R3,2) +
              sex +
              retired +
              shift +
              apoe_e4 +
              highest_qual +
              rcs(townsend_deprivation_index, knots_deprivation) +
              antidepressant_med +
              antipsychotic_med +
              insomnia_med +
              ethnicity +
              avg_total_household_income +
              smok_status,
            data = imp_long)

lrtest(mod, mod2)


# Model for risk ratios

mod <- glm(dem ~ 
             rcs(timegroup, knots_timegroup)*poly(R1,2) +
             rcs(timegroup, knots_timegroup)*poly(R2,2) +
             rcs(timegroup, knots_timegroup)*poly(R3,2) +
             sex +
             retired +
             shift +
             apoe_e4 +
             highest_qual +
             rcs(townsend_deprivation_index, knots_deprivation) +
             antidepressant_med +
             antipsychotic_med +
             insomnia_med +
             ethnicity +
             avg_total_household_income +
             smok_status,
           data = imp_long,
           family = "binomial")


## data for g-computation/standardisation

timegroup <- 3

setDT(imp)
imp[, id := .I]
imp_len <- nrow(imp)
imp_long <- imp[rep(seq_len(imp_len), each = timegroup)]
imp_long[, timegroup := rep(1:timegroup, imp_len)]


## Estimate substitution effects

short_sleep_inactive <-
  calc_substitution(short_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_inactivity"),
                    timegroup = timegroup)

short_sleep_light <-
  calc_substitution(short_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_light"),
                    timegroup = timegroup)

short_sleep_mvpa <-
  calc_substitution(short_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_mvpa"),
                    timegroup = timegroup)

avg_sleep_inactive <-
  calc_substitution(avg_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_inactivity"),
                    timegroup = timegroup)

avg_sleep_light <-
  calc_substitution(avg_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_light"),
                    timegroup = timegroup)

avg_sleep_mvpa <-
  calc_substitution(avg_sleep_geo_mean,
                    imp_long,
                    mod,
                    c("avg_sleep", "avg_mvpa"),
                    timegroup = timegroup)

full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
  full_join(short_sleep_mvpa, by = "offset") |>
  full_join(avg_sleep_inactive, by = "offset") |>
  full_join(avg_sleep_light, by = "offset") |>
  full_join(avg_sleep_mvpa, by = "offset")


## Plot 

plot_data <-
  pivot_longer(full_df,
               -offset,
               values_to = "risk",
               names_to = "Substitution") |>
  group_by(Substitution) |>
  mutate(ref_risk = ifelse(offset == 0, risk, NA_real_)) |>
  fill(ref_risk, .direction = "downup") |>
  mutate(risk_dif = risk - ref_risk,
         risk_ratio = risk / ref_risk) |>
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
  
  plot_data2 <- plot_data[plot_data$Substitution == sub &
                            plot_data$Reference == refcomp,]
  
  
  plot_data2 |>
    ggplot(aes(x = offset, y = risk_ratio)) +
    geom_line(colour = colour) +
    facet_wrap( ~ Substitution, nrow = 2) +
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
    geom_segment(aes(
      x = 1,
      y = -0.625,
      xend = 15,
      yend = -0.625
    ),
    arrow = arrow(length = unit(0.15, "cm"))) +
    geom_segment(aes(
      x = -1,
      y = -0.625,
      xend = -15,
      yend = -0.625
    ),
    arrow = arrow(length = unit(0.15, "cm"))) +
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
    coord_cartesian(ylim = c(0.33, 3),
                    expand = FALSE,
                    clip = "off") +
    cowplot::theme_cowplot() +
    theme(
      text = element_text(size = 12, family = "serif"),
      plot.margin = unit(c(1, 1, 4, 1), "lines"),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
}
  
# normal sleepers
p1 <- rr_plot("inactivity", "Normal sleepers", "#fc020f")
p2 <- rr_plot("light activity", "Normal sleepers", "#145e01")
p3 <- rr_plot("MVPA", "Normal sleepers", "#011869")

pnorm <-
  plot_grid(
    NULL,
    p1 + labs(x = "", title = "A") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    p2 + labs(x = "", title = "C") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    p3 + labs(x = "", title = "E") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    align = "vh",
    rel_heights = c(0.05,1,1,1),
    nrow = 4,
    labels = "Normal sleepers",
    hjust = -1
  )

# short sleepers
p4 <- rr_plot("inactivity", "Short sleepers", "#ff747b")
p5 <- rr_plot("light activity", "Short sleepers", "#6ed853")
p6 <- rr_plot("MVPA", "Short sleepers", "#708ff9")

pshort <-
  plot_grid(
    NULL,
    p4 + labs(x = "", title = "A") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    p5 + labs(x = "", title = "C") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    p6 + labs(x = "", title = "E") + theme(legend.position = "none",plot.title.position = "plot",plot.title = element_text(size=16)),
    align = "vh",
    rel_heights = c(0.05,1,1,1),
    nrow = 4,
    labels = "Normal sleepers",
    hjust = -1
  )

plot <- plot_grid(pnorm,
                  pshort,
                  nrow = 1)

ggsave(
  file.path(
    data_dir,
    "../../Papers/Substitution Analysis/Appendix_figures/first_3years.png"
  ),
  plot = plot,
  device = "png",
  bg = "white",
  width = 10,
  height = 12,
  dpi = 500
)
