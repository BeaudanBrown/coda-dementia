# Basic Linear Regression and isotemporal substitution with Compositional Explanatory Variables

# Packages and installs
library(tidyverse)
library(compositions)
library(data.table)
library(here)
library(dotenv)
library(rms)
library(survival)
library(ggplot2)
library(mvtnorm)

sleep_idx <- 1
inactive_idx <- 2
light_idx <- 3
vig_idx <- 4
sub_amount <- 120
mins_in_day <- 1440
hrs_in_day <- 24
short_sleep_hours <- 6
long_sleep_hours <- 9

source("dem_models.R")
source("mri_models.R")
source("ideal_comp.R")

inc <- -sub_amount:sub_amount / mins_in_day

mri_avg_sleep_comp <-
  mri_comp[mri_comp$mri_df.sleep_n >= short_sleep_hours / hrs_in_day &
           mri_comp$mri_df.sleep_n <= long_sleep_hours / hrs_in_day]

avg_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n >= short_sleep_hours / hrs_in_day &
           dem_comp$dem_df.sleep_n <= long_sleep_hours / hrs_in_day]
short_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n < short_sleep_hours / hrs_in_day]
long_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n > long_sleep_hours / hrs_in_day]

# generate_plot <- function(conditions, base_offset, mean_ilr) {
generate_plot <- function(title, substitution, cohorts, model) {
  # Create an empty data frame
  df <- data.frame()
  scaled_x <- inc * mins_in_day

  # Process each condition
  for (cohort in cohorts) {
    # Create offset matrix
    adjusted_comp <- cohort$base_comp
    mean_ilr <- cohort$ilr

    adjusted_comp[, substitution[1]] <- adjusted_comp[, substitution[1]] + inc

    other_idx <- setdiff(1:ncol(cohort$base_comp), substitution[1])

    for (i in other_idx) {
      adjusted_comp[, i] <-
        adjusted_comp[, i] - (inc * (adjusted_comp[, i] / rowSums(adjusted_comp[, -substitution[1]])))
    }

    # Calculate ilrs
    ilrs <- t(apply(adjusted_comp, 1, function(comp) ilr(acomp(comp), V = v)))

    # Get contrasts
    contrasts <- apply(ilrs, 1, function(new_ilr) {
      contrast_out <- contrast(
        model,
        list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
        list(R1 = mean_ilr[1], R2 = mean_ilr[2], R3 = mean_ilr[3])
      )
      return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
    })

    df <- rbind(df, data.frame(
      X = scaled_x,
      Y = contrasts[1, ],
      Y_lower = contrasts[2, ],
      Y_upper = contrasts[3, ],
      name = cohort$name
    ))
  }

  return(ggplot(df, aes(x = X)) +
    geom_line(aes(y = Y, color = name)) +
    geom_ribbon(aes(ymin = Y_lower, ymax = Y_upper, fill = name), alpha = 0.2) +
    facet_wrap(~name) +
    xlab("min/day realocated") +
    ylab("Delta log hazard") +
    geom_hline(yintercept = 0) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = "darkblue")))
}

build_cohort <- function(name, comp, v) {
  geo_mean <-
    apply(comp, 2, function(x) exp(mean(log(x))))
  avg_ilr <-
    ilr(acomp(geo_mean), V = v)
  base_comp <-
    matrix(rep(geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)

  return(list(name = name, base_comp = base_comp, ilr = avg_ilr))
}

cohorts <- list(
  build_cohort("1 Short Sleepers", short_sleep_comp, v),
  build_cohort("2 Average Sleepers", avg_sleep_comp, v),
  build_cohort("3 Long Sleepers", long_sleep_comp, v)
)

mri_cohorts <- list(
  list(
    name = "Average Sleepers",
    base_comp = mri_avg_offset,
    ilr = mri_avg_ilr
  )
)

sleep_inactive_plot <-
  generate_plot("Substituting sleep for inactivity", c(sleep_idx, inactive_idx), cohorts, dem_model)
sleep_light_plot <-
  generate_plot("Substituting sleep for light", c(sleep_idx, light_idx), cohorts, dem_model)
sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), cohorts, dem_model)

ggsave(filename = "plots/avg_sub_sleep_plot.png", plot = sub_sleep_avg_plot)
ggsave(filename = "plots/short_sleep_sub_sleep_plot.png", plot = short_sleep_sub_sleep_plot)
ggsave(filename = "plots/long_sleep_sub_sleep_plot.png", plot = long_sleep_sub_sleep_plot)
