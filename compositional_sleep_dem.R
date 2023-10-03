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

source("dem_models.R")
source("mri_models.R")
source("ideal_comp.R")

sleep_idx <- 1
inactive_idx <- 2
light_idx <- 3
vig_idx <- 4
sub_amount <- 120
mins_in_day <- 1440
hrs_in_day <- 24
short_sleep_hours <- 6
long_sleep_hours <- 9

mri_avg_sleep_comp <-
  mri_comp[mri_comp$mri_df.sleep_n >= short_sleep_hours / hrs_in_day & mri_comp$mri_df.sleep_n <= long_sleep_hours / hrs_in_day]

avg_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n >= short_sleep_hours / hrs_in_day & dem_comp$dem_df.sleep_n <= long_sleep_hours / hrs_in_day]
short_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n < short_sleep_hours / hrs_in_day]
long_sleep_comp <-
  dem_comp[dem_comp$dem_df.sleep_n > long_sleep_hours / hrs_in_day]

geo_mean <-
  apply(dem_comp, 2, function(x) exp(mean(log(x))))

avg_sleep_geo_mean <-
  apply(avg_sleep_comp, 2, function(x) exp(mean(log(x))))
mri_avg_sleep_geo_mean <-
  apply(mri_avg_sleep_comp, 2, function(x) exp(mean(log(x))))
short_sleep_geo_mean <-
  apply(short_sleep_comp, 2, function(x) exp(mean(log(x))))
long_sleep_geo_mean <-
  apply(long_sleep_comp, 2, function(x) exp(mean(log(x))))

avg_ilr <-
  ilr(acomp(avg_sleep_geo_mean), V = v)
mri_avg_ilr <-
  ilr(acomp(mri_avg_sleep_geo_mean), V = v)
short_sleep_mean_ilr <-
  ilr(acomp(short_sleep_geo_mean), V = v)
long_sleep_mean_ilr <-
  ilr(acomp(long_sleep_geo_mean), V = v)

inc <- -sub_amount:sub_amount / mins_in_day
avg_offset <-
  matrix(rep(avg_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
mri_avg_offset <-
  matrix(rep(mri_avg_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
short_sleep_base_offset <-
  matrix(rep(short_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)
long_sleep_base_offset <-
  matrix(rep(long_sleep_geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)

# List of replacement prefixes and their respective pairs of indexes
# These are then used to generate and plot the relevant contrasts

sub_sleep_conditions <- list(
  inactive_sleep = c(inactive_idx, sleep_idx),
  light_sleep = c(light_idx, sleep_idx),
  vig_sleep = c(vig_idx, sleep_idx)
)

sub_inactive_conditions <- list(
  sleep_inactive = c(sleep_idx, inactive_idx),
  light_inactive = c(light_idx, inactive_idx),
  vig_inactive = c(vig_idx, inactive_idx)
)

sub_light_conditions <- list(
  inactive_light = c(inactive_idx, light_idx),
  sleep_light = c(sleep_idx, light_idx),
  vig_light = c(vig_idx, light_idx)
)

sub_vig_conditions <- list(
  inactive_vig = c(inactive_idx, vig_idx),
  light_vig = c(light_idx, vig_idx),
  sleep_vig = c(sleep_idx, vig_idx)
)

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
    # adjusted_comp[, substitution[2]] <- adjusted_comp[, substitution[2]] - inc

    # filtered_indices <- any(adjusted_comp <= 0)
    # adjusted_comp <- adjusted_comp[filtered_indices,]
    # scaled_x <- scaled_x[filtered_indices] * mins_in_day

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

sub_sleep_flat_plot <-
  flat_sub(avg_offset, avg_ilr)
short_sleep_flat_plot <-
  flat_sub(short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_flat_plot <-
  flat_sub(long_sleep_base_offset, long_sleep_mean_ilr)

cohorts <- list(
  list(
    name = "1 Short Sleepers",
    base_comp = short_sleep_base_offset,
    ilr = short_sleep_mean_ilr
  ),
  list(
    name = "2 Average Sleepers",
    base_comp = avg_offset,
    ilr = avg_ilr
  ),
  list(
    name = "3 Long Sleepers",
    base_comp = long_sleep_base_offset,
    ilr = long_sleep_mean_ilr
  )
)

mri_cohorts <- list(
  list(
    name = "Average Sleepers",
    base_comp = mri_avg_offset,
    ilr = mri_avg_ilr
  )
)

sleep_inactive_plot <-
  generate_plot("Substituting sleep for inactivity", c(sleep_idx, 0), cohorts)
sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), cohorts)
sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), cohorts)

mri_sleep_mvpa_plot <-
  generate_plot("Substituting sleep for MVPA", c(sleep_idx, vig_idx), mri_cohorts, tbv_model)
mri_sleep_inactive_plot <-
  generate_plot("Substituting sleep for inactivity", c(sleep_idx, inactive_idx), mri_cohorts, tbv_model)
mri_sleep_light_plot <-
  generate_plot("Substituting sleep for light", c(sleep_idx, light_idx), mri_cohorts, tbv_model)

sub_sleep_avg_plot <-
  generate_plot(sub_sleep_conditions, avg_offset, avg_ilr)
short_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_sleep_plot <-
  generate_plot(sub_sleep_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_inactive_avg_plot <-
  generate_plot(sub_inactive_conditions, avg_offset, avg_ilr)
short_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_inactive_plot <-
  generate_plot(sub_inactive_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_light_avg_plot <-
  generate_plot(sub_light_conditions, avg_offset, avg_ilr)
short_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_light_plot <-
  generate_plot(sub_light_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

sub_vig_avg_plot <-
  generate_plot(sub_vig_conditions, avg_offset, avg_ilr)
short_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, short_sleep_base_offset, short_sleep_mean_ilr)
long_sleep_sub_vig_plot <-
  generate_plot(sub_vig_conditions, long_sleep_base_offset, long_sleep_mean_ilr)

ggsave(filename = "plots/flat_avg_plot.png", plot = sub_sleep_flat_plot, bg = "white")
ggsave(filename = "plots/flat_short_sleep_plot.png", plot = short_sleep_flat_plot, bg = "white")
ggsave(filename = "plots/flat_long_sleep_plot.png", plot = long_sleep_flat_plot, bg = "white")

ggsave(filename = "plots/avg_sub_sleep_plot.png", plot = sub_sleep_avg_plot)
ggsave(filename = "plots/short_sleep_sub_sleep_plot.png", plot = short_sleep_sub_sleep_plot)
ggsave(filename = "plots/long_sleep_sub_sleep_plot.png", plot = long_sleep_sub_sleep_plot)

ggsave(filename = "plots/avg_sub_inactive_plot.png", plot = sub_inactive_avg_plot)
ggsave(filename = "plots/short_sleep_sub_inactive_plot.png", plot = short_sleep_sub_inactive_plot)
ggsave(filename = "plots/long_sleep_sub_inactive_plot.png", plot = long_sleep_sub_inactive_plot)

ggsave(filename = "plots/avg_sub_light_plot.png", plot = sub_light_avg_plot)
ggsave(filename = "plots/short_sleep_sub_light_plot.png", plot = short_sleep_sub_light_plot)
ggsave(filename = "plots/long_sleep_sub_light_plot.png", plot = long_sleep_sub_light_plot)

ggsave(filename = "plots/avg_sub_vig_plot.png", plot = sub_vig_avg_plot)
ggsave(filename = "plots/short_sleep_sub_vig_plot.png", plot = short_sleep_sub_vig_plot)
ggsave(filename = "plots/long_sleep_sub_vig_plot.png", plot = long_sleep_sub_vig_plot)

