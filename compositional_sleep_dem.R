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
library(cowplot)

sleep_idx <- 1
inactive_idx <- 2
light_idx <- 3
vig_idx <- 4
mins_in_day <- 1440
hours_in_day <- 24

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

# Variables and source data
source_df <-
  fread(file.path(data_dir, "Accelerometery/Processed_GGIR/part5_personsumMM_output.csv"))

# Merge in outcome data
outcome_df <-
  fread(file.path(data_dir, "SRI/sri_data_may2023.csv"), stringsAsFactors = TRUE) |>
  as_tibble() |>
  select(eid, sex, BMI, age_assessment, time_to_dem, date_of_death, accel_date, date_acdem2, dem)

# Death from other causes censoring event
full_df <- left_join(outcome_df, source_df, by = "eid") |>
  mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |>
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, accel_date),
    competing == 1 ~ difftime(date_of_death, accel_date),
    TRUE ~ difftime("2022-01-01", accel_date)
  )) |>
  mutate(time_to_dem = as.integer(time_to_dem))

# TODO: Consult data dictionary to figure out why totals we are calculating don't match existing totals
full_df$awake_sleep <-
  full_df$dur_spt_wake_IN_min_pla +
  full_df$dur_spt_wake_LIG_min_pla +
  full_df$dur_spt_wake_MOD_min_pla +
  full_df$dur_spt_wake_VIG_min_pla

full_df$mins_worn <-
  full_df$dur_spt_sleep_min_pla +
  full_df$dur_day_total_IN_min_pla +
  full_df$dur_day_total_LIG_min_pla +
  full_df$dur_day_total_MOD_min_pla +
  full_df$dur_day_total_VIG_min_pla +
  full_df$awake_sleep

# Variables normalised to 1440 relative to their proportion of total wear time
full_df$sleep_n <-
  (full_df$dur_spt_sleep_min_pla / full_df$mins_worn) * mins_in_day
full_df$inactive_n <-
  ((full_df$dur_day_total_IN_min_pla + full_df$dur_spt_wake_IN_min_pla) / full_df$mins_worn) * mins_in_day
full_df$light_n <-
  ((full_df$dur_day_total_LIG_min_pla + full_df$dur_spt_wake_LIG_min_pla) / full_df$mins_worn) * mins_in_day
full_df$moderate_n <-
  ((full_df$dur_day_total_MOD_min_pla + full_df$dur_spt_wake_MOD_min_pla) / full_df$mins_worn) * mins_in_day
full_df$vigorous_n <-
  ((full_df$dur_day_total_VIG_min_pla + full_df$dur_spt_wake_VIG_min_pla) / full_df$mins_worn) * mins_in_day

full_df$mvpa_n <- full_df$moderate_n + full_df$vigorous_n

comp <- acomp(data.frame(full_df$sleep_n, full_df$inactive_n, full_df$light_n, full_df$mvpa_n))
geo_mean <-
  apply(comp, 2, function(x) exp(mean(log(x))))

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)
v <- gsi.buildilrBase(t(sbp))

base_ilr <-
  ilr(comp, V = v) |>
  setNames(c("R1", "R2", "R3"))
mean_ilr <-
  ilr(acomp(geo_mean), V = v)

model_data <- select(full_df, dem, time_to_dem, sex, age_assessment) |> cbind(base_ilr)

options(datadist = datadist(model_data))
model <-
  cph(
      Surv(time_to_dem, dem) ~
        rcs(R1, 3) + rcs(R2, 3) +
        rcs(R3, 3) +
        sex +
        rcs(age_assessment, 3),
      data = model_data
  )

inc <- -90:90 / 1440
scaled_x <- inc * 1440
base_offset <- matrix(rep(geo_mean, times = length(inc)), nrow = length(inc), byrow = TRUE)

# List of replacement prefixes and their respective pairs of indexes
# These are then used to generate and plot the relevant contrasts
conditions <- list(sleep_inactive = c(sleep_idx, inactive_idx),
                   sleep_light = c(sleep_idx, light_idx),
                   sleep_vig = c(sleep_idx, vig_idx))

# Initialize your contrast list
contrasts <- list()

# Process each condition
for (cond in names(conditions)) {
  # Create offset matrix
  offset <- base_offset
  offset[, conditions[[cond]][1]] <- base_offset[, conditions[[cond]][1]] + inc
  offset[, conditions[[cond]][2]] <- base_offset[, conditions[[cond]][2]] - inc

  # Calculate ilrs
  ilrs <- t(apply(offset, 1, function(comp) ilr(acomp(comp), V = v)))

  # Get contrasts
  contrasts[[cond]] <- apply(ilrs, 1, function(new_ilr) {
    contrast_out <- contrast(
      model,
      list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
      list(R1 = mean_ilr[1], R2 = mean_ilr[2], R3 = mean_ilr[3])
    )
    return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
  })
}

# Create an empty data frame
df <- data.frame()

# Iterate over conditions and add data to the dataframe
for (cond in names(conditions)) {
  temp_df <- data.frame(
    X = scaled_x,
    Y = contrasts[[cond]][1, ],
    Y_lower = contrasts[[cond]][2, ],
    Y_upper = contrasts[[cond]][3, ],
    condition = cond
  )
  df <- rbind(df, temp_df)
}

# Plot using ggplot
ggplot(df, aes(x = X)) +
  geom_line(aes(y = Y, color = condition)) +
  geom_ribbon(aes(ymin = Y_lower, ymax = Y_upper, fill = condition), alpha = 0.2) +
  xlab("min/day realocated") +
  ylab("Delta log hazard") +
  theme_cowplot()
