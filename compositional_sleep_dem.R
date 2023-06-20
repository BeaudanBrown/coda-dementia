# Basic Linear Regression and isotemporal substitution with Compositional Explanatory Variables

# Packages and installs
library(tidyverse)
library(compositions)
library(data.table)
library(here)
library(dotenv)
library(rms)
library(survival)

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
full_df$in_n <-
  ((full_df$dur_day_total_IN_min_pla + full_df$dur_spt_wake_IN_min_pla) / full_df$mins_worn) * mins_in_day
full_df$lig_n <-
  ((full_df$dur_day_total_LIG_min_pla + full_df$dur_spt_wake_LIG_min_pla) / full_df$mins_worn) * mins_in_day
full_df$mod_n <-
  ((full_df$dur_day_total_MOD_min_pla + full_df$dur_spt_wake_MOD_min_pla) / full_df$mins_worn) * mins_in_day
full_df$vig_n <-
  ((full_df$dur_day_total_VIG_min_pla + full_df$dur_spt_wake_VIG_min_pla) / full_df$mins_worn) * mins_in_day

# TODO: Figure out tf Lach mean by this
# ill add awake-sleep and sleep just for a neat 4 part composition
# we need to discuss actually how we will handle this though
# could also add these components to the day components
# or have a 5 part composition that describes the activity done in the sleep period
# to discuss!

full_df$mvpa_n <- full_df$mod_n + full_df$vig_n

contrast_substitution <- function(df, sbp, sub_vec) {
  comp <- acomp(data.frame(df$sleep_n, df$in_n, df$lig_n, df$mvpa_n))
  geo_mean <-
    apply(comp, 2, function(x) exp(mean(log(x))))

  v <- gsi.buildilrBase(t(sbp))

  base_ilr <-
    ilr(comp, V = v) |>
    setNames(c("R1", "R2", "R3"))
  mean_ilr <-
    ilr(acomp(geo_mean), V = v)
  new_ilr <-
    acomp(geo_mean + sub_vec) |>
    ilr(V = v)

  model_data <- select(df, dem, time_to_dem, sex, age_assessment) |> cbind(base_ilr)

  dd <- datadist(model_data)
  options(datadist = "dd")
  model <-
    cph(
        Surv(time_to_dem, dem) ~
          rcs(R1, 3) + rcs(R2, 3) +
          rcs(R3, 3) +
          sex +
          rcs(age_assessment, 3),
        data = model_data
    )

  contrast <- contrast(
    model,
    list(R1 = mean_ilr[1], R2 = mean_ilr[2], R3 = mean_ilr[3]),
    list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3])
  )
  return(contrast$Contrast)
}
