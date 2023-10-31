# Helper functions for fitting model and calculating risks
library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)
library(survival)
library(compositions)
library(mice)

# constants
mins_in_day <- 1440
sub_steps <- 4
sub_step_mins <- 15

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

## Define SBP

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)

v <- gsi.buildilrBase(t(sbp))

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

### Imputation ###
# set date variables to strings to avoid errors

boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)

# Matrix of variables to include in imputation model

predmat <- quickpred(boot_data,
  mincor = 0,
  exclude = c(
    "date_acdem2", "date_accel", "date_of_death",
    "avg_sleep", "avg_inactivity", "avg_light",
    "avg_mvpa"
  )
)
predmat["date_acdem2", ] <- 0
predmat["date_of_death", ] <- 0

# method for each imputed variable
imp_methods <- make.method(boot_data)
# exclude dates from being imputed
imp_methods["date_acdem2"] <- ""
imp_methods["date_of_death"] <- ""

fit_model <- function(imp, reg_formula) {
  imp_long <- survSplit(Surv(time = age_accel, event = dem, time2 = age_dem) ~ .,
    data = imp,
    cut = seq(
      from = min(imp$age_dem),
      to = max(imp$age_dem),
      length.out = 76
    ),
    episode = "timegroup", end = "age_end", event = "dem",
    start = "age_start"
  )

  model <- glm(reg_formula, data = imp_long, family = binomial)

  return(model)
}

calc_risk <- function(composition, stacked_data, model, timegroup) {
  ilr <- ilr(composition, V = v)

  ilr_data <-
    mutate(stacked_data,
           R1 = ilr[1], R2 = ilr[2], R3 = ilr[3])

  ilr_data$haz <-
    predict(model, newdata = ilr_data, type = "response")

  ilr_data <- ilr_data |>
    group_by(id) |>
    arrange(timegroup) |>
    mutate(risk = 1 - cumprod(1 - haz)) |>
    ungroup()

  risk <- ilr_data |>
    filter(timegroup == timegroup) |>
    summarise(mean = mean(risk))

  return(risk)
}

calc_substitution <- function(base_comp, imp_stacked, model, substitution, timegroup) {
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  sub_comps <- data.frame(matrix(rep(base_comp, length(inc)), nrow = length(inc), byrow = TRUE))
  colnames(sub_comps) <- c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
  sub_comps[, substitution[1]] <- sub_comps[, substitution[1]] + inc
  sub_comps[, substitution[2]] <- sub_comps[, substitution[2]] - inc
  sub_risks <- bind_rows(apply(sub_comps, 1, function(comp) calc_risk(acomp(comp), imp_stacked, model, timegroup)))

  return(data.frame(
    offset = inc * mins_in_day,
    risks = sub_risks
  ) |>
    setNames(c("offset",
               paste(substitution[2], deparse(substitute(base_comp)), sep = "_"))))
}

## Return results from substitutions

## reference compositions

all_comp <- acomp(boot_data[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

short_sleep_comp <-
  all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
short_sleep_geo_mean <-
  acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

avg_sleep_comp <-
  all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
avg_sleep_geo_mean <-
  acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

sub_results <- function(model, imp, timegroup) {
  # set up stacked dataset for g-computation (one row for each follow-up interval)
  imp_stacked <- do.call("rbind", replicate(timegroup, imp, simplify = FALSE))

  imp_stacked$timegroup <- rep(1:timegroup, nrow(imp))

  ## substitutions

  short_sleep_inactive <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  short_sleep_light <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked, model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  short_sleep_mvpa <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_mvpa"),
                      timegroup = timegroup)

  avg_sleep_inactive <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  avg_sleep_light <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  avg_sleep_mvpa <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_mvpa"),
                      timegroup = timegroup)

  full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset") |>
    pivot_longer(-offset, values_to = "risk", names_to = "Substitution") |>
    group_by(Substitution) |>
    mutate(ref_risk = ifelse(offset == 0, risk, NA_real_)) |>
    fill(ref_risk, .direction = "downup") |>
    mutate(risk_dif = risk - ref_risk,
           risk_ratio = risk / ref_risk) |>
    mutate(Reference = ifelse(str_detect(Substitution, "short_sleep"),
                              "Short sleepers", "Normal sleepers")) |>
    mutate(Substitution = str_remove(Substitution, "_short_sleep_geo_mean")) |>
    mutate(Substitution = str_remove(Substitution, "_avg_sleep_geo_mean")) |>
    mutate(Substitution = ifelse(Substitution == "avg_inactivity", "Inactivity",
                                 ifelse(Substitution == "avg_light", "Light activity", "MVPA")))

  return(full_df)
}


## bootstrap function

boot_fn <- function(data, indices, reg_formula, timegroup) {
  this_sample <- data[indices, ]

  # imputation
  imp <- mice(this_sample, m = 1,
              maxit = 1,
              predictorMatrix = predmat,
              methods = imp_methods)

  imp <- complete(imp)
  imp$id <- seq_len(nrow(imp))

  # fit model
  model <- fit_model(imp, reg_formula)

  # data for g-computation/standardisation
  imp_stacked <- do.call("rbind", replicate(timegroup, imp, simplify = FALSE))
  imp_stacked$timegroup <- rep(1:timegroup, nrow(imp))

  # substitutions for each reference comp
  short_sleep_inactive <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  short_sleep_light <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  short_sleep_mvpa <-
    calc_substitution(short_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_mvpa"),
                      timegroup = timegroup)

  avg_sleep_inactive <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  avg_sleep_light <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  avg_sleep_mvpa <-
    calc_substitution(avg_sleep_geo_mean,
                      imp_stacked,
                      model,
                      c("avg_sleep", "avg_mvpa"),
                      timegroup = timegroup)

  full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset")

  return(as.matrix(full_df))
}
