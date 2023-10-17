# load packages
library(mice)
library(tidyverse)
library(survival)
library(rms)
library(compositions)
library(data.table)
library(foreach)
library(doParallel)
library(boot)

mins_in_day <- 1440
sub_steps <- 4
sub_step_mins <- 15
short_sleep_hours <- 6
hrs_in_day <- 24

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

# method for each imputed variable
imp.methods <- make.method(boot_data)

# exclude dates from being imputed
predmat["date_acdem2", ] <- 0
predmat["date_of_death", ] <- 0
imp.methods["date_acdem2"] <- ""
imp.methods["date_of_death"] <- ""

## Impute ##

#imp <- mice(boot_data, predictorMatrix = predmat,
#             method = imp.methods, m = 1)

# saveRDS(imp, file.path(data_dir, "full_imp.rds"))

imp <- read_rds(file.path(data_dir, "full_imp.rds"))

imp <- complete(imp)

# Helper functions for fitting model and calculating risks

fit_model <- function(imp) {
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

  model <-
    glm(dem ~ rcs(timegroup, 5) + poly(R1, 2) + poly(R2, 2) + poly(R3, 2) +
          rcs(bp_syst_avg, 3) + sex + retired + shift + apoe_e4 + highest_qual +
          rcs(townsend_deprivation_index, 3) + antidepressant_med +
          antipsychotic_med + insomnia_med + ethnicity + avg_total_household_income +
          smok_status, family = binomial, data = imp_long)

  return(model)
}

calc_risk <- function(composition, stacked_data, model) {
  ilr <- ilr(composition, V = v)

  ilr_data <-
    mutate(stacked_data,
      R1 = ilr[1], R2 = ilr[2], R3 = ilr[3]
    )

  ilr_data$haz <-
    predict(model, newdata = ilr_data, type = "response")

  ilr_data <- ilr_data |>
    group_by(id) |>
    arrange(timegroup) |>
    mutate(risk = 1 - cumprod(1 - haz)) |>
    ungroup()

  risk <- ilr_data |>
    filter(timegroup == 76) |>
    summarise(mean = mean(risk))

  return(risk)
}

calc_substitution <- function(base_comp, imp_stacked, model, substitution) {
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  sub_comps <- data.frame(matrix(rep(base_comp, length(inc)), nrow = length(inc), byrow = TRUE))
  colnames(sub_comps) <- c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
  sub_comps[, substitution[1]] <- sub_comps[, substitution[1]] + inc
  sub_comps[, substitution[2]] <- sub_comps[, substitution[2]] - inc
  sub_risks <- bind_rows(apply(sub_comps, 1, function(comp) calc_risk(acomp(comp), imp_stacked, model)))

  return(data.frame(
    offset = inc * mins_in_day,
    risks = sub_risks
  ))
}

### FULL SAMPLE ###

imp$id <- seq_len(nrow(imp))

## fit model

model <- fit_model(imp)

# reference composition

all_comp <- acomp(imp[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

short_sleep_comp <-
  all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
short_sleep_geo_mean <-
  acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

avg_sleep_comp <-
  all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
avg_sleep_geo_mean <-
  acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

# stacked dataset for predictions (one row for each follow-up interval)

imp_stacked <- do.call("rbind", replicate(76, imp, simplify = FALSE))

imp_stacked$timegroup <- rep(1:76, nrow(imp))

# substitution

ref_mvpa <- calc_substitution(avg_sleep_geo_mean, imp_stacked, model, 
                              c("avg_sleep", "avg_mvpa"))

ref_mvpa$ratio <-
  ref_mvpa$mean / ref_mvpa$mean[which(ref_mvpa$offset == 0)]
ggplot(ref_mvpa, aes(x=offset, y=ratio)) + geom_line()

### Bootstrap ###

boot_fn <- function(data, indices) {
  this_sample <- data[indices, ]
  imp <- mice(this_sample, maxit = 1, m = 1, predictorMatrix = predmat)
  imp <- complete(imp)
  imp$id <- seq_len(nrow(imp))

  model <- fit_model(imp)

  imp_stacked <- do.call("rbind", replicate(76, imp, simplify = FALSE))

  imp_stacked$age_start <- as.integer(
    rep(
      seq(min(imp$age_dem),
          max(imp$age_dem),
          length.out = 76
        ),
      nrow(imp)
    )
  )

  all_comp <- acomp(imp[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

  short_sleep_comp <-
    all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
  short_sleep_geo_mean <-
    acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

  avg_sleep_comp <-
    all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
  avg_sleep_geo_mean <-
    acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

  short_sleep_inactive <- calc_substitution(short_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_inactivity"))
  short_sleep_light <- calc_substitution(short_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_light"))
  short_sleep_mvpa <- calc_substitution(short_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_mvpa"))

  avg_sleep_inactive <- calc_substitution(avg_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_inactivity"))
  avg_sleep_light <- calc_substitution(avg_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_light"))
  avg_sleep_mvpa <- calc_substitution(avg_sleep_geo_mean, imp_stacked, model, c("avg_sleep", "avg_mvpa"))

  full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
             full_join(short_sleep_mvpa, by = "offset") |>
             full_join(avg_sleep_inactive, by = "offset") |>
             full_join(avg_sleep_light, by = "offset") |>
             full_join(avg_sleep_mvpa, by = "offset") |>
             select(-offset)

  return(as.matrix(full_df))
}

ncpus <- 6
# small_sample <- boot_data[sample(1:nrow(boot_data), 5000), ]
# result <- boot(data = small_sample, statistic = boot_fn, R = ncpus, parallel = "multicore", ncpus = ncpus)
result <- boot(data = boot_data, statistic = boot_fn, R = ncpus, parallel = "multicore", ncpus = ncpus)



