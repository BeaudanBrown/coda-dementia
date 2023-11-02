source("dem_models.R")

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24
ncpus <- 5
bootstrap_iterations <- ncpus

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))
# set date variables to strings to avoid errors
boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)
## UNCOMMENT FOR TEST RUNS
# boot_data <- boot_data[1:1000, ]

# Run bootstrap for primary model
run_primary_bootstrap <- function() {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    reg_formula = primary_formula,
    output_name = "boot_primary.rds"
  )
}

# Run bootstrap for sensitivity analysis model 1
run_s1_bootstrap <- function(boot_data) {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    reg_formula = s1_formula,
    output_name = "boot_s1.rds"
  )
}

# Run bootstrap for sensitivity analysis model 2
run_s2_bootstrap <- function(boot_data) {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    reg_formula = s2_formula,
    output_name = "boot_s2.rds"
  )
}

# Run bootstrap for a particular model formula and timegroup target, outputting to a file
run_bootstrap <- function(boot_data, timegroup, reg_formula, output_name) {
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

  result <- boot(
    data = boot_data,
    statistic = bootstrap_substitutions_fn,
    reg_formula = reg_formula,
    timegroup = timegroup,
    predmat = predmat,
    imp_methods = imp_methods,
    short_sleep_geo_mean = short_sleep_geo_mean,
    avg_sleep_geo_mean = avg_sleep_geo_mean,
    R = bootstrap_iterations,
    parallel = "multicore",
    ncpus = ncpus
  )

  # Prepend timestamp to avoid accidental data loss
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  output_name_with_timestamp <- paste(timestamp, output_name, sep = "_")
  saveRDS(result, file.path(output_dir, output_name_with_timestamp))
}

# Ran each bootstrap iteration
# Imputes the data, fits the model for the given formula and predicts the risks for a
# variety of substitutions and cohorts
bootstrap_substitutions_fn <- function(
  data,
  indices,
  reg_formula,
  timegroup,
  predmat,
  imp_methods,
  short_sleep_geo_mean,
  avg_sleep_geo_mean
) {
  this_sample <- data[indices, ]

  # imputation
  imp <- mice(this_sample, m = 1,
              predictorMatrix = predmat,
              methods = imp_methods)

  imp <- complete(imp)
  imp_len <- nrow(imp)
  imp$id <- seq_len(imp_len)

  # fit model
  model <- fit_model(imp, reg_formula)

  # data for g-computation/standardisation
  imp <- do.call("rbind", replicate(timegroup, imp, simplify = FALSE))
  imp$timegroup <- rep(1:timegroup, imp_len)
  setDT(imp) # convert imp to a data.table

  # substitutions for each reference comp
  short_sleep_inactive <-
    calc_substitution(short_sleep_geo_mean,
                      imp,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  short_sleep_light <-
    calc_substitution(short_sleep_geo_mean,
                      imp,
                      model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  short_sleep_mvpa <-
    calc_substitution(short_sleep_geo_mean,
                      imp,
                      model,
                      c("avg_sleep", "avg_mvpa"),
                      timegroup = timegroup)

  avg_sleep_inactive <-
    calc_substitution(avg_sleep_geo_mean,
                      imp,
                      model,
                      c("avg_sleep", "avg_inactivity"),
                      timegroup = timegroup)

  avg_sleep_light <-
    calc_substitution(avg_sleep_geo_mean,
                      imp,
                      model,
                      c("avg_sleep", "avg_light"),
                      timegroup = timegroup)

  avg_sleep_mvpa <-
    calc_substitution(avg_sleep_geo_mean,
                      imp,
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


## Return results from substitutions

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
