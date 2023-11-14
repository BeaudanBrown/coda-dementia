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
  "renv"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

if (length(new_packages)) install.packages(new_packages)

lapply(list_of_packages, library, character.only = TRUE)

# constants
mins_in_day <- 1440
sub_steps <- 12
sub_step_mins <- 5

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


get_primary_formula <- function(data) {
  knots_timegroup <- quantile(data[["timegroup"]], c(0.05, 0.275, 0.5, 0.725, 0.95))
  knots_deprivation <- quantile(data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))

  primary_formula <- as.formula(dem ~ rcs(timegroup, knots_timegroup) +
      poly(R1, 2) +
      poly(R2, 2) +
      poly(R3, 2) +
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
      smok_status
  )

  return(primary_formula)
}


get_s1_formula <- function(data) {
  knots_timegroup <- quantile(data[["timegroup"]], c(0.05, 0.275, 0.5, 0.725, 0.95))
  knots_deprivation <- quantile(data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))
  knots_waso <- quantile(data[["avg_WASO"]], c(0.1, 0.5, 0.9))

  s1_formula <- as.formula(dem ~ rcs(timegroup, knots_timegroup) +
      poly(R1, 2) +
      poly(R2, 2) +
      poly(R3, 2) +
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
      smok_status +
      rcs(avg_WASO, knots_waso)
  )
}

get_s2_formula <- function(data) {
  knots_timegroup <- quantile(data[["timegroup"]], c(0.05, 0.275, 0.5, 0.725, 0.95))
  knots_deprivation <- quantile(data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))
  knots_bmi <- quantile(data[["BMI"]], c(0.1, 0.5, 0.9))
  knots_bp <- quantile(data[["bp_syst_avg"]], c(0.1, 0.5, 0.9))

  s2_formula <- as.formula(dem ~ rcs(timegroup, knots_timegroup) +
      poly(R1, 2) +
      poly(R2, 2) +
      poly(R3, 2) +
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
      smok_status +
      sick_disabled +
      prev_diabetes +
      prev_cancer +
      prev_mental_disorder +
      prev_nervous_system +
      prev_cvd +
      bp_med +
      rcs(BMI, knots_bmi) +
      rcs(bp_syst_avg, knots_bp)
  )
}

predict_composition_risk <- function(composition, stacked_data_table, model, timegroup) {
  ilr <- ilr(composition, V = v)

  stacked_data_table[, c("R1", "R2", "R3") := list(ilr[1], ilr[2], ilr[3])]

  stacked_data_table[, haz := predict(model, newdata = .SD, type = "response")]

  setkey(stacked_data_table, id, timegroup) # sort and set keys for efficient grouping and joining
  stacked_data_table[, risk := 1 - cumprod(1 - haz), by = id]

  risk <- stacked_data_table[timegroup == timegroup, .(mean = mean(risk))]

  return(risk)
}

calc_substitution <- function(base_comp, imp_stacked_dt, model, substitution, timegroup) {

  # The list of substitutions to be calculated in minutes
  inc <- -sub_steps:sub_steps * (sub_step_mins / mins_in_day)

  # Initialize a list to hold generated data.tables
  sub_comps_list <- vector("list", length(inc))

  # Loop over inc and create a sub_comps data table for each element
  for (i in seq_along(inc)) {
    # The list of compositions to be fed into the model after applying the substitutions
    sub_comps <- as.data.table(t(base_comp))
    setnames(sub_comps, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa"))

    sub_comps[, (substitution[1]) := .SD[[substitution[1]]] + inc[i]]
    sub_comps[, (substitution[2]) := .SD[[substitution[2]]] - inc[i]]

    # Store the data.table into list
    sub_comps_list[[i]] <- sub_comps
  }

  # Combine all data.tables in the list
  sub_comps <- rbindlist(sub_comps_list)

  # The risk predicted by the model for each of these composition
  sub_risks <-
    rbindlist(lapply(
                     seq_len(nrow(sub_comps)),
                     function(i) predict_composition_risk(acomp(sub_comps[i]), imp_stacked_dt, model, timegroup))
    )

  result <-
    setnames(
      data.table(offset = inc * mins_in_day, risks = sub_risks),
      c("offset", paste0(substitution[2], "_", deparse(substitute(base_comp))))
    )

  return(result)
}

fit_model <- function(imp, create_formula_fn) {
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

  model <- glm(create_formula_fn(imp_long), data = imp_long, family = binomial)

  return(strip_glm(model))
}

strip_glm <- function(cm) {
  cm$y <- c()
  cm$model <- c()

  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()

  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()

  return(cm)
}

## Impute ##

# imp <- mice(boot_data, predictorMatrix = predmat,
#            method = imp_methods, m = 1)

# saveRDS(imp, file.path(data_dir, "full_imp.rds"))

# imp <- read_rds(file.path(data_dir, "full_imp.rds"))

# imp <- complete(imp)
# imp$id <- seq_len(nrow(imp))

#### Primary model - bootstrap ####

## Get results for substitutions

# primary_model <- fit_model(imp, primary_formula)
# model_s1 <- fit_model(imp, s1_formula)
# model_s2 <- fit_model(imp, s2_formula)

#plot_data <- sub_results(model, imp, timegroup = 55)
# write_rds(plot_data, file.path(data_dir,"primary_mod_res.rds"))
#plot_data <- read_rds(file.path(data_dir, "primary_mod_res.rds"))

## Plot ##

# plot1 <- plot_data |>
#   ggplot(aes(x = offset, y = risk_ratio, colour = Reference)) +
#   geom_line() +
#   facet_wrap(~Substitution) +
#   cowplot::theme_cowplot() +
#   labs(x = "Time added to sleep", y = "Dementia cumulative incidence by age 75")
#
# plot1

#### Figure 1 ####

# combine bootstrap and full sample estimates

#### Sensitivity 1 - full sample ####

#
# ## Get substitution results
#
# plot_data_s1 <- sub_results(model_s1, imp)
#
# plot_s1 <-
#   plot_data_s1 |>
#   ggplot(aes(x = offset, y = risk, colour = Reference)) +
#   geom_line() +
#   facet_wrap(~Substitution) +
#   cowplot::theme_cowplot() +
#   labs(x = "Time added to sleep", y = "Dementia cumulative incidence by age 75") +
#   ylim(c(0, 0.05))
#
# #### Sensitivity 1 - bootstrap ####
#
# ncpus <- 6
# boot_out_s1 <- boot(
#   data = boot_data, statistic = boot_substitutions_fn, reg_formula = s1_formula,
#   R = ncpus, parallel = "multicore", ncpus = ncpus
# )
#
#
# #### Sensitivity 2 - full sample ####
#
#
# ## Get substitution results
#
# plot_data_s2 <- sub_results(model_s2, imp)
#
# #### Sensitivity 2 - bootstrap ####
#
# ncpus <- 6
# boot_out_s2 <- boot(
#   data = boot_data, statistic = boot_substitutions_fn, reg_formula = s2_formula,
#   R = ncpus, parallel = "multicore", ncpus = ncpus
# )
#
#
# #### Sensitivity 3 - full sample ####
#
# ## Create dataset with time since accelerometry as the timescale ##
#
# imp_long <- survSplit(Surv(time = time_to_dem, event = dem) ~ .,
#   data = imp,
#   cut = seq(
#     from = min(imp$time_to_dem),
#     to = max(imp$time_to_dem),
#     length.out = 18
#   ),
#   episode = "timegroup", end = "time_start", event = "dem",
#   start = "time_end"
# )
#
#
# model_s3 <-
#   glm(
#     dem ~ rcs(timegroup, 5) * (poly(R1, 2) + poly(R2, 2) + poly(R3, 2)) +
#       rcs(age_accel, 3) +
#       sex +
#       retired +
#       shift +
#       apoe_e4 +
#       highest_qual +
#       rcs(townsend_deprivation_index, 3) +
#       antidepressant_med +
#       antipsychotic_med +
#       insomnia_med +
#       ethnicity +
#       avg_total_household_income +
#       smok_status,
#     data = imp_long,
#     family = binomial
#   )
