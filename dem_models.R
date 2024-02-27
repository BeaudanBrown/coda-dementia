source("utils.R")

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

get_s3_formula <- function(data) {
  knots_timegroup <- quantile(data[["timegroup"]], c(0.05, 0.275, 0.5, 0.725, 0.95))
  knots_deprivation <- quantile(data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))
  knots_bmi <- quantile(data[["BMI"]], c(0.1, 0.5, 0.9))
  knots_bp <- quantile(data[["bp_syst_avg"]], c(0.1, 0.5, 0.9))
  
  s3_formula <- as.formula(dem ~ rcs(timegroup, knots_timegroup) +
                             (sex + retired + avg_total_household_income + 
                                smok_status)*poly(R1, 2) +
                             (sex + retired + avg_total_household_income + 
                                smok_status)*poly(R2, 2) +
                             (sex + retired + avg_total_household_income + 
                                smok_status)*poly(R3, 2) +
                             shift +
                             apoe_e4 +
                             highest_qual +
                             rcs(townsend_deprivation_index, knots_deprivation) +
                             antidepressant_med +
                             antipsychotic_med +
                             insomnia_med +
                             ethnicity
  )
}

predict_composition_risk <- 
  function(composition, stacked_data_table, model, timegroup, empirical) {
  
  if(isTRUE(empirical)){
    ilr <- ilr(composition, V = v)
    
    stacked_data_table[, c("R1", "R2", "R3") := list(ilr[1], ilr[2], ilr[3])]
    
    stacked_data_table[, haz := predict(model, newdata = .SD, type = "response")]
    
    setkey(stacked_data_table, id, timegroup) # sort and set keys for efficient grouping and joining
    stacked_data_table[, risk := 1 - cumprod(1 - haz), by = id]
    
    risk <- stacked_data_table[timegroup == timegroup, .(mean = mean(risk))]
    
    return(risk)
  } else {
    
    ilr <- ilr(composition, V = v)
    
    stacked_data_table[, c("R1", "R2", "R3") := list(ilr[1], ilr[2], ilr[3])]
    
    # shift covariates to match mean (continuous vars) or probability (categorical vars) 
    # of Schoeler et al pseudo-pop (see paper)
    
    for(i in unique(stacked_data_table$id)){
      stacked_data_table[id==i, sex := sample(c("female","male"), 1, prob = c(0.504,0.496))]
      stacked_data_table[id==i, retired := rbinom(1, 1, prob = 0.193)]
      stacked_data_table[id==i, avg_total_household_income := 
                           sample(c("<18","18-30","31-50","52-100",">100"), 1, 
                                  prob = c(0.264,0.141,0.205,0.145,0.435))]
      stacked_data_table[id==i, smok_status := 
                           sample(c("current","former","never"), 1, prob = c(0.208,0.359,0.433))]
    }

    stacked_data_table[, haz := predict(model, newdata = .SD, type = "response")]
    
    setkey(stacked_data_table, id, timegroup) # sort and set keys for efficient grouping and joining
    stacked_data_table[, risk := 1 - cumprod(1 - haz), by = id]
    
    risk <- stacked_data_table[timegroup == timegroup, .(mean = mean(risk))]
    
    return(risk)
  }
}

calc_substitution <- function(base_comp, imp_stacked_dt, model, substitution, timegroup, empirical) {
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

