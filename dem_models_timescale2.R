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


predict_composition_risk <- function(composition, stacked_data_table, model) {
  ilr <- ilr(composition, V = v)
  stacked_data_table[, c("R1", "R2", "R3") := list(ilr[1], ilr[2], ilr[3])]

  stacked_data_table[, haz :=  predict(model, newdata = .SD, type = "response")]
  print(stacked_data_table)
  
  return(stacked_data_table$haz)
  
}

calc_substitution <- function(base_comp, imp_stacked_dt, model, substitution) {

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
    rbindlist(seq_len(nrow(sub_comps)),
              function(i) predict_composition_risk(acomp(sub_comps[i]), imp_stacked_dt, model))

  result <-
    setnames(
      data.table(offset = inc * mins_in_day, risks = sub_risks),
      c("offset", paste0(substitution[2], "_", deparse(substitute(base_comp))))
    )

  return(result)
}


get_primary_formula <- function(data) {
  knots_timegroup <- quantile(data[["timegroup"]], c(0.05, 0.275, 0.5, 0.725, 0.95))
  knots_deprivation <- quantile(data[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))
  
  primary_formula <- as.formula(dem ~ 
                                  rcs(timegroup, knots_timegroup)*poly(R1, 2) +
                                  rcs(timegroup, knots_timegroup)*poly(R2, 2) +
                                  rcs(timegroup, knots_timegroup)*poly(R3, 2) + 
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

fit_model <- function(imp, create_formula_fn) {
  imp_long <- survSplit(Surv(time = time_to_dem, event = dem) ~ .,
                        data = imp,
                        cut = seq(
                          from = min(imp$age_dem),
                          to = max(imp$age_dem),
                          length.out = 34
                        ),
                        episode = "timegroup", end = "time_start", event = "dem",
                        start = "time_end"
  )
  
  model <- glm(get_primary_formula(imp_long), data = imp_long, family = binomial)
  
  return(model)
}

