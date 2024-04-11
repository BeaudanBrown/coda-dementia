source("dem_models.R")

# load packages
list_of_packages <- c(
  "mvtnorm"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]

if (length(new_packages)) install.packages(new_packages)

lapply(list_of_packages, library, character.only = TRUE)

generate_compositions <- function(lower, upper) {
  step_size <- 15

  # Round down for the lower bounds
  lower <- floor(lower / step_size) * step_size
  # Round up for the upper bounds
  upper <- ceiling(upper / step_size) * step_size

  # Create an empty data frame
  df <- expand.grid(
    sleep = seq(lower[1], upper[1], by = step_size),
    inactive = seq(lower[2], upper[2], by = step_size),
    light = seq(lower[3], upper[3], by = step_size),
    modvig = seq(lower[4], upper[4], by = step_size)
  )

  # Calculate total times
  df$total <- rowSums(df)

  # Filter out rows where total time is not 24 * 60
  df <- df[df$total == mins_in_day, ]
  rownames(df) <- NULL

  return(df)
}

generate_hazards <- function(comps, model, ref_row) {
  ilrs <- t(apply(comps[, c(1, 2, 3, 4)], 1, function(comp) ilr(acomp(comp), V = v)))

  # Convert contrasts into a data frame
  risks <- apply(ilrs, 1, function(new_ilr) {
    ref_row[c("R1", "R2", "R3")] <- new_ilr
    predict(model, newdata = ref_row, type = "response")
  })
  return(risks)
}

get_best_and_worst_comp <- function(df) {
  ## Load data
  predmat <- quickpred(df,
    mincor = 0,
    exclude = c(
      "avg_sleep", "avg_inactivity", "avg_light",
      "avg_mvpa", "eid", "time_to_dem"
    )
  )
  dem_df <- mice(df, m = 1, predictorMatrix = predmat, maxit = maxit)
  dem_df <- complete(dem_df)

  dem_comp_df <- as.data.frame(acomp(dem_df[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")]))
  dem_comp_model <-
    lm(cbind(dem_comp_df$avg_sleep, dem_comp_df$avg_inactivity, dem_comp_df$avg_light) ~ 1, data = dem_comp_df)
  vcv_mat <- vcov(dem_comp_model)
  dem_comp_df$dens <- apply(dem_comp_df[, -4], 1, dmvnorm, mean = coef(dem_comp_model), sigma = vcv_mat, log = TRUE)
  sorted_df <- dem_comp_df[order(dem_comp_df$dens, decreasing = TRUE), ]
  sorted_df[, 1:4] <- sorted_df[, 1:4] * 24
  most_common_comp <- acomp(sorted_df[1, 1:4])

  threshold <- quantile(dem_comp_df$dens, probs = c(0.025), na.rm = TRUE)

  quantiles <- apply(dem_comp_df, 2, function(column) {
    quantile(column, probs = c(0.025, 0.975), na.rm = TRUE)
  })

  # Convert into minutes
  quantiles_in_minutes <- quantiles * 24 * 60

  # Store in "lower" and "upper" vectors
  lower <- quantiles_in_minutes[1, ]
  upper <- quantiles_in_minutes[2, ]

  generated_comps <- generate_compositions(lower, upper)
  generated_comps <- generated_comps / mins_in_day
  generated_comps$dens <-
    apply(generated_comps[, -c(4, 5)], 1, dmvnorm, mean = coef(dem_comp_model), sigma = vcv_mat, log = TRUE)
  generated_comps <- generated_comps[generated_comps$dens > threshold, ]

  # fit model
  dem_model <- fit_model(dem_df, get_primary_formula)[["model_dem"]]
  ref_row <- as.data.frame(dem_df[1, ])
  ref_row$timegroup <- 55

  generated_comps$haz <- generate_hazards(generated_comps, dem_model, ref_row)
  generated_comps <- generated_comps[order(generated_comps$haz), ]
  best_comp <- acomp(generated_comps[1, c(1, 2, 3, 4)])
  worst_comp <- acomp(generated_comps[nrow(generated_comps), c(1, 2, 3, 4)])

  result <- as.data.frame(t(rbind(best_comp, worst_comp, most_common_comp)))
  colnames(result) <- c("best", "worst", "most_common")
  return(result)
}
