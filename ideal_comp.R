source("dem_models.R")

# load packages
list_of_packages <- c(
  "mvtnorm"
)

new_packages <- list_of_packages[
  !(list_of_packages %in% installed.packages()[, "Package"])
]

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
    sleep = seq(lower["avg_sleep"], upper["avg_sleep"], by = step_size),
    inactive = seq(lower["avg_inactivity"], upper["avg_inactivity"], by = step_size),
    light = seq(lower["avg_light"], upper["avg_light"], by = step_size),
    modvig = seq(lower["avg_mvpa"], upper["avg_mvpa"], by = step_size)
  )

  # Calculate total times
  df$total <- rowSums(df)

  # Filter out rows where total time is not 24 * 60
  df <- df[df$total == mins_in_day, ]
  rownames(df) <- NULL

  return(df)
}

generate_hazards <- function(comps, model, model_formula, ref_row) {
  ilrs <- t(apply(
    comps[, c(1, 2, 3, 4)], 1, function(comp) ilr(acomp(comp), V = v)
  ))

  # Convert contrasts into a data frame
  risks <- apply(ilrs, 1, function(new_ilr) {
    ref_row[c("R1", "R2", "R3")] <- new_ilr
    newx <- model.matrix(model_formula, ref_row)
    predict(model, newdata = newx, type = "response")
  })
  return(risks)
}

get_best_and_worst_comp <- function(df, timegroup_cuts, filename = "best_and_worst") {
  file_path <- file.path(output_dir, paste0(filename, ".rds"))
  if (file.exists(file_path)) {
    rds_data <- read_rds(file_path)
    return(rds_data)
  }

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

  dem_comp_df <- as.data.frame(
    dem_df[, c("avg_sleep", "avg_inactivity", "avg_mvpa", "avg_light")]
  )

  include_cols <- c("avg_sleep", "avg_inactivity", "avg_mvpa")
  mean_vec <- as.numeric(apply(dem_comp_df[, include_cols], 2, mean))
  vcv_mat <- cov(dem_comp_df[, include_cols])
  dem_comp_df$dens <- apply(
    dem_comp_df[, include_cols], 1, dmvnorm,
    mean = mean_vec, sigma = vcv_mat, log = TRUE
  )

  threshold <- quantile(dem_comp_df$dens, probs = c(0.025), na.rm = TRUE)

  quantiles <- apply(dem_comp_df, 2, function(column) {
    quantile(column, probs = c(0.025, 0.975), na.rm = TRUE)
  })

  # Store in "lower" and "upper" vectors
  lower <- quantiles[1, ]
  upper <- quantiles[2, ]

  generated_comps <- generate_compositions(lower, upper)
  generated_comps$dens <-
    apply(generated_comps[, c("sleep", "inactive", "modvig")], 1, dmvnorm,
      mean = mean_vec, sigma = vcv_mat, log = TRUE
    )
  generated_comps <- generated_comps[generated_comps$dens > threshold, ]

  # fit model
  models <- fit_model(
    dem_df, timegroup_cuts, get_primary_formula
  )
  dem_model <- models[["model_dem"]]
  model_formula <- models[["model_formula"]]
  ref_row <- as.data.frame(dem_df[1, ])
  # Arbitrary, since hazards will be equivalently ranked across timegroups
  ref_row$timegroup <- 55

  generated_comps$haz <- generate_hazards(
    generated_comps, dem_model, model_formula, ref_row
  )
  generated_comps <- generated_comps[order(generated_comps$haz), ]
  best_comp <- acomp(generated_comps[1, c(1, 2, 3, 4)])
  worst_comp <- acomp(generated_comps[nrow(generated_comps), c(1, 2, 3, 4)])
  typical_comp <- acomp(generated_comps[generated_comps$dens == max(generated_comps$dens), c(1, 2, 3, 4)])

  result <- as.data.frame(t(rbind(best_comp, worst_comp, typical_comp)))
  colnames(result) <- c("best", "worst", "typical")

  write_rds(result, file_path)
  return(result)
}
