dem_comp_df <- as.data.frame(dem_comp)
dem_comp_model <- lm(cbind(dem_df.sleep_n, dem_df.inactive_n, dem_df.light_n) ~ 1, data = dem_comp_df)
vcv_mat <- vcov(dem_comp_model)
dem_comp_df$dens <- apply(dem_comp_df[, -4], 1, dmvnorm, mean = coef(dem_comp_model), sigma = vcv_mat, log = TRUE)
sorted_df <- dem_comp_df[order(dem_comp_df$dens, decreasing = TRUE), ]
sorted_df[, 1:4] <- sorted_df[, 1:4] * 24

threshold <- quantile(dem_comp_df$dens, probs = c(0.025), na.rm = TRUE)

step_size <- 15

quantiles <- apply(dem_comp, 2, function(column) {
  quantile(column, probs = c(0.025, 0.975), na.rm = TRUE)
})

# Convert into minutes
quantiles_in_minutes <- quantiles * 24 * 60

# Store in "lower" and "upper" vectors
lower <- quantiles_in_minutes[1, ]
upper <- quantiles_in_minutes[2, ]

generate_compositions <- function(lower, upper) {
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

generated_comps <- generate_compositions(lower, upper)
generated_comps <- generated_comps / mins_in_day
generated_comps$dens <- apply(generated_comps[, -c(4, 5)], 1, dmvnorm, mean = coef(dem_comp_model), sigma = vcv_mat, log = TRUE)
generated_comps <- generated_comps[generated_comps$dens > threshold, ]

generate_hazards <- function(comps, base_comp, base_ilr) {
  ilrs <- t(apply(comps, 1, function(comp) ilr(acomp(comp), V = v)))

  # Convert contrasts into a data frame
  contrasts <- t(apply(ilrs, 1, function(new_ilr) {
    contrast_out <- contrast(
      dem_model,
      list(R1 = new_ilr[1], R2 = new_ilr[2], R3 = new_ilr[3]),
      list(R1 = base_ilr[1], R2 = base_ilr[2], R3 = base_ilr[3])
    )
    return(c(contrast_out$Contrast, contrast_out$Lower, contrast_out$Upper))
  }))
  contrasts_df <- as.data.frame(contrasts)
  colnames(contrasts_df) <- c("Contrast", "Lower", "Upper")

  # Add contrasts to comps data frame
  comps$Contrast <- contrasts_df$Contrast
  comps$Lower <- contrasts_df$Lower
  comps$Upper <- contrasts_df$Upper

  return(comps)
}

contrasts <- generate_hazards(generated_comps[, c(1, 2, 3, 4)], avg_sleep_geo_mean, dem_base_ilr)
sorted_contrasts <- contrasts[order(-contrasts$Contrast), ]
best_comp <- sorted_contrasts[nrow(sorted_contrasts), c(1, 2, 3, 4)]
worst_comp <- sorted_contrasts[1, c(1, 2, 3, 4)]
best_ilr <- ilr(acomp(best_comp))
worst_ilr <- ilr(acomp(worst_comp))

# Predicting average tbv for best and worst compositions
best_df <- mri_model_data
best_df$R1 <- rep(best_ilr[1], times = nrow(best_df))
best_df$R2 <- rep(best_ilr[2], times = nrow(best_df))
best_df$R3 <- rep(best_ilr[3], times = nrow(best_df))
best_results <- predict(tbv_model, newdata = best_df)
best_avg <- mean(best_results, na.rm = TRUE)

worst_df <- mri_model_data
worst_df$R1 <- rep(worst_ilr[1], times = nrow(worst_df))
worst_df$R2 <- rep(worst_ilr[2], times = nrow(worst_df))
worst_df$R3 <- rep(worst_ilr[3], times = nrow(worst_df))
worst_results <- predict(tbv_model, newdata = worst_df)
worst_avg <- mean(worst_results, na.rm = TRUE)

