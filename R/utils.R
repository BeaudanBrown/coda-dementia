Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

bootstrap_sample <- function(data) {
  out <- data[sample(nrow(data), nrow(data), TRUE)]
  out[, id := .I]
  return(out)
}

## Define SBP
sbp <- matrix(
  c(
    1,
    1,
    -1,
    -1,
    1,
    -1,
    0,
    0,
    0,
    0,
    1,
    -1
  ),
  ncol = 4,
  byrow = TRUE
)

v <- compositions::gsi.buildilrBase(t(sbp))


## strip unneccessary model components

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

save_plot <- function(plot, file_path) {
  ggsave(
    file_path,
    plot = plot,
    device = "svg",
    width = 10,
    height = 12
  )
}

ordinal_to_numeric <- function(data) {
  ## set binary/ordinal factor variables to numeric to allow use of pmm
  data$ethnicity <- as.numeric(data$ethnicity)
  data$avg_total_household_income <- ifelse(
    data$avg_total_household_income == "<18",
    0,
    ifelse(
      data$avg_total_household_income == "18-30",
      1,
      ifelse(data$avg_total_household_income == "31-50", 2, 3)
    )
  )
  data$highest_qual <- ifelse(
    data$highest_qual == "O",
    0,
    ifelse(
      data$highest_qual == "NVQ",
      1,
      ifelse(data$highest_qual == "A", 2, 3)
    )
  )
  data$smok_status <- as.numeric(data$smok_status) - 1
  data$bp_med <- as.numeric(data$bp_med) - 1
  data$any_cvd <- as.numeric(data$any_cvd) - 1
  data$chronotype <- as.numeric(data$chronotype) - 1
  data$apoe_e4 <- as.numeric(data$apoe_e4) - 1

  data$freq_depressed_twoweeks <- ifelse(
    data$freq_depressed_twoweeks == "not at all",
    0,
    ifelse(
      data$freq_depressed_twoweeks == "several days",
      1,
      ifelse(data$freq_depressed_twoweeks == "more than half", 2, 3)
    )
  )
  data$diagnosed_diabetes <- as.numeric(data$diagnosed_diabetes) - 1

  data
}
