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
