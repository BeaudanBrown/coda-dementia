#### Primary

## Outcome model
SL.glm.Q <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  # make model matrix
  X <- cbind(
    Y = X$Y,
    model.matrix(get_primary_outcome_formula(X), data = X)[, -1]
  ) |>
    as.data.frame()
  fit.glm <- glm(
    Y ~ .,
    data = X,
    family = family,
    weights = obsWeights,
    model = model
  )
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  newX <- model.matrix(get_primary_outcome_formula(newX), data = newX)[,
    -1
  ] |>
    as.data.frame()
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm.Q"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glm.Q <- function(object, newdata, ...) {
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  newdata <- model.matrix(
    get_primary_outcome_formula(newdata),
    data = newdata
  )[,
    -1
  ] |>
    as.data.frame()

  pred <- predict(
    object = object$object,
    newdata = newdata,
    type = "response"
  )
  pred
}

## Treatment model

SL.glm.g <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  # make model matrix
  X <- cbind(
    Y = X$Y,
    model.matrix(get_primary_treatment_formula(X), data = X)[, -1]
  ) |>
    as.data.frame()
  fit.glm <- glm(
    Y ~ .,
    data = X,
    family = family,
    weights = obsWeights,
    model = model
  )
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  newX <- model.matrix(get_primary_treatment_formula(newX), data = newX)[,
    -1
  ] |>
    as.data.frame()
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm.g"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glm.g <- function(object, newdata, ...) {
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  newdata <- model.matrix(
    get_primary_treatment_formula(newdata),
    data = newdata
  )[,
    -1
  ] |>
    as.data.frame()

  pred <- predict(
    object = object$object,
    newdata = newdata,
    type = "response"
  )
  pred
}
