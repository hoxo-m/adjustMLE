#' adjust MLE
#'
#' @param fit glm object.
#' @param quiet logical.
#'
#' @export
adjustMLE <- function(fit, method = c("SLOE", "ProbeFrontier"), quiet = FALSE) {
  if (!is_binomial(fit$family)) stop("must glm(family = binomial())")
  if (length(fit$x) == 0) stop("must glm(x = TRUE)")
  if (length(fit$y) == 0) stop("must glm(y = TRUE)")

  method <- match.arg(method)

  X <- fit$x
  X <- X[, colnames(X) != "(Intercept)"]
  y <- fit$y

  p <- ncol(X)
  n <- NROW(y)
  kappa <- p / n

  if (!quiet) message("Checking MLE exists...", appendLF = FALSE)
  result <- detectseparation::detect_separation(X, y, family = binomial())
  is_separable <- result$outcome
  if (is_separable) stop("MLE does not exist.")
  if (!quiet) message("ok")

  if (method == "SLOE") {
    kappa_hat <- gamma_hat <- NULL

    if (!quiet) message("Estimating eta...", appendLF = FALSE)
    eta2_hat <- estimate_eta_square(fit, X)
    if (!quiet) message("done")

    if (!quiet) message("Solving equation (2.6)...", appendLF = FALSE)
    solution <- solve_equation_2_6(eta2_hat, kappa = kappa)$x
    if (!quiet) message("done")
  } else if (method == "ProbeFrontier") {
    if (!quiet) message("Searching kappa_hat...")
    kappa_hat <- search_kappa_hat(X, y, quiet = quiet)

    if (!quiet) message("Estimating gamma...", appendLF = FALSE)
    gamma_hat <- estimate_gamma(kappa_hat)$root
    if (!quiet) message("done")

    if (!quiet) message("Solving equation [5]...", appendLF = FALSE)
    solution <- solve_equation_5(gamma_hat, kappa = kappa)$x
    if (!quiet) message("done")
  }

  alpha_hat <- solution[1]
  sigma_squared_hat <- solution[2]
  lambda_hat <- solution[3]

  # equation [11]
  factor_for_chi_squared <- kappa * sigma_squared_hat / lambda_hat

  # change fit
  fit$coefficients <- fit$coefficients / alpha_hat

  family <- binomial()
  eta <- predict(fit, newdata = data.frame(X))
  mu <- predict(fit, newdata = data.frame(X), type = "response")
  fit$residuals <- (y - mu) / family$mu.eta(eta)

  fit$fitted.values <- mu
  fit$linear.predictors <- eta

  fit$deviance <- sum(family$dev.resids(y, mu, wt = 1))

  fit$aic <- family$aic(y, n, mu, wt = 1, fit$deviance) + 2 * fit$rank

  fit$parameters <- list(
    alpha_hat = alpha_hat, sigma_squared_hat = sigma_squared_hat,
    lambda_hat = lambda_hat, factor_for_chi_squared = factor_for_chi_squared,
    kappa = kappa, kappa_hat = kappa_hat, gamma_hat = gamma_hat
  )

  keep <- c(
    "parameters", "coefficients", "residuals", "fitted.values", "rank",
    "family", "linear.predictors", "deviance", "aic", "null.deviance", "iter",
    "prior.weights", "df.residual", "df.null", "converged", "boundary",
    "call", "y"
  )
  fit_adj <- structure(fit[keep], class = c("adjustMLE", class(fit)))
  fit_adj
}

