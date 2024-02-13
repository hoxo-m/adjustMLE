#' Estimate corrupted signal strength
#'
#' @param fit glm object
#' @param X input data
#'
#' @return estimated eta squared
estimate_eta_square <- function(fit, X) {
  d <- ncol(X)
  n <- nrow(X)

  g_d <- logistic_derivative

  logit <- predict(fit, type = "link")
  resid <- residuals(fit, type = "response")

  # inverse Hessian matrix (equation (3.3))
  H_factor <- g_d(logit)
  H <- matrix(0, d, d)
  for (i in seq_len(n)) {
    H <- H - H_factor[i] * t(X[i, , drop = FALSE]) %*% X[i, , drop = FALSE]
  }
  H_inv <- solve(H)
  # H_inv <- -vcov(fit)[-1, -1, drop = FALSE]

  U <- apply(X, 1, function(row) t(row) %*% H_inv %*% matrix(row))

  S <- logit + resid * U / (1 + H_factor * U)

  # variance of S (equation (3.4))
  eta2_hat <- mean(S^2) - mean(S)^2
  eta2_hat
}
