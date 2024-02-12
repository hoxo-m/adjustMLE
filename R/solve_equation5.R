#' Solve the system of nonlinear equations [5].
#' For details, see section 4.b.
#'
#' @param gamma_hat double
#' @param kappa double
#' @param quiet logical
#' @param initial_values numeric vector
#' @param tol double
#'
#' @return alpha, sigma squared, and lambda
#'
solve_equation_5 <- function(gamma_hat, kappa) {
  # covariance matrix (equation [6])
  gamma2 <- gamma_hat^2
  Sigma <- function(alpha, sigma2, kappa) {
    cov <- -alpha * gamma2
    sigma <- c(gamma2, cov, cov, alpha^2 * gamma2 + kappa * sigma2)
    matrix(sigma, nrow = 2, ncol = 2)
  }

  solve_SE_system(Sigma, kappa)
}

