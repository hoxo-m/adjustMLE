#' Solve the system of non-linear equations with gamma.
#' Equation (2.6) from Yadlowsky et al. (2021).
#' For details, see section 2.2.
#'
#' @param eta_hat numeric. estimated corrupted signal strength
#' @param kappa numeric. aspect ratio d/n
#' @param initial_values numeric vector of three elements
#' @param tol numeric. tolerance of integral
#'
#' @return solution contains x = [alpha, sigma_squared, lambda]
solve_the_system_of_equations_with_eta <- function(
    eta2_hat, kappa, initial_values = c(2, 1, 1), tol = 0.08) {

  # covariance matrix (equation (3.1))
  Sigma <- function(alpha, sigma2, kappa) {
    cov <- -(eta2_hat - kappa * sigma2) / alpha
    gamma2 <- -cov/alpha
    sigma <- c(gamma2, cov, cov, eta2_hat)
    matrix(sigma, nrow = 2, ncol = 2)
  }

  solve_the_system_of_equations(Sigma, kappa, initial_values, tol)
}
