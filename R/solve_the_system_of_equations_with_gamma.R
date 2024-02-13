#' Solve the system of non-linear equations with gamma.
#' Equation [5] as presented by Sur and Candès (2018).
#' For details, see section 4.b.
#'
#' @param gamma_hat numeric. estimated signal strength
#' @param kappa numeric. dimensionality p/n
#' @param initial_values numeric vector of three elements
#' @param tol numeric. tolerance of integral
#'
#' @return solution contains x = [alpha, sigma_squared, lambda]
solve_the_system_of_equations_with_gamma <- function(
    gamma_hat, kappa, initial_values = c(2, 1, 1), tol = 0.08) {

  # covariance matrix (equation [6])
  gamma2 <- gamma_hat^2
  Sigma <- function(alpha, sigma2, kappa) {
    cov <- -alpha * gamma2
    sigma <- c(gamma2, cov, cov, alpha^2 * gamma2 + kappa * sigma2)
    matrix(sigma, nrow = 2, ncol = 2)
  }

  solve_the_system_of_equations(Sigma, kappa, initial_values, tol)
}

