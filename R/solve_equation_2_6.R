solve_equation_2_6 <- function(eta2_hat, kappa) {
  # equation (3.1)
  Sigma <- function(alpha, sigma2, kappa) {
    cov <- -(eta2_hat - kappa * sigma2) / alpha
    gamma2 <- -cov/alpha
    sigma <- c(gamma2, cov, cov, eta2_hat)
    # gamma2 <- eta2_hat - kappa * sigma2
    # cov <- - alpha * gamma2
    # sigma <- c(gamma2, cov, cov, eta2_hat)
    matrix(sigma, nrow = 2, ncol = 2)
  }

  solve_SE_system(Sigma, kappa = kappa)
}
