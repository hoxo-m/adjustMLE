#' Estimate gamma by solving equation [4]
#'
#' @param kappa_hat double
#'
#' @return estimated gamma
#'
#' @export
estimate_gamma <- function(kappa_hat) {
  # expected value in equation [4]
  objective <- function(t, gamma) {
    # probability density function of the standard normal distribution
    pdf_std_normal <- dnorm

    # probability density function of the random variable V in Theorem 1
    pdf_v <- function(v) {
      2 * logistic(gamma * v) * pdf_std_normal(v)
    }

    integrand <- function(zv) {
      z <- zv[1]
      v <- zv[2]
      f <- function(z, v, t) max(z - t * v, 0) ^ 2
      f(z, v, t) * pdf_std_normal(z) * pdf_v(v)
    }

    # integrate to compute expected value
    result <- cubature::hcubature(integrand, c(-10, -10), c(10, 10), tol = 1e-3)
    result$integral
  }

  # inverse of g_MLE (equation [4])
  inverse_g_MLE <- function(gamma) {
    # minimize expected value
    solution <- optimize(objective, gamma = gamma, c(0, 20), tol = 1e-2)
    solution$objective
  }

  # solve equation to estimate gamma
  solution <- uniroot(function(gamma) inverse_g_MLE(gamma) - kappa_hat,
                      tol = 1e-06, interval = c(0, 50))

  solution
}


