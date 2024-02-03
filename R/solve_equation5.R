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
solve_equation_5 <- function(gamma_hat, kappa,
                             initial_values = c(2, 1, 1), tol = 0.08) {
  # sigmoid function and its derivative and integral
  rho <- function(t) log(1 + exp(t))
  rho_d <- logistic
  rho_dd <- function(t) exp(t) / ((1 + exp(t))^2)

  # proximal mapping operator (equation [7])
  prox <- function(z, lambda) {
    objective <- function(t) lambda * rho(t) + (t - z)^2 / 2
    interval <- c(-abs(z) - lambda - 1, abs(z) + lambda + 1)
    solution <- optimize(objective, interval = interval)
    solution$minimum
  }

  # covariance matrix (equation [6])
  gamma2 <- gamma_hat^2
  Sigma <- function(alpha, sigma2, kappa) {
    cov <- -alpha * gamma2
    sigma <- c(gamma2, cov, cov, alpha^2 * gamma2 + kappa * sigma2)
    matrix(sigma, nrow = 2, ncol = 2)
  }

  # generate probability density function for 2D normal distribution (mean = 0)
  generate_pdf_normal_2d <- function(Sigma) {
    var1 <- Sigma[1, 1]
    var2 <- Sigma[2, 2]
    s1 <- sqrt(var1)
    s2 <- sqrt(var2)
    inv_s1s2 <- 2 / (s1 * s2)
    correlation <- Sigma[1, 2] / (s1 * s2)
    correlation <- correlation |> min(1) |> max(-1)
    term1 <- 1 / (2 * pi * s1 * s2 * sqrt(1 - correlation^2))
    term2_1 <- - 1 / (2 * (1 - correlation^2))
    function(x, y) {
      term2_2 <- x^2/var1 + y^2/var2 - correlation*x*y*inv_s1s2
      term2 <- exp(term2_1 * term2_2)
      term1 * term2
    }
  }

  equation1 <- function(alpha, sigma2, lambda, pdf_normal_2d) {
    integrand <- function(q) {
      value <- (2 * rho_d(q[1]) * (lambda * rho_d(prox(q[2], lambda))) ^ 2)
      value * pdf_normal_2d(q[1], q[2])
    }
    result <- cubature::pcubature(integrand, c(-10, -10), c(10, 10), tol = tol)
    result$integral
  }

  equation2 <- function(alpha, sigma2, lambda, pdf_normal_2d) {
    integrand <- function(q) {
      value <- rho_d(q[1]) * q[1] * lambda * rho_d(prox(q[2], lambda))
      value * pdf_normal_2d(q[1], q[2])
    }
    result <- cubature::pcubature(integrand, c(-10, -10), c(10, 10), tol = tol)
    result$integral
  }

  equation3 <- function(alpha, sigma2, lambda, pdf_normal_2d) {
    integrand <- function(q) {
      value <- 2 * rho_d(q[1]) / (1 + lambda * rho_dd(prox(q[2], lambda)))
      value * pdf_normal_2d(q[1], q[2])
    }
    result <- cubature::pcubature(integrand, c(-10, -10), c(10, 10), tol = tol)
    result$integral
  }

  # the system of nonlinear equations [5]
  equations <- function(params) {
    alpha <- params[1]
    sigma2 <- params[2]
    lambda <- params[3]

    pdf_normal_2d <- generate_pdf_normal_2d(Sigma(alpha, sigma2, kappa))

    eq1 <- equation1(alpha, sigma2, lambda, pdf_normal_2d) - sigma2 * kappa^2
    eq2 <- equation2(alpha, sigma2, lambda, pdf_normal_2d)
    eq3 <- equation3(alpha, sigma2, lambda, pdf_normal_2d) - (1 - kappa)

    c(eq1, eq2, eq3)
  }

  # solve the equations [5]
  control <- list(ftol = 1e-6)
  solution <- nleqslv::nleqslv(initial_values, equations, control = control)

  if (solution$termcd != 1L) {
    warning("nleqslv() termcd=", solution$termcd, " : ", solution$message)
  }

  solution
}

