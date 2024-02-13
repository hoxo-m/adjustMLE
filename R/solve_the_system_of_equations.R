#' Solve the system of three non-linear equations.
#' Equation [5] from Sur and Candès (2018).
#' Equation (2.6) from Yadlowsky et al. (2021).
#'
#' @param Sigma function to calculate covariance matrix
#' @param kappa numeric. called "dimensionality" or "aspect ratio"
#' @param initial_values numeric vector of three elements
#' @param tol numeric. tolerance of integral
#'
#' @return solution contains x = [alpha, sigma_squared, lambda]
solve_the_system_of_equations <- function(Sigma, kappa,
                                          initial_values = c(2, 1, 1),
                                          tol = 0.08) {
  # sigmoid function and its derivative and integral
  rho <- logistic_integral
  rho_d <- logistic
  rho_dd <- logistic_derivative

  # proximal mapping operator (equation [7])
  prox <- function(z, lambda) {
    objective <- function(t) lambda * rho(t) + (t - z)^2 / 2
    interval <- c(-abs(z) - lambda - 1, abs(z) + lambda + 1)
    solution <- optimize(objective, interval = interval)
    solution$minimum
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

  # the system of non-linear equations
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

  # solve the equations
  control <- list(ftol = 1e-6)
  solution <- nleqslv::nleqslv(initial_values, equations, control = control)

  if (solution$termcd != 1L) {
    warning("nleqslv() termcd=", solution$termcd, " : ", solution$message)
  }

  solution
}
