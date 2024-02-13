#' Logistic function (a.k.a. sigmoid function)
logistic <- plogis

#' Derivative of logistic function
logistic_derivative <- dlogis

#' Integral of logistic function
logistic_integral <- function(t) log(1 + exp(t))

#' Check whether binomial family
is_binomial <- function(family) {
  identical(family, binomial(), ignore.environment = TRUE)
}
