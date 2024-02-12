#' Logistic function (sigmoid function)
#'
#' @export
logistic <- plogis

#' Derivative of logistic function
#'
#' @export
logistic_derivative <- dlogis

#' Derivative of logistic function
#'
#' @export
logistic_integral <- function(t) log(1 + exp(t))

# Check whether binomial family
is_binomial <- function(family) {
  identical(family, binomial(), ignore.environment = TRUE)
}
