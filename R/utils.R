#' Logistic function (sigmoid function)
#'
#' @export
logistic <- function(t) 1 / (1 + exp(-t))

#' Check whether binomial family
is_binomial <- function(family) {
  identical(family, binomial(), ignore.environment = TRUE)
}
