#' Search kappa_hat by Probing the MLE Frontier.
#' For details, see section 5.a.
#'
#' @param X data.frame or matrix
#' @param y vector or matrix
#' @param quiet logical
#' @param B integer
#' @param kappa_step double
#' @param seed integer
#' @param solver character
#' @param n_cores integer
#'
#' @return kappa_hat
search_kappa_hat <- function(X, y, quiet = FALSE, B = 50L, kappa_step = 1e-3,
                             seed = NULL, solver = "lpsolve",
                             n_cores = parallel::detectCores() - 1L) {
  p <- ncol(X)
  n <- nrow(X)

  # setting for parallel computation
  if (n_cores <= 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession)
  }

  # initialize
  kappa_min <- p / n
  kappa_max <- 0.5
  pi_hat_min <- 0
  pi_hat_max <- 1

  if (!is.null(seed)) set.seed(seed)

  # binary search until kappa_max - kappa_min <= kappa_step
  n_search <- ceiling(log2(kappa_max - kappa_min) - log2(kappa_step))

  if (!quiet) pb <- txtProgressBar(style = 3, min = 0, max = n_search * B)
  is_speed_over_accuracy <- TRUE
  for (i in seq_len(n_search)) {
    kappa_j <- (kappa_max - kappa_min) / 2 + kappa_min
    n_j <- round(p / kappa_j)

    n_trials <- 0L
    n_success <- 0L
    while (n_trials < B) {
      is_separatable_list <- list()
      S <- if (is_speed_over_accuracy) n_cores else B
      for (s in seq_len(S)) {
        # subsample
        index <- sample(seq_len(n), n_j)
        X_subsample <- X[index, , drop = FALSE]
        y_subsample <- y[index]

        # check whether MLE exists
        is_separatable_list[[s]] <- future::future(seed = TRUE, {
          control <- list(solver = solver)
          result <- detectseparation::detect_separation(
            X_subsample, y_subsample, family = binomial(), control = control)
          result$outcome
        })
      }
      is_separatable_list <- Map(future::value, is_separatable_list)
      n_success <- n_success + sum(unlist(is_separatable_list))
      n_trials <- n_trials + S

      # If there is a high probability that pi_hat < 0.5 or pi_hat > 0.5,
      # prioritize speed over accuracy.
      if (is_speed_over_accuracy) {
        binom <- binom.test(n_success, n_trials)
        if (binom$conf.int[1] > 0.5 || binom$conf.int[2] < 0.5) {
          break
        }
        if (!quiet) setTxtProgressBar(pb, (i-1) * B + n_trials)
      }
    }
    # Once over B, keep to prioritize accuracy.
    if (n_trials >= B) {
      is_speed_over_accuracy <- FALSE
    }
    # percentage of MLE that does not exist
    pi_hat <- n_success / n_trials

    # narrow down the search range
    if (pi_hat < 0.5) {
      kappa_min <- kappa_j
      pi_hat_min <- pi_hat
    } else {
      kappa_max <- kappa_j
      pi_hat_max <- pi_hat
    }
    if (!quiet) setTxtProgressBar(pb, i * B)
  }
  if (!quiet) close(pb)

  # linear interpolation
  ratio <- (0.5 - pi_hat_min) / (pi_hat_max - pi_hat_min)
  kappa_hat <- ratio * (kappa_max - kappa_min) + kappa_min

  kappa_hat
}
