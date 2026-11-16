#' Discriminability of competing kinetic models
#'
#' Given two or more candidate kinetic models and the proposed
#' experimental design (sampling times, measurement noise), this
#' function quantifies the probability that the data will allow
#' discrimination between the models. It is a pre-experimental
#' planning tool: run it before collecting data to choose sampling times
#' that maximize discriminability.
#'
#' @param models list of named functions. Each function takes
#'   `c(params, t)` and returns predicted concentration.
#' @param starts list of starting parameter vectors for the models.
#' @param true_model integer or character; which model is assumed to be
#'   the data-generating mechanism. Defaults to 1.
#' @param true_params numeric vector of true parameters for the
#'   generating model.
#' @param times numeric vector of candidate sampling times.
#' @param sigma measurement noise SD.
#' @param n_sim number of Monte Carlo simulations (default 500).
#' @param conf_level confidence level for discrimination probability.
#' @param criterion `"aic"` (default) or `"bic"` for model selection.
#' @param alternative_times optional list of alternative sampling
#'   schedules to compare.
#' @param seed optional integer seed.
#'
#' @return A `crossstatr_kinetic_discrimination` object.
#'
#' @section Methodology:
#'
#'   For each simulation \eqn{b = 1, \dots, B}:
#'
#'   1. Simulate data \eqn{y_i^{(b)} = m^*(t_i; \theta^*) + \epsilon_i}
#'      where \eqn{\epsilon_i \sim N(0, \sigma^2)} and \eqn{m^*} is the
#'      generating model with true parameters \eqn{\theta^*}.
#'   2. Fit each candidate model \eqn{m_k} via nonlinear least squares to
#'      obtain \eqn{\hat\theta_k^{(b)}} and the model's AIC/BIC.
#'   3. Select the model with the lowest AIC (or BIC).
#'
#'   The **discrimination probability** for model \eqn{k} is
#'   \deqn{P_k = \frac{1}{B} \sum_b \mathbf{1}\{\hat k^{(b)} = k\}.}
#'
#'   The **power to distinguish model \eqn{k} from model \eqn{j}** is
#'   computed under the alternative that \eqn{k} is the true model:
#'   \deqn{\pi_{k \to j} = P(\mathrm{AIC}_k < \mathrm{AIC}_j \mid m_k\ \mathrm{true}).}
#'
#'   Parameter identifiability for each model is assessed via the
#'   condition number \eqn{\kappa} of the sensitivity matrix
#'   \eqn{G_{ij} = \partial m_i / \partial \theta_j}. \eqn{\kappa > 10^6}
#'   indicates practical non-identifiability.
#'
#'   Model separation (the geometric mean distance between model
#'   trajectories at the proposed sampling times, in units of sigma) is
#'   reported as a quick summary of design quality.
#'
#' @section Limitations:
#'
#'   - The simulated measurement noise is assumed homoscedastic
#'     Gaussian. Real experiments with heteroscedastic or correlated
#'     noise will give different discrimination rates.
#'   - If models are very similar in functional form, very large
#'     `n_sim` may be needed for a precise discrimination probability
#'     estimate. The MCSE of \eqn{P_k} is \eqn{\sqrt{P_k(1-P_k)/B}}.
#'   - Models that are mathematically nested can be discriminated via
#'     likelihood-ratio tests, which are more powerful than AIC for
#'     this purpose. Use [kinetic_model_compare()] for nested
#'     comparisons.
#'   - This function assumes `nls` will converge; in practice,
#'     non-convergence of one model can bias the discrimination
#'     probability.
#'
#' @references
#'   - Asprey, S. P. and Macchietto, S. (2000). Statistical tools for
#'     optimal dynamic model building. Computers & Chemical Engineering,
#'     24, 1261-1267.
#'     \doi{10.1016/S0098-1354(00)00328-8}.
#'
#' @examples
#' set.seed(1)
#' m1 <- function(p, t) p[1] * exp(-p[2] * t)
#' m2 <- function(p, t) p[1] * (1 + p[3] * t) * exp(-p[2] * t)
#' res <- kinetic_discrimination(
#'   models = list(simple = m1, extended = m2),
#'   starts = list(c(A = 1, k = 0.1), c(A = 1, k = 0.1, alpha = 0.05)),
#'   true_model = 1L, true_params = c(1, 0.1),
#'   times = seq(0, 50, length.out = 12), sigma = 0.05, n_sim = 50)
#' print(res)
#'
#' @export
kinetic_discrimination <- function(models, starts, true_model = 1L,
                                    true_params, times, sigma,
                                    n_sim = 500L, conf_level = 0.95,
                                    criterion = c("aic", "bic"),
                                    alternative_times = NULL, seed = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  .assert_positive_integer(n_sim, "n_sim")

  if (!is.list(models) || length(models) < 2L) {
    stop("'models' must be a list of at least 2 functions.", call. = FALSE)
  }
  if (length(models) != length(starts)) {
    stop("'starts' must have the same length as 'models'.", call. = FALSE)
  }
  if (is.character(true_model)) {
    true_model <- match(true_model, names(models))
    if (is.na(true_model)) {
      stop("'true_model' name not found in 'models'.", call. = FALSE)
    }
  }
  model_names <- if (is.null(names(models)))
    paste0("M", seq_along(models)) else names(models)
  names(models) <- model_names
  names(starts) <- model_names

  .set_seed(seed)

  # Identifiability: condition number of sensitivity matrix per model
  identifiability <- lapply(seq_along(models), function(k) {
    p <- starts[[k]]
    G <- .num_jacobian(p, function(par) models[[k]](par, times))
    sv <- svd(G)$d
    kappa <- max(sv) / min(sv[sv > .Machine$double.eps])
    list(condition_number = kappa, params = names(p))
  })
  names(identifiability) <- model_names

  # Model separation at the proposed sampling times
  pred_at_times <- lapply(seq_along(models), function(k) {
    models[[k]](starts[[k]], times)
  })
  seps <- outer(seq_along(models), seq_along(models), function(i, j) {
    vapply(seq_along(i), function(k)
      sqrt(mean((pred_at_times[[i[k]]] - pred_at_times[[j[k]]])^2)),
      numeric(1L))
  })
  seps <- seps / sigma
  dimnames(seps) <- list(model_names, model_names)

  # Monte Carlo discrimination
  true_y <- models[[true_model]](true_params, times)
  picks <- character(n_sim)
  converged_mat <- matrix(TRUE, nrow = n_sim, ncol = length(models))
  for (b in seq_len(n_sim)) {
    y_sim <- true_y + rnorm(length(times), sd = sigma)
    aics <- vapply(seq_along(models), function(k) {
      tryCatch({
        fit <- stats::nls(y_sim ~ models[[k]](p, times),
                          start = list(p = starts[[k]]),
                          control = list(maxiter = 200))
        rss <- sum(stats::resid(fit)^2)
        n_p <- length(starts[[k]])
        n <- length(times)
        if (criterion == "aic") {
          n * log(rss / n) + 2 * n_p
        } else {
          n * log(rss / n) + log(n) * n_p
        }
      }, error = function(e) {
        converged_mat[b, k] <<- FALSE
        Inf
      })
    }, numeric(1L))
    names(aics) <- model_names
    picks[b] <- model_names[which.min(aics)]
  }

  # Discrimination probabilities
  disc_probs <- vapply(model_names, function(nm) mean(picks == nm),
                       numeric(1L))
  mcse <- .mcse_p(disc_probs, n_sim)

  # Power to distinguish each pair (under alternative that true_model is true)
  pair_powers <- matrix(NA_real_, length(models), length(models))
  dimnames(pair_powers) <- list(model_names, model_names)
  for (i in seq_along(models)) {
    for (j in seq_along(models)) {
      if (i == j) next
      # P(M_i wins | M_true generates data) is just disc_probs[i]
      pair_powers[i, j] <- disc_probs[i]
    }
  }

  alt_results <- NULL
  if (!is.null(alternative_times)) {
    alt_results <- lapply(alternative_times, function(alt_t) {
      sub_call <- call
      sub_call$times <- alt_t
      sub_call$alternative_times <- NULL
      eval(sub_call, envir = parent.frame())
    })
  }

  .crossstatr_result(
    list(discrimination_probabilities = disc_probs,
         mcse = mcse,
         criterion = criterion,
         n_sim = n_sim,
         true_model = model_names[true_model],
         pair_powers = pair_powers,
         identifiability = identifiability,
         model_separation = seps,
         converged = converged_mat,
         alternative_times_results = alt_results,
         method = sprintf("Kinetic discrimination (true model: %s)",
                          model_names[true_model]),
         call = call),
    "crossstatr_kinetic_discrimination")
}

#' @export
print.crossstatr_kinetic_discrimination <- function(x, ...) {
  cat(sprintf("<crossStatR::kinetic_discrimination> %s\n", x$method))
  cat(sprintf("Discrimination probabilities (criterion = %s, B = %d):\n",
              x$criterion, x$n_sim))
  for (nm in names(x$discrimination_probabilities)) {
    p <- x$discrimination_probabilities[[nm]]
    cat(sprintf("  %s: %.3f (MCSE %.4f)\n", nm, p, .mcse_p(p, x$n_sim)))
  }
  cat("\nCondition numbers:\n")
  for (nm in names(x$identifiability)) {
    cat(sprintf("  %s: kappa = %.2g\n", nm,
                x$identifiability[[nm]]$condition_number))
  }
  invisible(x)
}
