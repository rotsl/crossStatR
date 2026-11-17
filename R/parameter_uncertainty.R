#' Propagate parameter uncertainty to derived quantities
#'
#' Given a parameter estimate and its variance-covariance matrix, propagate
#' uncertainty through a function `func(params) -> g(params)` to obtain a
#' point estimate, standard error, and confidence interval for `g`. Three
#' methods are offered: the delta method (linearisation), Monte Carlo
#' simulation from a multivariate normal, and a nonparametric bootstrap that
#' resamples the original data and refits.
#'
#' The delta method is fast and asymptotically exact, but only accurate when
#' `func` is well-approximated by a linear function in the neighbourhood of
#' the parameter estimate. Monte Carlo captures nonlinearity but assumes
#' the parameter estimate is approximately Gaussian (the same assumption
#' that motivates the t-based CI of the original fit). Bootstrap makes the
#' fewest distributional assumptions but requires the original data and a
#' refit function.
#'
#' @param params numeric vector of parameter estimates.
#' @param cov variance-covariance matrix of `params` (square, same length).
#' @param func function taking `params` and returning a numeric scalar or
#'   vector. The function defines the derived quantity `g`.
#' @param method `"delta"`, `"monte_carlo"`, or `"bootstrap"`.
#' @param n_sim integer; number of Monte Carlo or bootstrap replicates.
#' @param conf_level confidence level for the interval.
#' @param df positive degrees of freedom for delta-method intervals. The
#'   default `Inf` uses a normal critical value; residual degrees of freedom
#'   from a fitted model can be supplied for a Student-t approximation.
#' @param seed optional integer seed.
#' @param original_data required for `method = "bootstrap"`: a list with
#'   elements `data` (the original data frame) and `refit` (a function
#'   taking a data subset and returning a parameter vector of the same
#'   length as `params`). If `NULL` and `method = "bootstrap"`, the
#'   function falls back to the delta method with a warning.
#'
#' @return A `crossstatr_param_uncert` object with elements:
#'   \item{estimate}{numeric vector (or scalar) of point estimates of `g`.}
#'   \item{se}{standard error(s) of `g`.}
#'   \item{ci_low, ci_high}{confidence interval bounds.}
#'   \item{method}{method label.}
#'   \item{jacobian}{the Jacobian \eqn{J = \partial g / \partial\theta}
#'     (only for `"delta"`); `NULL` otherwise.}
#'   \item{sim_samples}{matrix of simulated `g` values
#'     (rows = replicates, columns = output components); `NULL` for delta.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Delta method (Oehlert, 1992). For a
#'   differentiable function \eqn{g: \mathbb{R}^p \to \mathbb{R}^q}, the
#'   first-order Taylor expansion around \eqn{\hat{\theta}} gives
#'   \deqn{g(\hat{\theta}) \approx g(\theta) + J(\theta)'(\hat{\theta} - \theta),}
#'   where \eqn{J_{ij} = \partial g_i / \partial \theta_j} is computed
#'   numerically by central finite differences. The variance of
#'   \eqn{g(\hat{\theta})} is approximated by
#'   \deqn{\mathrm{Cov}(g(\hat{\theta})) \approx J\,\mathrm{Cov}(\hat{\theta})\,J'.}
#'   Delta-method intervals use a normal critical value by default. A finite
#'   `df` requests a Student-t approximation with those degrees of freedom.
#'
#'   Monte Carlo: draw \eqn{\theta^{(b)} \sim \mathrm{MVN}(\hat{\theta},
#'   \hat{\Sigma})} for `b = 1, ..., n_sim`, compute \eqn{g^{(b)} =
#'   g(\theta^{(b)})}, and report the sample mean, sample SD, and the
#'   empirical quantiles of the simulated distribution.
#'
#'   Bootstrap: for each replicate `b`, draw `n` rows of `original_data`
#'   with replacement, refit the model via `refit(data_b)`, and compute
#'   \eqn{g^{(b)} = g(\hat{\theta}^{(b)})}. Percentile intervals are
#'   reported. Refit failures are dropped.
#'
#' @section Limitations:
#'
#'   - The delta method can be severely biased for highly nonlinear `func`
#'     or near parameter boundaries. Always compare against Monte Carlo if
#'     in doubt.
#'   - Monte Carlo assumes the parameter estimate is approximately Gaussian.
#'     For small samples or constrained parameters this is optimistic.
#'   - The bootstrap requires the original data and a working refit
#'     function. If `refit` fails on a substantial fraction of resamples,
#'     the resulting interval is biased. Refit-failure rate is not currently
#'     reported.
#'   - The covariance matrix `cov` is taken as-is; it is not symmetrised
#'     or projected to the nearest positive-definite matrix.
#'
#' @references
#'   - Oehlert, G. W. (1992). A Note on the Delta Method. The American
#'     Statistician, 46, 27-29.
#'     \doi{10.1080/00031305.1992.10475842}.
#'
#'   - Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
#'     Bootstrap. Chapman & Hall.
#'     \doi{10.1007/978-1-4899-4541-9}.
#'
#' @examples
#' set.seed(1)
#' # Fit a simple linear model y = a + b*x and propagate to g(a,b) = -a/b
#' x <- 1:20
#' y <- 2 + 0.5 * x + rnorm(x, sd = 1)
#' fit <- lm(y ~ x)
#' params <- coef(fit)
#' covmat <- vcov(fit)
#' g <- function(p) -p[1] / p[2]  # x-intercept
#'
#' # Delta method
#' r_delta <- parameter_uncertainty(params, covmat, g, method = "delta")
#' print(r_delta)
#'
#' # Monte Carlo
#' r_mc <- parameter_uncertainty(params, covmat, g, method = "monte_carlo",
#'                               n_sim = 1000, seed = 1)
#' print(r_mc)
#'
#' @export
parameter_uncertainty <- function(params, cov, func,
                                    method = c("delta", "monte_carlo",
                                               "bootstrap"),
                                    n_sim = 1000L, conf_level = 0.95,
                                    seed = NULL, original_data = NULL,
                                    df = Inf) {
  call <- match.call()
  method <- match.arg(method)
  .assert_numeric_vector(params, "params", min_length = 1L)
  if (!is.matrix(cov) || !is.numeric(cov)) {
    stop("'cov' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(cov) != length(params) || ncol(cov) != length(params)) {
    stop("'cov' must be a square matrix of size length(params).",
         call. = FALSE)
  }
  if (any(!is.finite(cov)) || !isSymmetric(cov) ||
      min(eigen(cov, symmetric = TRUE, only.values = TRUE)$values) < -1e-10) {
    stop("'cov' must be finite, symmetric and positive semidefinite.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be positive or Inf.", call. = FALSE)
  }
  if (!is.function(func)) {
    stop("'func' must be a function.", call. = FALSE)
  }
  .assert_positive_integer(n_sim, "n_sim")
  .assert_probability(conf_level, "conf_level")

  if (method == "bootstrap" && is.null(original_data)) {
    warning("method='bootstrap' requires original_data; falling back to delta method.",
            call. = FALSE)
    method <- "delta"
  }

  if (method == "delta") {
    est <- func(params)
    if (!is.numeric(est) || length(est) == 0L) {
      stop("'func' must return a numeric scalar or vector.", call. = FALSE)
    }
    est <- as.numeric(est)
    J <- .num_jacobian(params, func)
    # Cov(g) = J Cov(theta) J'
    var_g <- J %*% cov %*% t(J)
    se <- sqrt(pmax(diag(var_g), 0))
    if (is.infinite(df)) {
      tcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
    } else {
      tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = max(1L, df))
    }
    ci_lo <- est - tcrit * se
    ci_hi <- est + tcrit * se
    out <- list(estimate = est,
                se = as.numeric(se),
                ci_low = as.numeric(ci_lo),
                ci_high = as.numeric(ci_hi),
                method = "delta",
                jacobian = J,
                sim_samples = NULL,
                conf_level = conf_level,
                call = call)
  } else if (method == "monte_carlo") {
    old <- .set_seed(seed)
    on.exit(.restore_seed(old), add = TRUE)
    samples <- .mvn_sample(n_sim, params, cov)
    g_list <- apply(samples, 1L, function(p) as.numeric(func(p)))
    if (!is.matrix(g_list)) {
      g_list <- matrix(g_list, nrow = 1L)
    }
    g_samples <- t(g_list)
    est <- colMeans(g_samples, na.rm = TRUE)
    se <- apply(g_samples, 2L, .safe_sd)
    q_lo <- apply(g_samples, 2L, stats::quantile,
                  probs = (1 - conf_level) / 2, na.rm = TRUE)
    q_hi <- apply(g_samples, 2L, stats::quantile,
                  probs = 1 - (1 - conf_level) / 2, na.rm = TRUE)
    out <- list(estimate = as.numeric(est),
                se = as.numeric(se),
                ci_low = as.numeric(q_lo),
                ci_high = as.numeric(q_hi),
                method = "monte_carlo",
                jacobian = NULL,
                sim_samples = g_samples,
                conf_level = conf_level,
                call = call)
  } else {
    # Bootstrap
    if (!is.list(original_data) ||
        is.null(original_data$data) || is.null(original_data$refit)) {
      stop("For bootstrap, original_data must be a list with elements 'data' and 'refit'.",
           call. = FALSE)
    }
    data <- original_data$data
    refit <- original_data$refit
    if (!is.data.frame(data)) {
      stop("original_data$data must be a data frame.", call. = FALSE)
    }
    if (!is.function(refit)) {
      stop("original_data$refit must be a function.", call. = FALSE)
    }
    old <- .set_seed(seed)
    on.exit(.restore_seed(old), add = TRUE)
    n <- nrow(data)
    g0 <- tryCatch(as.numeric(func(params)), error = function(e) NULL)
    if (is.null(g0)) stop("'func' failed at the supplied params.", call. = FALSE)
    n_g <- length(g0)
    g_samples <- matrix(NA_real_, nrow = n_sim, ncol = n_g)
    n_failed <- 0L
    for (b in seq_len(n_sim)) {
      idx <- sample.int(n, replace = TRUE)
      boot_par <- tryCatch(refit(data[idx, , drop = FALSE]),
                          error = function(e) NULL)
      if (is.null(boot_par)) {
        n_failed <- n_failed + 1L
        next
      }
      gb <- tryCatch(as.numeric(func(boot_par)), error = function(e) NULL)
      if (is.null(gb) || length(gb) != n_g) {
        n_failed <- n_failed + 1L
        next
      }
      g_samples[b, ] <- gb
    }
    if (n_failed > 0L) {
      warning(sprintf("%d/%d bootstrap refits failed; results use the remainder.",
                      n_failed, n_sim), call. = FALSE)
    }
    est <- colMeans(g_samples, na.rm = TRUE)
    se <- apply(g_samples, 2L, .safe_sd)
    q_lo <- apply(g_samples, 2L, stats::quantile,
                  probs = (1 - conf_level) / 2, na.rm = TRUE)
    q_hi <- apply(g_samples, 2L, stats::quantile,
                  probs = 1 - (1 - conf_level) / 2, na.rm = TRUE)
    out <- list(estimate = as.numeric(est),
                se = as.numeric(se),
                ci_low = as.numeric(q_lo),
                ci_high = as.numeric(q_hi),
                method = "bootstrap",
                jacobian = NULL,
                sim_samples = g_samples,
                n_failed = n_failed,
                conf_level = conf_level,
                call = call)
  }

  out$tidy <- data.frame(
    component = seq_along(out$estimate),
    estimate = out$estimate,
    se = out$se,
    ci_low = out$ci_low,
    ci_high = out$ci_high,
    stringsAsFactors = FALSE
  )

  .crossstatr_result(out, "crossstatr_param_uncert")
}

# Sample from a multivariate normal. Uses eigendecomposition so that
# non-PD matrices (with small negative eigenvalues) still produce a result.
# @noRd
.mvn_sample <- function(n, mu, cov) {
  p <- length(mu)
  e <- eigen(cov, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), p)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  out <- z %*% t(L)
  out <- sweep(out, 2L, mu, "+")
  out
}

#' @export
print.crossstatr_param_uncert <- function(x, ...) {
  cat(sprintf("<crossStatR::parameter_uncertainty> method = %s\n", x$method))
  cat(sprintf("estimate = %s\n",
              paste(sprintf("%.4g", x$estimate), collapse = ", ")))
  cat(sprintf("se       = %s\n",
              paste(sprintf("%.4g", x$se), collapse = ", ")))
  cat(sprintf("%g%% CI:  [%s]\n",
              100 * x$conf_level,
              paste(sprintf("[%.4g, %.4g]", x$ci_low, x$ci_high),
                    collapse = ", ")))
  if (x$method == "delta") {
    cat("(Delta method; linearised at the supplied params.)\n")
  }
  invisible(x)
}
