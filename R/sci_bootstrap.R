#' Bootstrap resampling with multiple CI types
#'
#' Nonparametric, parametric, cluster, and stratified bootstrap. CIs are
#' returned as percentile, normal-approximation, or BCa intervals. The
#' resampling is implemented directly so that the user can pass any
#' statistic function without depending on the `boot` package's calling
#' conventions, although the package follows Efron & Tibshirani (1993) and Efron
#' (1987) for the BCa formula.
#'
#' The `parallel` argument exists for API symmetry with `boot::boot` but is a
#' no-op here. The documentation does not recommend using it in this function; if the analyst needs
#' parallelism wrap the call in `parallel::mclapply` over `n_replicates`
#' batches yourself.
#'
#' @param data a numeric vector or data.frame. Passed to `statistic` after
#'   resampling.
#' @param statistic a function that takes a resampled dataset of the same
#'   shape as `data` and returns a numeric scalar. Functions returning longer
#'   vectors are not supported; wrap them to extract a single statistic.
#' @param n_replicates integer; number of bootstrap replicates (default 1000).
#' @param type `"nonparametric"` (default), `"parametric"`, `"cluster"`, or
#'   `"stratified"`.
#' @param cluster optional cluster IDs for `type = "cluster"`. Must be the
#'   same length as `data` (or `nrow(data)` if a data.frame).
#' @param strata optional strata for `type = "stratified"`. Same length
#'   requirement as `cluster`.
#' @param ci_type `"perc"` (default), `"bca"`, or `"normal"`.
#' @param conf_level confidence level for the CI (default 0.95).
#' @param seed optional integer seed for reproducibility.
#' @param parallel logical; accepted but ignored. Included only so that the
#'   signature mirrors `boot::boot`.
#' @param parametric_sampler optional function of a single argument `n`
#'   returning a synthetic dataset (vector or data.frame of `n` rows) for
#'   `type = "parametric"`. Required when `type = "parametric"`.
#'
#' @return A `crossstatr_bootstrap` object with elements:
#'   \item{statistic}{the statistic computed on the original data.}
#'   \item{replicates}{numeric vector of bootstrap replicates.}
#'   \item{ci}{named numeric length-2 confidence interval.}
#'   \item{n_replicates}{number of replicates actually used.}
#'   \item{type}{resampling type.}
#'   \item{ci_type}{CI type.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Nonparametric bootstrap: draw a resample of size `n` (with replacement)
#'   from the empirical distribution, recompute the statistic, repeat.
#'
#'   Parametric bootstrap: simulate a new dataset from a model fitted to the
#'   original data, supplied via `parametric_sampler`. The sampler is called
#'   `n_replicates` times.
#'
#'   Cluster bootstrap: resample whole clusters (with replacement) and stack
#'   the observations within the sampled clusters back into a dataset of
#'   size `n` (Davison and Hinkley, 1997).
#'
#'   Stratified bootstrap: resample (with replacement) independently within
#'   each stratum, preserving the per-stratum sample sizes.
#'
#'   Percentile CI: empirical \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles of
#'   the bootstrap distribution.
#'
#'   Normal CI: \eqn{\hat{\theta} \pm z_{1-\alpha/2}\,\hat{\sigma}_{boot}},
#'   where \eqn{\hat{\sigma}_{boot}} is the standard deviation of the
#'   replicates.
#'
#'   BCa (Efron, 1987): the interval is the percentile interval with
#'   adjusted quantiles. Let \eqn{\hat{\theta}} be the original estimate and
#'   \eqn{\hat{\theta}^*_b} the bootstrap replicates. The bias-correction is
#'   \deqn{z_0 = \Phi^{-1}\!\Bigl(\#\{\hat{\theta}^*_b < \hat{\theta}\}/B\Bigr).}
#'   The acceleration is obtained by jackknife:
#'   \deqn{a = \frac{\sum_{i}(\bar{\theta}_{(\cdot)} - \theta_{(i)})^3}
#'                {6\,\bigl[\sum_{i}(\theta_{(i)} - \bar{\theta}_{(\cdot)})^2\bigr]^{3/2}},}
#'   where \eqn{\theta_{(i)}} is the leave-one-out statistic. The adjusted
#'   quantile levels are
#'   \deqn{\alpha_1 = \Phi\!\left(z_0 + \frac{z_0 + z_{\alpha/2}}
#'                              {1 - a(z_0 + z_{\alpha/2})}\right), \qquad
#'         \alpha_2 = \Phi\!\left(z_0 + \frac{z_0 + z_{1-\alpha/2}}
#'                              {1 - a(z_0 + z_{1-\alpha/2})}\right).}
#'
#' @section Limitations:
#'
#'   - BCa is implemented for `type = "nonparametric"` only. For other types,
#'     the function falls back to the percentile interval and emits a
#'     warning. This is because the jackknife acceleration requires
#'     leave-one-out resampling of the original data, which is not well
#'     defined for the parametric or cluster cases without further
#'     specification.
#'   - The statistic must return a finite numeric scalar. Statistics
#'     returning vectors, lists, or non-numeric values will error.
#'   - For very small `n` the bootstrap distribution can be discrete and the
#'     percentile interval has poor coverage; consider `ci_type = "bca"` and
#'     inspecting the replicates directly.
#'   - The `parallel` argument is a no-op. The package does not run any background
#'     processes.
#'
#' @references
#'   - Efron, B. (1979). Bootstrap Methods: Another Look at the Jackknife.
#'     The Annals of Statistics, 7, 1-26.
#'     \doi{10.1214/aos/1176344552}.
#'
#'   - Efron, B. (1987). Better Bootstrap Confidence Intervals. Journal of
#'     the American Statistical Association, 82, 171-185.
#'     \doi{10.1080/01621459.1987.10478410}.
#'
#'   - Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
#'     Bootstrap. Chapman & Hall.
#'     \doi{10.1007/978-1-4899-4541-9}.
#'
#'   - Davison, A. C. and Hinkley, D. V. (1997). Bootstrap Methods and
#'     their Application. Cambridge University Press.
#'     \doi{10.1017/CBO9780511802843}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(30, mean = 5)
#' mean_stat <- function(d) mean(d)
#' result <- sci_bootstrap(x, mean_stat, n_replicates = 500, seed = 42)
#' print(result)
#'
#' # data.frame input
#' df <- data.frame(y = rnorm(40), g = rep(c("A","B"), each = 20))
#' stat <- function(d) diff(tapply(d$y, d$g, mean))
#' sci_bootstrap(df, stat, type = "stratified", strata = df$g,
#'               ci_type = "bca", n_replicates = 300, seed = 7)
#'
#' @export
sci_bootstrap <- function(data, statistic, n_replicates = 1000L,
                           type = c("nonparametric", "parametric",
                                    "cluster", "stratified"),
                           cluster = NULL, strata = NULL,
                           ci_type = c("perc", "bca", "normal"),
                           conf_level = 0.95, seed = NULL, parallel = FALSE,
                           parametric_sampler = NULL) {
  call <- match.call()
  type <- match.arg(type)
  ci_type <- match.arg(ci_type)
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(n_replicates, "n_replicates")
  .assert_flag(parallel, "parallel")
  if (!is.function(statistic)) {
    stop("'statistic' must be a function.", call. = FALSE)
  }

  is_df <- is.data.frame(data)
  n <- if (is_df) nrow(data) else length(data)
  if (n < 2L) stop("'data' must have at least 2 observations.", call. = FALSE)

  if (type == "parametric" && is.null(parametric_sampler)) {
    stop("'parametric_sampler' is required when type = 'parametric'.",
         call. = FALSE)
  }
  if (type == "cluster" && is.null(cluster)) {
    stop("'cluster' is required when type = 'cluster'.", call. = FALSE)
  }
  if (type == "stratified" && is.null(strata)) {
    stop("'strata' is required when type = 'stratified'.", call. = FALSE)
  }
  if (!is.null(cluster) && length(cluster) != n) {
    stop("'cluster' must have length n (or nrow(data)).", call. = FALSE)
  }
  if (!is.null(strata) && length(strata) != n) {
    stop("'strata' must have length n (or nrow(data)).", call. = FALSE)
  }

  original_stat <- .safe_statistic(statistic, data)

  # Seed handling
  old_seed <- .set_seed(seed)
  on.exit(.restore_seed(old_seed), add = TRUE)

  repls <- numeric(n_replicates)
  for (b in seq_len(n_replicates)) {
    resample <- .boot_resample(data, n, type, cluster, strata,
                              parametric_sampler)
    repls[b] <- .safe_statistic(statistic, resample)
  }

  # CI
  ci <- .boot_ci_full(repls, original_stat, ci_type, conf_level,
                       data, statistic, type)

  out <- .crossstatr_result(
    list(statistic = original_stat,
         replicates = repls,
         ci = ci,
         n_replicates = as.integer(n_replicates),
         type = type,
         ci_type = ci_type,
         conf_level = conf_level,
         call = call),
    "crossstatr_bootstrap"
  )
  out
}

# Evaluate the statistic and assert a finite numeric scalar was returned.
# @noRd
.safe_statistic <- function(statistic, data) {
  val <- tryCatch(statistic(data), error = function(e) {
    stop(sprintf("'statistic' failed: %s", conditionMessage(e)),
         call. = FALSE)
  })
  if (!is.numeric(val) || length(val) != 1L) {
    stop("'statistic' must return a numeric scalar (length-1 numeric).",
         call. = FALSE)
  }
  as.numeric(val)
}

# Generate a single bootstrap resample of the same shape as `data`.
# @noRd
.boot_resample <- function(data, n, type, cluster, strata,
                           parametric_sampler) {
  is_df <- is.data.frame(data)
  if (type == "nonparametric") {
    idx <- sample.int(n, size = n, replace = TRUE)
    if (is_df) data[idx, , drop = FALSE] else data[idx]
  } else if (type == "parametric") {
    parametric_sampler(n)
  } else if (type == "cluster") {
    idx <- .resample_clusters(cluster)
    if (length(idx) == 0L) idx <- sample.int(n, size = n, replace = TRUE)
    if (is_df) data[idx, , drop = FALSE] else data[idx]
  } else { # stratified
    if (is_df) {
      out <- do.call(rbind, lapply(split(seq_len(n), strata), function(ii) {
        data[sample(ii, size = length(ii), replace = TRUE), , drop = FALSE]
      }))
      out
    } else {
      out <- unlist(lapply(split(seq_len(n), strata), function(ii) {
        data[sample(ii, size = length(ii), replace = TRUE)]
      }), use.names = FALSE)
      out
    }
  }
}

# Compute a confidence interval on the bootstrap replicates.
# @noRd
.boot_ci_full <- function(repls, original_stat, ci_type, conf_level,
                          data, statistic, type) {
  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha / 2)
  if (ci_type == "normal") {
    se <- stats::sd(repls, na.rm = TRUE)
    ci <- c(original_stat - zcrit * se, original_stat + zcrit * se)
    names(ci) <- paste0(c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2) * 100,
                       "%")
    return(ci)
  }
  if (ci_type == "perc") {
    ci <- as.numeric(stats::quantile(repls,
                                    probs = c(alpha / 2, 1 - alpha / 2),
                                    na.rm = TRUE))
    names(ci) <- paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
    return(ci)
  }
  # BCa: requires jackknife on the original data, so only sensible for
  # nonparametric.
  if (type != "nonparametric") {
    warning("BCa CI is only supported for type = 'nonparametric'; ",
            "falling back to percentile.", call. = FALSE)
    ci <- as.numeric(stats::quantile(repls,
                                    probs = c(alpha / 2, 1 - alpha / 2),
                                    na.rm = TRUE))
    names(ci) <- paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
    return(ci)
  }
  ci <- .bca_interval(repls, original_stat, data, statistic, conf_level)
  names(ci) <- paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
  ci
}

# BCa interval (Efron 1987). Returns a length-2 numeric.
# @noRd
.bca_interval <- function(repls, original_stat, data, statistic, conf_level) {
  alpha <- 1 - conf_level
  n <- if (is.data.frame(data)) nrow(data) else length(data)
  # Jackknife leave-one-out
  jack <- vapply(seq_len(n), function(i) {
    if (is.data.frame(data)) {
      d_i <- data[-i, , drop = FALSE]
    } else {
      d_i <- data[-i]
    }
    .safe_statistic(statistic, d_i)
  }, numeric(1L))
  jbar <- mean(jack, na.rm = TRUE)
  diffs <- jack - jbar
  num <- sum(diffs^3)
  denom <- 6 * (sum(diffs^2))^1.5
  a_hat <- if (denom == 0 || !is.finite(denom)) 0 else num / denom

  # Bias correction z0
  prop_less <- mean(repls < original_stat, na.rm = TRUE)
  # Guard against prop = 0 or 1 which give infinite z0
  prop_less <- min(max(prop_less, 1e-8), 1 - 1e-8)
  z0 <- stats::qnorm(prop_less)

  z_lo <- stats::qnorm(alpha / 2)
  z_hi <- stats::qnorm(1 - alpha / 2)
  adj1 <- stats::pnorm(z0 + (z0 + z_lo) / (1 - a_hat * (z0 + z_lo)))
  adj2 <- stats::pnorm(z0 + (z0 + z_hi) / (1 - a_hat * (z0 + z_hi)))
  as.numeric(stats::quantile(repls, probs = c(adj1, adj2), na.rm = TRUE))
}

#' @export
print.crossstatr_bootstrap <- function(x, ...) {
  cat(sprintf("<crossStatR::bootstrap> type = %s, ci_type = %s\n",
              x$type, x$ci_type))
  cat(sprintf("statistic (original) = %.4f\n", x$statistic))
  cat(sprintf("replicates: n = %d, mean = %.4f, sd = %.4f\n",
              length(x$replicates), mean(x$replicates),
              stats::sd(x$replicates)))
  cat(sprintf("%g%% CI: [%.4f, %.4f]\n", 100 * x$conf_level,
              x$ci[1L], x$ci[2L]))
  invisible(x)
}
