#' Curve similarity between two curves on possibly different x grids
#'
#' Quantifies the similarity of two curves (x, y) sampled on possibly
#' different x grids. When `align = TRUE`, both curves are linearly
#' interpolated onto a common grid (the overlap of the two x ranges) before
#' a distance or correlation metric is computed. If `align = FALSE`, the
#' curves must already share the same x coordinates.
#'
#' Several distance definitions are offered: the root-mean-square (L2), the
#' mean absolute (L1), the maximum absolute (max), the dynamic time warping
#' (DTW) cost, the Pearson correlation, and the Frobenius norm. For the
#' correlation method the package returns a Fisher-z confidence interval; for the
#' distance-based methods the package returns a percentile bootstrap interval obtained
#' by resampling the paired residuals.
#'
#' Equivalence testing is offered as a one-sided statement. If the user
#' supplies an `equivalence_threshold` (delta) the function reports whether
#' the observed distance falls below it. The package deliberately do not claim
#' equivalence when the observed distance exceeds delta. A small observed
#' distance supports equivalence; a large observed distance is silent on
#' equivalence until a formal two-one-sided-test (TOST) or interval-based
#' equivalence procedure is performed. Absence of evidence is not evidence
#' of absence.
#'
#' @param x1,y1 numeric vectors giving the first curve. `x1` should be
#'   suitable for `stats::approx` (typically increasing).
#' @param x2,y2 numeric vectors giving the second curve.
#' @param method character; one of `"l2"`, `"l1"`, `"max"`, `"dtw"`,
#'   `"correlation"`, `"frobenius"`.
#' @param align logical; if `TRUE` (default), interpolate both curves onto a
#'   common grid. If `FALSE`, `x1` and `x2` must be identical.
#' @param conf_level confidence level for the bootstrap or Fisher-z interval.
#' @param n_boot integer; number of bootstrap replicates for the distance CI.
#' @param seed optional integer seed for the bootstrap.
#' @param equivalence_threshold optional non-negative scalar; if supplied,
#'   the function reports whether `distance < equivalence_threshold`.
#'
#' @return A `crossstatr_curve_similarity` object with elements:
#'   \item{distance}{named numeric scalar: the chosen distance / correlation.}
#'   \item{ci_low, ci_high}{numeric bounds of the confidence interval.}
#'   \item{method}{method label.}
#'   \item{n_points}{number of grid points used.}
#'   \item{equivalence_threshold}{the threshold if supplied, else `NULL`.}
#'   \item{equivalent}{logical; `TRUE` if `distance < equivalence_threshold`,
#'     `NA` if no threshold was supplied.}
#'   \item{conf_level}{confidence level used.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   L2 distance (root-mean-square):
#'   \deqn{d_{L2} = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(y_{1i} - y_{2i})^2}.}
#'
#'   L1 distance (mean absolute):
#'   \deqn{d_{L1} = \frac{1}{n}\sum_{i=1}^{n}|y_{1i} - y_{2i}|.}
#'
#'   Maximum absolute deviation:
#'   \deqn{d_{\max} = \max_{i}|y_{1i} - y_{2i}|.}
#'
#'   Frobenius norm (unnormalised):
#'   \deqn{d_{F} = \sqrt{\sum_{i=1}^{n}(y_{1i} - y_{2i})^2}.}
#'
#'   DTW cost (Sakoe and Chiba, 1978) uses the recursion
#'   \deqn{D[i, j] = |y_{1i} - y_{2j}| + \min\{D[i-1, j], D[i, j-1], D[i-1, j-1]\},}
#'   with `D[1, 1] = |y_{11} - y_{21}|`. The reported cost is `D[n, m]`.
#'   No window constraint is imposed, so the algorithm is O(nm) in both
#'   time and memory.
#'
#'   Pearson correlation: \eqn{r = \mathrm{cor}(y_1, y_2)} with Fisher-z
#'   confidence interval on \eqn{z = \mathrm{atanh}(r)} of half-width
#'   \eqn{t_{1-\alpha/2,\,n-3}/\sqrt{n-3}} (Fisher, 1921).
#'
#'   Bootstrap CI: the paired residuals `e_i = y_{1i} - y_{2i}` are resampled
#'   with replacement `n_boot` times and the distance is recomputed on each
#'   bootstrap sample. The percentile interval at level `conf_level` is
#'   reported. DTW is not bootstrapped because time-ordering is integral to
#'   the metric.
#'
#' @section Limitations:
#'
#'   - Distance intervals use independent residual resampling rather than
#'     paired resampling of repeated experimental curves. Dependence within
#'     and between curves is not preserved by that bootstrap.
#'
#'   - Linear interpolation onto a common grid assumes the true curve is
#'     piecewise linear between sampled x points. Coarse sampling reduces
#'     the fidelity of the interpolated values.
#'   - The DTW implementation is O(nm) in time and memory; with large n it
#'     can be slow and memory-hungry.
#'   - Bootstrap CIs are percentile-based and require a reasonable sample
#'     size. For `n < 10` they can be unstable.
#'   - Failure to detect a difference is NOT proof of equivalence. The
#'     equivalence flag is one-directional: a small observed distance
#'     supports equivalence; a large observed distance is silent on
#'     equivalence until a formal TOST or interval-inclusion procedure is
#'     applied. See Wellek (2010) for equivalence testing methodology.
#'
#' @references
#'   - Sakoe, H. and Chiba, S. (1978). Dynamic programming algorithm
#'     optimization for spoken word recognition. IEEE Transactions on
#'     Acoustics, Speech, and Signal Processing, 26, 43-49.
#'     \doi{10.1109/TASSP.1978.1163055}.
#'
#'   - Fisher, R. A. (1921). On the probable error of a coefficient of
#'     correlation deduced from a small sample. Metron, 1, 3-32.
#'     \url{https://digital.library.adelaide.edu.au/items/002ad8fb-c23c-407b-8a89-f036d8da6030}.
#'
#'   - Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
#'     Bootstrap. Chapman & Hall.
#'     \doi{10.1007/978-1-4899-4541-9}.
#'
#'   - Wellek, S. (2010). Testing Statistical Hypotheses of Equivalence and
#'     Noninferiority. 2nd ed. Chapman & Hall/CRC.
#'     \doi{10.1201/EBK1439808184}.
#'
#' @examples
#' set.seed(1)
#' t1 <- seq(0, 5, length.out = 50)
#' t2 <- seq(0, 5, length.out = 80)
#' y1 <- sin(t1) + rnorm(t1, sd = 0.05)
#' y2 <- sin(t2 + 0.05) + rnorm(t2, sd = 0.05)
#' result <- curve_similarity(t1, y1, t2, y2, method = "l2")
#' print(result)
#'
#' # DTW handles unequal sampling and small shifts
#' curve_similarity(t1, y1, t2, y2, method = "dtw")
#'
#' # Correlation with Fisher-z CI
#' curve_similarity(t1, y1, t2, y2, method = "correlation")
#'
#' # Equivalence statement against a tolerance
#' curve_similarity(t1, y1, t2, y2, method = "l2",
#'                  equivalence_threshold = 0.2)
#'
#' @export
curve_similarity <- function(x1, y1, x2, y2,
                              method = c("l2", "l1", "max", "dtw",
                                         "correlation", "frobenius"),
                              align = TRUE, conf_level = 0.95,
                              n_boot = 500L, seed = NULL,
                              equivalence_threshold = NULL) {
  call <- match.call()
  method <- match.arg(method)
  .assert_numeric_vector(x1, "x1", min_length = 2L)
  .assert_numeric_vector(y1, "y1", min_length = 2L)
  .assert_numeric_vector(x2, "x2", min_length = 2L)
  .assert_numeric_vector(y2, "y2", min_length = 2L)
  if (length(x1) != length(y1)) {
    stop("'x1' and 'y1' must have equal length.", call. = FALSE)
  }
  if (length(x2) != length(y2)) {
    stop("'x2' and 'y2' must have equal length.", call. = FALSE)
  }
  .assert_flag(align, "align")
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(n_boot, "n_boot")
  if (!is.null(equivalence_threshold)) {
    if (!is.numeric(equivalence_threshold) ||
        length(equivalence_threshold) != 1L ||
        !is.finite(equivalence_threshold) ||
        equivalence_threshold < 0) {
      stop("'equivalence_threshold' must be a non-negative finite scalar.",
           call. = FALSE)
    }
  }

  if (align) {
    lo <- max(min(x1), min(x2))
    hi <- min(max(x1), max(x2))
    if (lo >= hi) {
      stop("Curves do not have an overlapping x range.", call. = FALSE)
    }
    n_target <- min(length(x1), length(x2))
    target_x <- seq(lo, hi, length.out = n_target)
    yy1 <- .interp_to_grid(x1, y1, target_x)
    yy2 <- .interp_to_grid(x2, y2, target_x)
    grid_x <- target_x
  } else {
    if (length(x1) != length(x2) ||
        !isTRUE(all.equal(as.numeric(x1), as.numeric(x2)))) {
      stop("With align=FALSE, 'x1' and 'x2' must be identical.",
           call. = FALSE)
    }
    yy1 <- y1; yy2 <- y2
    grid_x <- x1
  }
  n <- length(yy1)

  d <- .curve_distance(yy1, yy2, method)

  if (method == "correlation") {
    if (is.finite(d) && abs(d) < 1 && n > 3L) {
      z <- atanh(d)
      se_z <- 1 / sqrt(n - 3)
      tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = n - 3)
      ci_lo <- tanh(z - tcrit * se_z)
      ci_hi <- tanh(z + tcrit * se_z)
    } else {
      ci_lo <- ci_hi <- d
    }
  } else if (method == "dtw") {
    ci_lo <- ci_hi <- NA_real_
  } else {
    old <- .set_seed(seed)
    on.exit(.restore_seed(old), add = TRUE)
    resids <- yy1 - yy2
    boot_vals <- replicate(n_boot, {
      e_star <- sample(resids, replace = TRUE)
      .resid_distance(e_star, method, n)
    })
    ci_lo <- stats::quantile(boot_vals, (1 - conf_level) / 2, na.rm = TRUE)
    ci_hi <- stats::quantile(boot_vals, 1 - (1 - conf_level) / 2, na.rm = TRUE)
    names(ci_lo) <- names(ci_hi) <- NULL
  }

  equivalent <- if (!is.null(equivalence_threshold)) d < equivalence_threshold else NA

  tidy_df <- data.frame(
    method = method,
    distance = d,
    ci_low = unname(ci_lo),
    ci_high = unname(ci_hi),
    n_points = n,
    stringsAsFactors = FALSE
  )

  .crossstatr_result(
    list(distance = d,
         ci_low = unname(ci_lo),
         ci_high = unname(ci_hi),
         method = method,
         n_points = n,
         aligned = align,
         grid = grid_x,
         equivalence_threshold = equivalence_threshold,
         equivalent = equivalent,
         conf_level = conf_level,
         tidy = tidy_df,
         call = call),
    "crossstatr_curve_similarity"
  )
}

# Compute the chosen distance / correlation for two equal-length y vectors.
# @noRd
.curve_distance <- function(y1, y2, method) {
  switch(method,
    "l2" = sqrt(mean((y1 - y2)^2)),
    "l1" = mean(abs(y1 - y2)),
    "max" = max(abs(y1 - y2)),
    "frobenius" = sqrt(sum((y1 - y2)^2)),
    "correlation" = stats::cor(y1, y2, method = "pearson"),
    "dtw" = .dtw_distance(y1, y2),
    stop("Unknown method.", call. = FALSE)
  )
}

# Compute the distance directly on a residual vector e = y1 - y2.
# @noRd
.resid_distance <- function(e, method, n) {
  switch(method,
    "l2" = sqrt(mean(e^2)),
    "l1" = mean(abs(e)),
    "max" = max(abs(e)),
    "frobenius" = sqrt(sum(e^2)),
    stop("Bootstrap not supported for this method.", call. = FALSE)
  )
}

# Basic DTW with no window constraint. O(nm) memory.
# Returns the accumulated cost at D[n, m].
# @noRd
.dtw_distance <- function(y1, y2) {
  n <- length(y1); m <- length(y2)
  D <- matrix(0, nrow = n, ncol = m)
  D[1L, 1L] <- abs(y1[1L] - y2[1L])
  if (n > 1L) {
    for (i in 2:n) D[i, 1L] <- abs(y1[i] - y2[1L]) + D[i - 1L, 1L]
  }
  if (m > 1L) {
    for (j in 2:m) D[1L, j] <- abs(y1[1L] - y2[j]) + D[1L, j - 1L]
  }
  if (n > 1L && m > 1L) {
    for (i in 2:n) {
      for (j in 2:m) {
        D[i, j] <- abs(y1[i] - y2[j]) +
          min(D[i - 1L, j], D[i, j - 1L], D[i - 1L, j - 1L])
      }
    }
  }
  D[n, m]
}

#' @export
print.crossstatr_curve_similarity <- function(x, ...) {
  cat(sprintf("<crossStatR::curve_similarity> method = %s\n", x$method))
  cat(sprintf("distance = %.4g  (n = %d points, aligned = %s)\n",
              x$distance, x$n_points, x$aligned))
  if (is.finite(x$ci_low) && is.finite(x$ci_high)) {
    cat(sprintf("%g%% CI: [%.4g, %.4g]\n",
                100 * x$conf_level, x$ci_low, x$ci_high))
  } else {
    cat("CI: not available for this method.\n")
  }
  if (!is.null(x$equivalence_threshold)) {
    cat(sprintf("Equivalence threshold = %.4g; equivalent = %s\n",
                x$equivalence_threshold,
                if (isTRUE(x$equivalent)) "TRUE" else "FALSE"))
    cat("  (failure to exceed the threshold is not proof of equivalence;\n")
    cat("   a formal TOST or interval-inclusion test is required.)\n")
  }
  invisible(x)
}
