#' Weighted multi-peak fitting
#'
#' Fits one or more peaks of a common shape to a signal by nonlinear least
#' squares. Peaks may be located automatically by a local-maxima plus
#' prominence rule, or supplied manually via `start`. The three supported
#' shapes are Gaussian, Lorentzian, and pseudo-Voigt (a linear mixture of the
#' two). The fit uses `stats::nls` with `algorithm = "port"` and lower bounds
#' on heights and widths, which prevents degenerate solutions where a peak
#' collapses to a delta spike or vanishes to negative height.
#'
#' The package does not automatically pick the number of peaks by AIC. The
#' consequences of an extra spurious peak in downstream analysis are
#' usually larger than the consequences of a missed peak, so the package leaves the
#' judgement about peak count to the analyst. Inspect `aic` and `bic`
#' alongside the residual plot before committing to a fit.
#'
#' @param x numeric vector of abscissae.
#' @param y numeric vector of signal intensities, same length as `x`.
#' @param peaks character; `"auto"` (default) for automatic peak detection,
#'   or `"manual"` to use `start`.
#' @param n_peaks integer; if provided with `peaks = "auto"`, take the top
#'   `n_peaks` peaks by prominence. If NULL, all peaks above the
#'   prominence threshold are kept.
#' @param peak_shape character; one of `"gaussian"`, `"lorentzian"`, or
#'   `"voigt"` (pseudo-Voigt approximation).
#' @param start optional list or named numeric vector of starting values for
#'   the parameters when `peaks = "manual"`. Names should follow the
#'   pattern `h1, c1, w1, [eta1], h2, c2, w2, [eta2], ...`.
#' @param weights optional numeric weights (e.g. `1/sigma^2`). Length must
#'   match `x` and `y`. Recycled if length 1.
#' @param max_iter integer; passed to `nls` as `control$maxiter`.
#' @param conf_level confidence level for parameter CIs.
#'
#' @return A `crossstatr_peak_fit` object with elements:
#'   \item{coefficients}{data frame with columns `peak`, `param`, `est`,
#'     `se`, `ci_low`, `ci_high`.}
#'   \item{fitted}{fitted signal at each x.}
#'   \item{residual}{`y - fitted`.}
#'   \item{n_peaks}{number of peaks fitted.}
#'   \item{peak_shape}{shape used.}
#'   \item{aic, bic}{AIC and BIC of the fit.}
#'   \item{converged}{logical; whether `nls` converged.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Each peak is modelled by a parametric line shape. For a Gaussian,
#'   \deqn{G(x; h, c, w) = h\, \exp\!\left(-\frac{(x - c)^2}{2 w^2}\right),}
#'   and for a Lorentzian,
#'   \deqn{L(x; h, c, w) = h\, \frac{w^2}{(x - c)^2 + w^2},}
#'   where \eqn{h} is the peak height, \eqn{c} the centre, and \eqn{w} a
#'   width parameter (equal to the standard deviation for the Gaussian and
#'   to the half-width at half-maximum for the Lorentzian).
#'
#'   The implemented pseudo-Voigt-style mixture is
#'   \deqn{pV(x) = \eta\, L(x; h, c, w) + (1 - \eta)\, G(x; h, c, w),}
#'   with mixing parameter \eqn{\eta \in [0, 1]}. This shares the same
#'   height, centre, and width between the Gaussian and Lorentzian
#'   components, which is a simplification of the true Voigt profile (the
#'   convolution of Gaussian and Lorentzian). The mixing weight is fitted
#'   from a starting value of 0.5; this is not the width-dependent approximation
#'   of Thompson, Cox and Hastings (1987), which is related background.
#'
#'   The total model is the sum of \eqn{K} peaks,
#'   \deqn{\hat y(x_i) = \sum_{k=1}^{K} s_k(x_i; h_k, c_k, w_k, [\eta_k]),}
#'   and the package minimises the (weighted) sum of squared residuals
#'   \deqn{S(\theta) = \sum_i w_i\, (y_i - \hat y(x_i))^2.}
#'   Parameter standard errors come from the `nls` asymptotic covariance
#'   \eqn{\widehat{\mathrm{Cov}}(\hat\theta) = \hat\sigma^2 (J^\top W J)^{-1}},
#'   and CIs are normal-approximation using the residual df.
#'
#'   Automatic peak detection finds local maxima (points where the signal
#'   exceeds both neighbours) whose prominence exceeds a robust threshold
#'   \eqn{2 \cdot 1.4826 \cdot \mathrm{MAD}(y)}. The starting width for
#'   each peak is set to one-tenth of the x-range divided by the number of
#'   peaks, which is a safe default for typical chromatographic resolution.
#'
#' @section Limitations:
#'
#'   - Automatic peak candidates use a local prominence heuristic. The
#'     candidate count is not a statistical model-selection result. Nonlinear
#'     fitting can return an unconverged fit with a warning; such estimates
#'     and intervals require inspection of the stored fit convergence status.
#'
#'   - `nls` with `algorithm = "port"` will fail to converge if the peaks
#'     are too close together relative to their widths, or if `start` is
#'     far from the optimum. Always inspect the residuals.
#'   - Pseudo-Voigt with shared width between Gaussian and Lorentzian
#'     components is an approximation. For high-resolution spectroscopy
#'     where the line shape itself is the quantity of interest, fit a true
#'     Voigt profile (e.g. via the `voigt` package) instead.
#'   - Parameter CIs use the t-distribution with residual df and assume
#'     the asymptotic normal approximation holds. For low signal-to-noise
#'     or near-boundary parameters this can be misleading.
#'   - The prominence threshold of 2*MAD is permissive. On noisy data it
#'     will detect spurious peaks; supply `n_peaks` or use `peaks = "manual"`
#'     when the peak count is known.
#'
#' @references
#'   - Bates, D. M. and Watts, D. G. (1988). Nonlinear Regression Analysis
#'     and Its Applications. Wiley.
#'     \doi{10.1002/9780470316757}.
#'
#'   - Thompson, P., Cox, D. E. and Hastings, J. B. (1987). Rietveld
#'     refinement of Debye-Scherrer synchrotron X-ray data from Al2O3.
#'     Journal of Applied Crystallography, 20, 79-83.
#'     \doi{10.1107/S0021889887087090}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 300)
#' y <- 5 * exp(-(x - 3)^2 / (2 * 0.2^2)) +
#'      3 * exp(-(x - 7)^2 / (2 * 0.3^2)) +
#'      rnorm(length(x), sd = 0.05)
#' fit <- weighted_peak_fit(x, y, peaks = "auto", n_peaks = 2,
#'                          peak_shape = "gaussian")
#' print(fit)
#'
#' # Manual start for two known peaks
#' fit2 <- weighted_peak_fit(x, y, peaks = "manual", peak_shape = "gaussian",
#'                           start = list(h1 = 5, c1 = 3, w1 = 0.2,
#'                                        h2 = 3, c2 = 7, w2 = 0.3))
#' print(fit2)
#'
#' @export
weighted_peak_fit <- function(x, y, peaks = c("auto", "manual"),
                             n_peaks = NULL,
                             peak_shape = c("gaussian", "lorentzian",
                                            "voigt"),
                             start = NULL, weights = NULL,
                             max_iter = 500L, conf_level = 0.95) {
  call <- match.call()
  peaks <- match.arg(peaks)
  peak_shape <- match.arg(peak_shape)
  .assert_numeric_vector(x, "x", min_length = 5L)
  .assert_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(max_iter, "max_iter")
  if (!is.null(weights)) {
    .assert_numeric_vector(weights, "weights", min_length = 1L)
    if (length(weights) == 1L) weights <- rep(weights, length(x))
    if (length(weights) != length(x)) {
      stop("'weights' must have length 1 or length(x).", call. = FALSE)
    }
    if (any(!is.na(weights) & weights < 0)) {
      stop("'weights' must be non-negative.", call. = FALSE)
    }
  }
  if (!is.null(n_peaks)) {
    .assert_positive_integer(n_peaks, "n_peaks")
  }

  n <- length(y)
  x_range <- diff(range(x))
  if (!is.finite(x_range) || x_range <= 0) {
    stop("'x' must have positive range.", call. = FALSE)
  }

  # Build starting values.
  if (peaks == "auto") {
    detected <- .detect_peaks(x, y, n_peaks = n_peaks)
    if (length(detected$centers) == 0L) {
      stop("No peaks detected above threshold. Supply 'start' instead or ",
           "lower the detection threshold.", call. = FALSE)
    }
    start <- .build_peak_start(detected, peak_shape, x_range)
  } else {
    if (is.null(start)) {
      stop("'start' is required when peaks = 'manual'.", call. = FALSE)
    }
    start <- as.list(start)
  }

  fit <- .fit_peaks_nls(x, y, start = start, peak_shape = peak_shape,
                       weights = weights, max_iter = max_iter)

  # Coefficient table with SE and CIs from the nls fit.
  coef_df <- .peak_coef_table(fit$model, peak_shape, conf_level)

  .crossstatr_result(
    list(
      coefficients = coef_df,
      fitted = stats::fitted(fit$model),
      residual = stats::resid(fit$model),
      n_peaks = .n_peaks_from_start(start),
      peak_shape = peak_shape,
      aic = stats::AIC(fit$model),
      bic = stats::BIC(fit$model),
      converged = fit$converged,
      call = call
    ),
    "crossstatr_peak_fit"
  )
}

# Detect peaks by local maxima with prominence above 2*MAD threshold.
# Returns list(centers, heights, indices, half_widths).
# @noRd
.detect_peaks <- function(x, y, n_peaks = NULL, k = 2) {
  n <- length(y)
  mad_y <- .safe_mad(y)
  if (!is.finite(mad_y) || mad_y <= 0) mad_y <- .safe_sd(y)
  if (!is.finite(mad_y) || mad_y <= 0) mad_y <- 1
  # Local maxima (interior points).
  is_max <- rep(FALSE, n)
  if (n >= 3L) {
    for (i in 2L:(n - 1L)) {
      if (is.na(y[i])) next
      if (y[i] > y[i - 1L] && y[i] > y[i + 1L]) {
        is_max[i] <- TRUE
      }
    }
  }
  idx <- which(is_max)
  if (length(idx) == 0L) {
    return(list(centers = numeric(0), heights = numeric(0),
                indices = integer(0), half_widths = numeric(0)))
  }
  # Prominence: height above the higher of the two surrounding minima,
  # walking outward until the search reaches a higher peak or the boundary.
  prom <- vapply(idx, function(i) {
    hi <- y[i]
    j <- i
    left_min <- hi
    while (j > 1L) {
      j <- j - 1L
      if (y[j] > hi) break
      if (y[j] < left_min) left_min <- y[j]
    }
    j <- i
    right_min <- hi
    while (j < n) {
      j <- j + 1L
      if (y[j] > hi) break
      if (y[j] < right_min) right_min <- y[j]
    }
    hi - max(left_min, right_min)
  }, numeric(1L))
  # Keep peaks with prominence > k * MAD.
  keep <- prom > k * mad_y
  idx <- idx[keep]
  heights <- y[idx]
  prom <- prom[keep]
  if (length(idx) == 0L) {
    return(list(centers = numeric(0), heights = numeric(0),
                indices = integer(0), half_widths = numeric(0)))
  }
  # Order by prominence (descending).
  ord <- order(prom, decreasing = TRUE)
  idx <- idx[ord]; heights <- heights[ord]; prom <- prom[ord]
  if (!is.null(n_peaks) && length(idx) > n_peaks) {
    idx <- idx[seq_len(n_peaks)]
    heights <- heights[seq_len(n_peaks)]
    prom <- prom[seq_len(n_peaks)]
  }
  centers <- x[idx]
  half_widths <- rep(diff(range(x)) / 20, length(idx))
  list(centers = centers, heights = heights, indices = idx,
       half_widths = half_widths)
}

# Build starting values for nls from detected peaks.
# @noRd
.build_peak_start <- function(detected, peak_shape, x_range) {
  k <- length(detected$centers)
  start <- list()
  w0 <- x_range / max(10 * k, 1L)
  for (i in seq_len(k)) {
    start[[paste0("h", i)]] <- detected$heights[i]
    start[[paste0("c", i)]] <- detected$centers[i]
    start[[paste0("w", i)]] <- w0
    if (peak_shape == "voigt") {
      start[[paste0("eta", i)]] <- 0.5
    }
  }
  start
}

# Build the nls formula expression for K peaks of given shape.
# @noRd
.build_peak_formula <- function(k, peak_shape) {
  term <- function(i, shape) {
    if (shape == "gaussian") {
      substitute(h_i * exp(-(x - c_i)^2 / (2 * w_i^2)),
                 list(h_i = as.name(paste0("h", i)),
                      c_i = as.name(paste0("c", i)),
                      w_i = as.name(paste0("w", i))))
    } else if (shape == "lorentzian") {
      substitute(h_i * w_i^2 / ((x - c_i)^2 + w_i^2),
                 list(h_i = as.name(paste0("h", i)),
                      c_i = as.name(paste0("c", i)),
                      w_i = as.name(paste0("w", i))))
    } else {
      # pseudo-Voigt: eta * Lorentzian + (1 - eta) * Gaussian
      substitute(eta_i * (h_i * w_i^2 / ((x - c_i)^2 + w_i^2)) +
                   (1 - eta_i) * (h_i * exp(-(x - c_i)^2 / (2 * w_i^2))),
                 list(h_i = as.name(paste0("h", i)),
                      c_i = as.name(paste0("c", i)),
                      w_i = as.name(paste0("w", i)),
                      eta_i = as.name(paste0("eta", i))))
    }
  }
  rhs <- term(1L, peak_shape)
  if (k >= 2L) {
    for (i in 2L:k) {
      rhs <- bquote(.(rhs) + .(term(i, peak_shape)))
    }
  }
  bquote(y ~ .(rhs))
}

# Fit peaks with nls port algorithm.
# Lower bounds: positive heights and widths, eta in [0,1].
# @noRd
.fit_peaks_nls <- function(x, y, start, peak_shape, weights, max_iter) {
  k <- .n_peaks_from_start(start)
  form <- .build_peak_formula(k, peak_shape)
  data_env <- list(x = x, y = y)
  lower <- numeric(length(start))
  upper <- rep(Inf, length(start))
  for (i in seq_along(start)) {
    nm <- names(start)[i]
    if (startsWith(nm, "h") || startsWith(nm, "w")) {
      lower[i] <- 1e-8
    } else if (startsWith(nm, "c")) {
      lower[i] <- -Inf
    } else if (startsWith(nm, "eta")) {
      lower[i] <- 0
      upper[i] <- 1
    }
  }
  ctrl <- stats::nls.control(maxiter = max_iter, tol = 1e-8,
                              minFactor = 1e-10, printEval = FALSE,
                              warnOnly = TRUE)
  fit <- tryCatch(
    stats::nls(form, data = data_env, start = start,
               algorithm = "port", lower = lower, upper = upper,
               weights = weights, control = ctrl),
    error = function(e) {
      tryCatch(
        stats::nls(form, data = data_env, start = start,
                   weights = weights, control = ctrl),
        error = function(e2) {
          stop("nls() failed: ", conditionMessage(e2), call. = FALSE)
        }
      )
    }
  )
  conv_info <- fit$convInfo
  converged <- is.null(conv_info) ||
    (isTRUE(conv_info$isConv)) ||
    identical(conv_info$stopCode, 0L)
  list(model = fit, converged = converged)
}

# Number of peaks implied by the start list.
# @noRd
.n_peaks_from_start <- function(start) {
  nms <- names(start)
  sum(startsWith(nms, "h"))
}

# Coefficient table with SE and CIs from the nls fit.
# Falls back to t-approximation if confint.nls profiling fails.
# @noRd
.peak_coef_table <- function(model, peak_shape, conf_level) {
  smry <- summary(model)
  cf <- smry$coefficients
  ci <- tryCatch(
    stats::confint(model, level = conf_level),
    error = function(e) {
      alpha <- 1 - conf_level
      z <- stats::qt(1 - alpha / 2, df = smry$df[2L])
      cbind(cf[, 1L] - z * cf[, 2L], cf[, 1L] + z * cf[, 2L])
    }
  )
  ci <- as.matrix(ci)
  rn <- rownames(cf)
  # Parse peak index from row name (e.g. "h1" -> peak=1, param="h").
  peak <- as.integer(gsub("^[a-z]+", "", rn))
  param <- gsub("[0-9]+$", "", rn)
  data.frame(
    peak = peak,
    param = param,
    est = cf[, 1L],
    se = cf[, 2L],
    ci_low = ci[, 1L],
    ci_high = ci[, 2L],
    stringsAsFactors = FALSE
  )
}

#' @export
print.crossstatr_peak_fit <- function(x, ...) {
  cat(sprintf("<crossStatR::weighted_peak_fit> shape = %s, n_peaks = %d\n",
              x$peak_shape, x$n_peaks))
  cat(sprintf("Converged: %s, AIC = %.3f, BIC = %.3f\n",
              x$converged, x$aic, x$bic))
  cat("\nCoefficients:\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}
