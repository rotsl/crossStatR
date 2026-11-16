#' Adaptive smoothing for noisy signals
#'
#' Smooths a noisy signal with a choice of three backends: Savitzky-Golay
#' polynomial smoothing, smoothing spline, and lowess. For Savitzky-Golay,
#' the window size is chosen automatically by minimising leave-one-out
#' PRESS over a candidate grid, subject to a peak-preservation constraint
#' that vetoes windows which reduce the global maximum by more than 5%.
#'
#' Smoothing is a bias-variance tradeoff. Aggressive smoothing reduces
#' variance (noise) at the cost of bias: peak heights shrink, peak positions
#' can shift, and narrow peaks may disappear entirely. The package documents this
#' carefully because the bias propagates into downstream inference. A
#' smoothed peak height is not a peak height, it is a biased estimate of
#' one; the same applies to peak widths, peak areas, and any test statistic
#' computed on smoothed data. Whenever possible, fit a model directly to
#' the raw noisy data and let the likelihood handle the noise.
#'
#' @param x numeric vector of abscissae.
#' @param y numeric vector of signal intensities, same length as `x`.
#' @param method character; one of `"sgolay"` (default), `"spline"`, or
#'   `"lowess"`.
#' @param window integer; window size for `"sgolay"`. Must be odd. If NULL,
#'   chosen adaptively from the candidate set `c(5, 7, 9, 11, 15, 21)`.
#' @param poly_order integer; polynomial order for `"sgolay"`. Must be
#'   less than `window`.
#' @param spar numeric; smoothing parameter for `"spline"`, passed to
#'   `stats::smooth.spline`. If NULL, the spline chooses its own.
#' @param span numeric in (0, 1]; span for `"lowess"`.
#' @param preserve_peaks logical; if TRUE (default), the adaptive search
#'   for `sgolay` vetoes windows that reduce the global maximum of `y` by
#'   more than 5%.
#'
#' @return A `crossstatr_smoothing` object with elements:
#'   \item{x, y}{inputs.}
#'   \item{smoothed}{smoothed signal.}
#'   \item{method}{smoothing method used.}
#'   \item{window}{window size, if `method = "sgolay"`; `NA_integer_`
#'     otherwise.}
#'   \item{spar}{spline spar, if `method = "spline"`; `NA_real_`
#'     otherwise.}
#'   \item{span}{lowess span, if `method = "lowess"`; `NA_real_`
#'     otherwise.}
#'   \item{press}{cross-validated PRESS of the chosen smoothing
#'     parameter.}
#'   \item{peak_height_loss}{fraction of the global maximum lost to
#'     smoothing.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Savitzky-Golay smoothing (Savitzky & Golay, 1964; Steinier, Termonia &
#'   Deltour, 1972) fits a polynomial of order \eqn{m} in a sliding window
#'   of size \eqn{2h + 1} and replaces the central point with the value of
#'   the polynomial at the centre. Equivalently, the smoothed value at
#'   point \eqn{i} is
#'   \deqn{\hat y_i = \sum_{j=-h}^{h} c_j\, y_{i+j},}
#'   where the coefficient vector \eqn{c} is the first row of
#'   \eqn{(A^\top A)^{-1} A^\top} and \eqn{A} is the Vandermonde design
#'   matrix of the polynomial in the window. The package implements this directly
#'   by solving for the polynomial at each window and reading off the
#'   central fitted value, with reflection padding at the boundaries.
#'
#'   Smoothing spline: `stats::smooth.spline(x, y, spar = spar)` minimises
#'   the penalised residual sum of squares
#'   \deqn{\sum_i (y_i - f(x_i))^2 + \lambda \int [f''(t)]^2\, dt,}
#'   with `spar` controlling the tradeoff.
#'
#'   Lowess (Cleveland, 1979) performs locally weighted polynomial
#'   regression with tricube weights, controlled by `span`.
#'
#'   The cross-validated PRESS for a smoothing parameter \eqn{\theta} is
#'   \deqn{\mathrm{PRESS}(\theta) = \sum_i (y_i - \hat y^{(-i)}_{\theta}(x_i))^2,}
#'   where \eqn{\hat y^{(-i)}_{\theta}} is the smoothed signal computed
#'   without point \eqn{i}. For linear smoothers this can be computed
#'   exactly via the hat matrix \eqn{S_\theta}:
#'   \deqn{\hat y^{(-i)}_{\theta}(x_i) = \frac{\hat y_i - S_{ii} y_i}{1 - S_{ii}},}
#'   where \eqn{S_{ii}} is estimated by a finite-difference probe of the
#'   smoother (perturb each \eqn{y_i} by a small \eqn{\epsilon} and read
#'   off the resulting change in \eqn{\hat y_i}). The peak-preservation
#'   constraint vetoes any candidate whose smoothed global maximum is less
#'   than 95% of the raw global maximum.
#'
#' @section Limitations:
#'
#'   - Smoothing biases downstream inference. A smoothed peak height is a
#'     downward-biased estimate of the true peak height; the bias is
#'     largest for narrow peaks and aggressive smoothing. The same applies
#'     to peak widths (upward-biased) and peak areas (downward-biased for
#'     narrow peaks). Fit models on the raw data when the analyst can.
#'   - The adaptive PRESS criterion tends to pick small windows because it
#'     rewards predictive accuracy, which is dominated by the noise level.
#'     The peak-preservation constraint partially mitigates this, but the analyst
#'     should confirm the choice is reasonable for the application.
#'   - Savitzky-Golay with `poly_order` close to `window - 1` will
#'     overfit. The documentation recommends `poly_order <= 4` and `window >= 5`.
#'   - For non-uniformly sampled `x`, `sgolay` uses the index position
#'     (not the x value) as the abscissa, which introduces a small bias if
#'     `x` is very unevenly spaced. Use `method = "spline"` in that case.
#'
#' @references
#'   - Savitzky, A. and Golay, M. J. E. (1964). Smoothing and
#'     Differentiation of Data by Simplified Least Squares Procedures.
#'     Analytical Chemistry, 36, 1627-1639.
#'     \doi{10.1021/ac60214a047}.
#'
#'   - Steinier, J., Termonia, Y. and Deltour, J. (1972). Smoothing and
#'     differentiation of data by simplified least square procedure.
#'     Analytical Chemistry, 44, 1906-1909.
#'     \doi{10.1021/ac60319a045}.
#'
#'   - Cleveland, W. S. (1979). Robust Locally Weighted Regression and
#'     Smoothing Scatterplots. Journal of the American Statistical
#'     Association, 74, 829-836.
#'     \doi{10.1080/01621459.1979.10481038}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 200)
#' y <- 5 * exp(-(x - 5)^2 / 0.1) + rnorm(length(x), sd = 0.5)
#' fit <- adaptive_smoothing(x, y, method = "sgolay", poly_order = 3L)
#' print(fit)
#'
#' fit2 <- adaptive_smoothing(x, y, method = "spline")
#' print(fit2)
#'
#' fit3 <- adaptive_smoothing(x, y, method = "lowess", span = 0.2)
#' print(fit3)
#'
#' @export
adaptive_smoothing <- function(x, y, method = c("sgolay", "spline", "lowess"),
                              window = NULL, poly_order = 3L, spar = NULL,
                              span = 0.5, preserve_peaks = TRUE) {
  call <- match.call()
  method <- match.arg(method)
  .assert_numeric_vector(x, "x", min_length = 5L)
  .assert_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  .assert_positive_integer(poly_order, "poly_order")
  .assert_flag(preserve_peaks, "preserve_peaks")
  if (!is.null(spar)) {
    if (!is.numeric(spar) || length(spar) != 1L ||
        !is.finite(spar) || spar < 0 || spar > 1) {
      stop("'spar' must be a single numeric in [0, 1] or NULL.",
           call. = FALSE)
    }
  }
  .assert_probability(span, "span")

  if (method == "sgolay") {
    candidates_window <- c(5L, 7L, 9L, 11L, 15L, 21L)
    candidates_window <- candidates_window[candidates_window <= length(x)]
    if (length(candidates_window) == 0L) {
      stop("'x' too short for Savitzky-Golay (need at least 5 points).",
           call. = FALSE)
    }
    if (!is.null(window)) {
      .assert_positive_integer(window, "window")
      if (window %% 2L == 0L) {
        stop("'window' must be odd for Savitzky-Golay.", call. = FALSE)
      }
      if (poly_order >= window) {
        stop("'poly_order' must be less than 'window'.", call. = FALSE)
      }
      chosen_window <- window
      smoothed <- .sgolay_smooth(y, chosen_window, poly_order)
      press <- .press_sgolay(y, chosen_window, poly_order)
    } else {
      grid <- data.frame(
        window = candidates_window,
        press = NA_real_,
        peak_height_loss = NA_real_,
        stringsAsFactors = FALSE
      )
      for (i in seq_along(candidates_window)) {
        w <- candidates_window[i]
        if (poly_order >= w) next
        smoothed_i <- .sgolay_smooth(y, w, poly_order)
        grid$press[i] <- .press_sgolay(y, w, poly_order)
        grid$peak_height_loss[i] <- .peak_loss(y, smoothed_i)
      }
      chosen_window <- .pick_window(grid, preserve_peaks)
      smoothed <- .sgolay_smooth(y, chosen_window, poly_order)
      press <- grid$press[grid$window == chosen_window]
    }
    out_window <- chosen_window
    out_spar  <- NA_real_
    out_span  <- NA_real_
  } else if (method == "spline") {
    smoothed <- .spline_smooth(x, y, spar)
    press <- .press_spline(x, y, spar)
    out_window <- NA_integer_
    out_spar  <- if (is.null(spar)) NA_real_ else spar
    out_span  <- NA_real_
  } else {
    # lowess
    smoothed <- .lowess_smooth(x, y, span)
    press <- .press_lowess(x, y, span)
    out_window <- NA_integer_
    out_spar  <- NA_real_
    out_span  <- span
  }

  peak_loss <- .peak_loss(y, smoothed)

  .crossstatr_result(
    list(
      x = x,
      y = y,
      smoothed = smoothed,
      method = method,
      window = out_window,
      spar = out_spar,
      span = out_span,
      press = press,
      peak_height_loss = peak_loss,
      call = call
    ),
    "crossstatr_smoothing"
  )
}

# Direct Savitzky-Golay smoothing: for each point, fit a polynomial of
# order poly_order to a window of size `window` centred at that point,
# return the fitted value at the centre. Reflection padding at boundaries.
# @noRd
.sgolay_smooth <- function(y, window, poly_order) {
  if (window %% 2L == 0L) stop("window must be odd", call. = FALSE)
  if (poly_order >= window) stop("poly_order must be < window", call. = FALSE)
  half <- (window - 1L) %/% 2L
  k <- (-half:half)
  A <- outer(k, 0L:poly_order, function(a, b) a^b)
  # Coefficient that picks the constant term: first row of (A'A)^-1 A'.
  M <- tryCatch(
    solve(crossprod(A), t(A)),
    error = function(e) MASS::ginv(crossprod(A)) %*% t(A)
  )
  kernel <- M[1L, ]
  n <- length(y)
  out <- numeric(n)
  # Pad with reflections at the boundaries.
  if (half > 0L) {
    y_pad <- c(rev(y[1L:half]), y, rev(y[(n - half + 1L):n]))
  } else {
    y_pad <- y
  }
  for (i in seq_len(n)) {
    win <- y_pad[i + (0L:(window - 1L))]
    out[i] <- sum(win * kernel)
  }
  out
}

# Spline smoothing wrapper. If spar is NULL, stats::smooth.spline picks.
# @noRd
.spline_smooth <- function(x, y, spar) {
  args <- list(x = x, y = y, all.knots = FALSE)
  if (!is.null(spar)) args$spar <- spar
  fit <- do.call(stats::smooth.spline, args)
  as.numeric(stats::predict(fit, x)$y)
}

# Lowess smoothing wrapper. Aligns back to input order.
# @noRd
.lowess_smooth <- function(x, y, span) {
  ord <- order(x)
  lo <- stats::lowess(x[ord], y[ord], f = span)
  out <- numeric(length(x))
  out[ord] <- lo$y
  out
}

# PRESS for Savitzky-Golay: hat-matrix approach via finite-difference probe.
# @noRd
.press_sgolay <- function(y, window, poly_order) {
  n <- length(y)
  eps <- 1e-6
  base <- .sgolay_smooth(y, window, poly_order)
  s_diag <- numeric(n)
  for (i in seq_len(n)) {
    yp <- y; yp[i] <- yp[i] + eps
    sp <- .sgolay_smooth(yp, window, poly_order)
    s_diag[i] <- (sp[i] - base[i]) / eps
  }
  # Leave-one-out prediction: (y_i - S_ii y_i) / (1 - S_ii)
  s_diag <- pmin(pmax(s_diag, 1e-10), 1 - 1e-10)
  loo_pred <- (base - s_diag * y) / (1 - s_diag)
  sum((y - loo_pred)^2, na.rm = TRUE)
}

# PRESS for smoothing spline.
# @noRd
.press_spline <- function(x, y, spar) {
  n <- length(y)
  args <- list(x = x, y = y, all.knots = FALSE)
  if (!is.null(spar)) args$spar <- spar
  fit <- do.call(stats::smooth.spline, args)
  smoothed <- as.numeric(stats::predict(fit, x)$y)
  eps <- 1e-6
  s_diag <- numeric(n)
  for (i in seq_len(n)) {
    yp <- y; yp[i] <- yp[i] + eps
    args2 <- args; args2$y <- yp
    fp <- do.call(stats::smooth.spline, args2)
    sp <- as.numeric(stats::predict(fp, x)$y)
    s_diag[i] <- (sp[i] - smoothed[i]) / eps
  }
  s_diag <- pmin(pmax(s_diag, 1e-10), 1 - 1e-10)
  loo_pred <- (smoothed - s_diag * y) / (1 - s_diag)
  sum((y - loo_pred)^2, na.rm = TRUE)
}

# PRESS for lowess. Same finite-difference approach.
# @noRd
.press_lowess <- function(x, y, span) {
  n <- length(y)
  smoothed <- .lowess_smooth(x, y, span)
  eps <- 1e-6
  s_diag <- numeric(n)
  for (i in seq_len(n)) {
    yp <- y; yp[i] <- yp[i] + eps
    sp <- .lowess_smooth(x, yp, span)
    s_diag[i] <- (sp[i] - smoothed[i]) / eps
  }
  s_diag <- pmin(pmax(s_diag, 1e-10), 1 - 1e-10)
  loo_pred <- (smoothed - s_diag * y) / (1 - s_diag)
  sum((y - loo_pred)^2, na.rm = TRUE)
}

# Peak height loss: fraction of raw max lost in smoothing.
# @noRd
.peak_loss <- function(y, smoothed) {
  raw_max <- max(y, na.rm = TRUE)
  sm_max <- max(smoothed, na.rm = TRUE)
  if (!is.finite(raw_max) || raw_max <= 0) return(0)
  max(0, (raw_max - sm_max) / raw_max)
}

# Pick best window from grid by min PRESS subject to peak preservation.
# @noRd
.pick_window <- function(grid, preserve_peaks) {
  ok <- is.finite(grid$press)
  if (preserve_peaks) {
    ok <- ok & (grid$peak_height_loss <= 0.05)
  }
  if (!any(ok)) ok <- is.finite(grid$press)
  if (!any(ok)) ok <- rep(TRUE, nrow(grid))
  sub <- grid[ok, , drop = FALSE]
  best_idx <- which.min(sub$press)[1L]
  as.integer(sub$window[best_idx])
}

#' @export
print.crossstatr_smoothing <- function(x, ...) {
  cat(sprintf("<crossStatR::adaptive_smoothing> method = %s\n", x$method))
  if (x$method == "sgolay") {
    cat(sprintf("window = %d\n", x$window))
  } else if (x$method == "spline") {
    cat(sprintf("spar = %s\n", if (is.na(x$spar)) "(auto)" else sprintf("%.3f", x$spar)))
  } else {
    cat(sprintf("span = %.3f\n", x$span))
  }
  cat(sprintf("PRESS = %.4f\n", x$press))
  cat(sprintf("Peak height loss = %.2f%%\n", 100 * x$peak_height_loss))
  invisible(x)
}
