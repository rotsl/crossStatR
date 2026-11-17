#' Characterize measurement noise
#'
#' Decomposes a sampled signal into a smooth trend and a residual noise
#' component, then summarises the noise: overall standard deviation, lag-1
#' through lag-10 autocorrelation, heteroscedasticity (correlation of
#' absolute residuals with the fitted trend), a power spectral density
#' estimate via FFT with a dominant-frequency pick, and signal-dependent
#' noise (residual SD within bins of the fitted signal level).
#'
#' Two trend-removal strategies are available. `method = "residual"` fits a
#' polynomial of `degree` and takes the residuals. `method = "difference"`
#' uses successive differences `y[i+1] - y[i]` as a noise proxy, which is
#' appropriate when the signal is smooth relative to the sampling rate so
#' that differences are dominated by noise rather than signal curvature.
#'
#' @param x numeric vector of abscissae (time, wavelength, etc.).
#' @param y numeric vector of measurements (same length as `x`).
#' @param method `"residual"` (default) or `"difference"`.
#' @param degree integer; polynomial degree for `method = "residual"`.
#' @param conf_level confidence level used for the heteroscedasticity test
#'   and the dominant-frequency CI label.
#'
#' @return A `crossstatr_noise` object with elements:
#'   \item{noise_sd}{overall standard deviation of the residual series.}
#'   \item{acf}{named numeric vector of autocorrelations at lags 1 to 10.}
#'   \item{heteroscedasticity_r}{Pearson correlation of `|residual|` with
#'     the fitted signal.}
#'   \item{heteroscedasticity_p}{two-sided p-value of the correlation test.}
#'   \item{dominant_freq}{frequency (in cycles per unit of `x`) with the
#'     largest PSD component (excluding DC).}
#'   \item{psd}{data frame of `freq` and `power` (positive frequencies only).}
#'   \item{signal_dependent}{data frame of `bin`, `mean_signal`,
#'     `sd_noise` summarising residual SD within bins of the fitted trend.}
#'   \item{method}{method label.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   For `method = "residual"` The package fits \eqn{y_i = \beta_0 + \sum_{k=1}^{d}
#'   \beta_k x_i^k + \epsilon_i} by ordinary least squares. The residual
#'   series \eqn{e_i = y_i - \hat{y}_i} is taken as the noise realisation.
#'
#'   For `method = "difference"` The package uses \eqn{e_i = y_{i+1} - y_i}. If the
#'   underlying signal is locally linear and the noise has variance
#'   \eqn{\sigma_\epsilon^2}, then \eqn{\mathrm{Var}(e_i) \approx
#'   2\sigma_\epsilon^2}, so the original noise SD is approximately
#'   \eqn{\mathrm{sd}(e)/\sqrt{2}}.
#'
#'   Autocorrelation at lag \eqn{k}:
#'   \deqn{\rho_k = \frac{\mathrm{Cov}(e_i, e_{i+k})}{\mathrm{Var}(e_i)}.}
#'
#'   Heteroscedasticity is tested by the Pearson correlation between
#'   \eqn{|e_i|} and \eqn{\hat{y}_i}. This is an exploratory correlation diagnostic,
#'   not the White or studentized Breusch-Pagan heteroscedasticity test.
#'
#'   The power spectral density is computed via the periodogram
#'   \deqn{P(\omega) = \frac{1}{n}\left|\sum_{i=1}^{n}e_i e^{-i 2\pi
#'   \omega x_i / T}\right|^2,}
#'   with `T` the span of `x`. The dominant non-DC frequency is reported.
#'
#'   For signal-dependent noise the fitted values are binned into
#'   approximately 10 equal-probability bins, and the residual SD is
#'   computed within each bin.
#'
#' @section Limitations:
#'
#'   - The polynomial trend is global; for sharply varying signals it can
#'     over-fit in some regions and under-fit in others. A smoothing
#'     spline or local polynomial regression is preferable for non-smooth
#'     trends; this function does not yet offer that option.
#'   - The difference method assumes the signal is approximately linear
#'     between consecutive samples. If the sampling is sparse relative to
#'     signal curvature, the differences will be dominated by signal, not
#'     noise, and `noise_sd` will be inflated.
#'   - The periodogram has high variance as a spectral estimator; the
#'     dominant-frequency pick is noisy for short series. Use Welch's
#'     method or a multitaper estimator for serious spectral analysis.
#'   - Heteroscedasticity as tested here only flags a monotone relation
#'     between residual magnitude and signal level; it does not detect
#'   other patterns (e.g. variance depending on time).
#'
#' @references
#'   - Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (2008). Time Series
#'     Analysis: Forecasting and Control. 4th ed. Wiley.
#'     \doi{10.1002/9781118619193}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 200)
#' y <- sin(x) + rnorm(x, sd = 0.1 + 0.05 * abs(sin(x)))
#' result <- noise_characterize(x, y, method = "residual", degree = 5L)
#' print(result)
#' head(result$psd)
#' head(result$signal_dependent)
#'
#' # Difference method
#' noise_characterize(x, y, method = "difference")
#'
#' @export
noise_characterize <- function(x, y,
                                method = c("residual", "difference"),
                                degree = 3L, conf_level = 0.95) {
  call <- match.call()
  method <- match.arg(method)
  .assert_numeric_vector(x, "x", min_length = 5L)
  .assert_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have equal length.", call. = FALSE)
  }
  .assert_positive_integer(degree, "degree")
  .assert_probability(conf_level, "conf_level")

  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L) {
    stop("Need at least 5 finite (x, y) pairs.", call. = FALSE)
  }
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n <= degree + 1L) {
    stop(sprintf("n (%d) must exceed degree + 1 (%d) for residual fit.",
                 n, degree + 1L), call. = FALSE)
  }

  dt <- mean(diff(x), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) dt <- 1

  if (method == "residual") {
    fit <- stats::lm(y ~ poly(x, degree, raw = TRUE))
    fitted_vals <- stats::fitted(fit)
    resids <- stats::resid(fit)
    method_label <- sprintf("Polynomial regression (degree %d) residuals",
                            degree)
    resids_x <- x
  } else {
    fitted_vals <- (y[-n] + y[-1L]) / 2
    resids <- diff(y)
    resids_x <- (x[-n] + x[-1L]) / 2
    method_label <- "Successive differences (y[i+1] - y[i])"
    n <- n - 1L
  }

  noise_sd <- .safe_sd(resids)

  # Autocorrelation at lags 1..10 (or fewer for short series)
  max_lag <- min(10L, n - 1L)
  acf_vals <- tryCatch(
    stats::acf(resids, plot = FALSE, lag.max = max_lag,
               na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (!is.null(acf_vals) && length(acf_vals$acf) > 1L) {
    v <- as.vector(acf_vals$acf)[-1L]
    acf_named <- setNames(v, paste0("lag", seq_along(v)))
  } else {
    acf_named <- setNames(rep(NA_real_, max_lag),
                          paste0("lag", seq_len(max_lag)))
  }

  # Heteroscedasticity: correlation of |resid| with fitted
  if (length(fitted_vals) == length(resids) &&
      stats::sd(abs(resids)) > 0 && stats::sd(fitted_vals) > 0) {
    ct <- stats::cor.test(abs(resids), fitted_vals, method = "pearson")
    het_r <- unname(ct$estimate)
    het_p <- ct$p.value
  } else {
    het_r <- NA_real_
    het_p <- NA_real_
  }

  # PSD via FFT (periodogram)
  n_fft <- stats::nextn(n, 2L)
  e_centered <- resids - mean(resids, na.rm = TRUE)
  e_padded <- c(e_centered, rep(0, n_fft - n))
  fft_res <- stats::fft(e_padded)
  psd_all <- (Mod(fft_res)^2) / n  # normalise by n, not n_fft
  n_pos <- n_fft %/% 2L
  freqs <- seq(0, 1, length.out = n_fft + 1L)[1:n_pos] / dt
  psd_pos <- psd_all[1:n_pos]
  if (length(psd_pos) > 2L) {
    dom_idx <- which.max(psd_pos[-1L]) + 1L
    dominant_freq <- freqs[dom_idx]
  } else {
    dominant_freq <- NA_real_
  }
  psd_df <- data.frame(freq = freqs, power = psd_pos,
                       stringsAsFactors = FALSE)

  # Signal-dependent noise: bin fitted values, compute SD of residuals
  n_bins <- max(3L, min(10L, floor(n / 5)))
  qs <- stats::quantile(fitted_vals, probs = seq(0, 1, length.out = n_bins + 1L),
                        na.rm = TRUE)
  qs <- unique(qs)
  if (length(qs) < 3L) {
    sig_dep <- data.frame(bin = 1L, mean_signal = mean(fitted_vals, na.rm = TRUE),
                         sd_noise = noise_sd, n = as.integer(n),
                         stringsAsFactors = FALSE)
  } else {
    bin_idx <- cut(fitted_vals, breaks = qs, include.lowest = TRUE)
    sig_dep <- do.call(rbind, lapply(seq_along(levels(bin_idx)), function(b) {
      sel <- bin_idx == levels(bin_idx)[b] & !is.na(bin_idx)
      if (sum(sel) < 2L) return(NULL)
      data.frame(bin = b,
                 mean_signal = mean(fitted_vals[sel], na.rm = TRUE),
                 sd_noise = .safe_sd(resids[sel]),
                 n = as.integer(sum(sel)),
                 stringsAsFactors = FALSE)
    }))
  }

  tidy_df <- data.frame(
    noise_sd = noise_sd,
    heteroscedasticity_r = het_r,
    heteroscedasticity_p = het_p,
    dominant_freq = dominant_freq,
    method = method_label,
    stringsAsFactors = FALSE
  )

  .crossstatr_result(
    list(noise_sd = noise_sd,
         acf = acf_named,
         heteroscedasticity_r = het_r,
         heteroscedasticity_p = het_p,
         dominant_freq = dominant_freq,
         psd = psd_df,
         signal_dependent = sig_dep,
         method = method_label,
         fitted = fitted_vals,
         residuals = resids,
         conf_level = conf_level,
         tidy = tidy_df,
         call = call),
    "crossstatr_noise"
  )
}

#' @export
print.crossstatr_noise <- function(x, ...) {
  cat(sprintf("<crossStatR::noise_characterize> %s\n", x$method))
  cat(sprintf("noise_sd = %.4g\n", x$noise_sd))
  cat("Autocorrelation (lags 1-10):")
  ac <- x$acf
  cat(sprintf(" %s\n",
              paste(sprintf("%.3f", ac), collapse = " ")))
  cat(sprintf("heteroscedasticity r = %.3f (p = %g)\n",
              x$heteroscedasticity_r, x$heteroscedasticity_p))
  cat(sprintf("dominant frequency = %.4g\n", x$dominant_freq))
  cat("\nSignal-dependent noise:\n")
  print(x$signal_dependent, row.names = FALSE)
  invisible(x)
}
