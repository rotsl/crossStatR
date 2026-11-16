#' Statistical evidence assessment for signal peaks
#'
#' Quantifies the evidence that a candidate feature in a measured signal
#' is a true peak rather than a noise fluctuation. The function reports
#' an exploratory evidence summary (likelihood ratio and
#' an exploratory resampling tail fraction); it does not report Bayesian
#' posterior probabilities because doing so would require specifying a
#' prior, which the package does not presume.
#'
#' @param x numeric vector of x-values (e.g. wavelength, m/z, voltage).
#' @param y numeric vector of signal intensities.
#' @param peak_index integer vector of candidate peak positions (indices
#'   into `x` and `y`). If NULL, peaks are auto-detected via local maxima
#'   with prominence above a noise threshold.
#' @param baseline numeric vector of baseline estimates (same length as
#'   `y`); if NULL, a simple running-minimum baseline is used.
#' @param noise_sd estimate of noise SD. If NULL, estimated from the MAD
#'   of the baseline-subtracted signal in flat regions.
#' @param conf_level confidence level for the peak amplitude CI.
#' @param method `"lr"` (default, likelihood ratio) or `"snr"` (signal-to-
#'   noise threshold with Bonferroni correction).
#' @param n_perm number of permutations for the empirical null
#'   distribution of the peak height (default 999). Set to 0 to use the
#'   Gaussian null approximation.
#' @param seed optional integer seed.
#'
#' @return A `crossstatr_peak_evidence` object with per-peak evidence
#'   table.
#'
#' @section Methodology:
#'
#'   Candidate amplitude is the baseline-subtracted value
#'   \deqn{A = y_{i_0} - b_{i_0}.}
#'   The reported Gaussian likelihood-style score is
#'   \deqn{\Lambda = \exp(A^2 / (2\sigma^2)).}
#'   This alternative-to-null ratio is not the null-to-alternative ratio
#'   used in the conventional Wilks statistic. No chi-square calibration
#'   or exact false-positive guarantee is established for this score.
#'
#'   For resampling, a smoothing spline is retained, its residuals are
#'   shuffled, and a maximum local contrast across the supplied candidate
#'   indices is compared with the observed amplitude. The reported tail
#'   fraction uses a plus-one correction. The candidate set is not
#'   re-detected within each resample, and the observed and resampled
#'   amplitude definitions differ. Local simulations found poor null
#'   calibration and low detection rates; the result is exploratory.
#'
#'   The Gaussian SNR threshold uses \eqn{z_{1-\alpha/(2m)}} for
#'   \eqn{m} candidates. With `n_perm = 0`, the implemented tail formula
#'   divides by \eqn{m}; it is not a valid Bonferroni-adjusted p-value.
#'   Peak-location uncertainty is a threshold-crossing half-width, not a
#'   calibrated confidence interval. Wilks (1938) and Currie (1968) provide
#'   background on likelihood inference and analytical detection; neither
#'   validates this implementation.
#'
#' @section Limitations:
#'
#'   - The likelihood-ratio form assumes a known noise SD; using an
#'     estimated SD inflates the apparent evidence in small samples.
#'     Neither the Gaussian nor permutation output has established
#'     calibration for small samples.
#'   - Overlapping peaks are not modelled jointly; if two candidate
#'     peaks are within one FWHM of each other, the LR for each is
#'     optimistic. Use [weighted_peak_fit()] for proper deconvolution.
#'   - Auto-detected peaks at the signal boundaries are unreliable
#'     because the baseline estimate is poorly constrained there.
#'   - This is a frequentist evidence measure, not a posterior
#'     probability. If the analyst require Bayesian posterior probabilities,
#'     specify the prior explicitly and use a different tool.
#'
#' @references
#'   - Wilks, S. S. (1938). The Large-Sample Distribution of the Likelihood
#'     Ratio for Testing Composite Hypotheses. The Annals of Mathematical
#'     Statistics, 9, 60-62.
#'     \doi{10.1214/aoms/1177732360}.
#'
#'   - Currie, L. A. (1968). Limits for qualitative detection and
#'     quantitative determination. Application to radiochemistry.
#'     Analytical Chemistry, 40, 586-593.
#'     \doi{10.1021/ac60259a007}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(400, 800, by = 0.5)
#' y <- dnorm(x, mean = 550, sd = 8) * 5 + dnorm(x, mean = 700, sd = 5) * 3 +
#'      0.1 + rnorm(length(x), sd = 0.05)
#' res <- peak_evidence(x, y)
#' print(res)
#'
#' @export
peak_evidence <- function(x, y, peak_index = NULL, baseline = NULL,
                          noise_sd = NULL, conf_level = 0.95,
                          method = c("lr", "snr"),
                          n_perm = 999L, seed = NULL) {
  call <- match.call()
  .assert_numeric_vector(x, "x", min_length = 5L)
  .assert_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have equal length.", call. = FALSE)
  }
  method <- match.arg(method)
  .assert_probability(conf_level, "conf_level")

  # Baseline: simple running-min if not supplied
  if (is.null(baseline)) {
    k <- 7L  # half-window
    baseline <- stats::approx(x, .running_min(y, k), xout = x,
                              rule = 2)$y
  }

  y_bc <- y - baseline

  # Noise SD from MAD of low-amplitude regions
  if (is.null(noise_sd)) {
    noise_sd <- stats::mad(y_bc, constant = 1.4826)
    if (!is.finite(noise_sd) || noise_sd <= 0) {
      noise_sd <- stats::sd(y_bc)
    }
  }
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("Could not estimate a positive noise SD.", call. = FALSE)
  }

  # Detect peaks if not given
  if (is.null(peak_index)) {
    peak_index <- .detect_local_maxima(y_bc, noise_sd)
  }
  if (length(peak_index) == 0L) {
    return(.crossstatr_result(
      list(peaks = data.frame(index = integer(0), x = numeric(0),
                              amplitude = numeric(0), snr = numeric(0),
                              lr = numeric(0), log_lr = numeric(0),
                              fpp_perm = numeric(0),
                              ci_low = numeric(0), ci_high = numeric(0),
                              location_uncertainty = numeric(0)),
           method = "Peak evidence (none detected)",
           noise_sd = noise_sd,
           n_perm = n_perm,
           conf_level = conf_level,
           call = call),
      "crossstatr_peak_evidence"))
  }

  m <- length(peak_index)
  alpha <- 1 - conf_level

  # Compute amplitude, SNR, LR, FPP per peak
  amps <- vapply(peak_index, function(idx) {
    # Local baseline = min of two side windows
    lo <- max(1L, idx - 5L); hi <- min(length(y), idx + 5L)
    side_b <- min(mean(baseline[max(1L, idx - 10L):(idx - 1L)]),
                  mean(baseline[(idx + 1L):min(length(y), idx + 10L)]))
    y_bc[idx]  # use baseline-subtracted value
  }, numeric(1L))

  # Bonferroni-corrected SNR threshold
  snr_thresh <- stats::qnorm(1 - alpha / (2 * m))

  # Permutation null: shuffle residuals, recompute max amplitude in
  # each candidate's neighbourhood.
  if (n_perm > 0L) {
    .set_seed(seed)
    # Build residual by subtracting a smooth signal estimate
    smooth <- stats::smooth.spline(x, y_bc, spar = 0.5)$y
    resid <- y_bc - smooth
    perm_max_amp <- replicate(n_perm, {
      sp <- sample(resid)
      yc <- smooth + sp
      max(yc[peak_index] - pmin(
        approx(x, yc, xout = x[peak_index] - 5, rule = 2)$y,
        approx(x, yc, xout = x[peak_index] + 5, rule = 2)$y),
        na.rm = TRUE)
    })
  } else {
    perm_max_amp <- NULL
  }

  # LR: lambda = exp(A^2 / (2 sigma^2)); log LR = A^2/(2 sigma^2)
  log_lr <- amps^2 / (2 * noise_sd^2)
  lr <- exp(log_lr)

  # FPP via permutation, or Gaussian approximation
  if (n_perm > 0L) {
    fpp <- vapply(amps, function(A) {
      (1 + sum(perm_max_amp >= A)) / (n_perm + 1)
    }, numeric(1L))
  } else {
    # Gaussian tail probability with Bonferroni correction
    fpp <- 2 * (1 - stats::pnorm(abs(amps) / noise_sd)) / m
    fpp <- pmin(fpp, 1)
  }

  # CI on amplitude: A ± z_{1-alpha/2} sigma
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci_low  <- amps - zcrit * noise_sd
  ci_high <- amps + zcrit * noise_sd

  # Location uncertainty: width at which SNR drops below snr_thresh
  loc_unc <- vapply(peak_index, function(idx) {
    # Find indices where y_bc drops below thresh*sigma around peak
    left  <- idx; right <- idx
    while (left > 1L && y_bc[left] > snr_thresh * noise_sd) left <- left - 1L
    while (right < length(y_bc) && y_bc[right] > snr_thresh * noise_sd) right <- right + 1L
    (x[right] - x[left]) / 2
  }, numeric(1L))

  peaks_df <- data.frame(
    index = peak_index,
    x = x[peak_index],
    amplitude = amps,
    snr = amps / noise_sd,
    lr = lr,
    log_lr = log_lr,
    fpp_perm = fpp,
    ci_low = ci_low,
    ci_high = ci_high,
    location_uncertainty = loc_unc,
    signif = fpp < alpha
  )

  .crossstatr_result(
    list(peaks = peaks_df,
         method = if (method == "lr")
           "Peak evidence via likelihood ratio + permutation FPP"
         else sprintf("Peak evidence via SNR > z_{1-a/(2m)} = %.3f (Bonferroni, m=%d)",
                       snr_thresh, m),
         noise_sd = noise_sd,
         n_perm = n_perm,
         conf_level = conf_level,
         call = call),
    "crossstatr_peak_evidence")
}

.running_min <- function(y, half_window) {
  n <- length(y)
  out <- y
  k <- half_window
  for (i in seq_len(n)) {
    lo <- max(1L, i - k); hi <- min(n, i + k)
    out[i] <- min(y[lo:hi])
  }
  out
}

.detect_local_maxima <- function(y, noise_sd) {
  n <- length(y)
  if (n < 3L) return(integer(0))
  is_max <- (y[2L:(n - 1L)] > y[1L:(n - 2L)]) &
            (y[2L:(n - 1L)] >= y[3L:n])
  idx <- which(is_max) + 1L
  # Filter by prominence
  idx <- idx[y[idx] > 3 * noise_sd]
  idx
}

#' @export
print.crossstatr_peak_evidence <- function(x, ...) {
  cat(sprintf("<crossStatR::peak_evidence> %s\n", x$method))
  cat(sprintf("Noise SD: %.4g, n_perm: %d, conf_level: %.0f%%\n",
              x$noise_sd, x$n_perm, 100 * x$conf_level))
  if (nrow(x$peaks) == 0L) {
    cat("No peaks detected.\n")
  } else {
    cat("\nPer-peak evidence:\n")
    print(x$peaks, row.names = FALSE)
  }
  invisible(x)
}
