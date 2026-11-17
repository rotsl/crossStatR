#' Reliability analysis for analytical detection
#'
#' Produces concentration-dependent detection reliability curves:
#' probability of detection (PD), false-positive probability (FPR), and
#' false-negative probability (FNR), with bootstrap confidence bands. This
#' is the proper frequentist companion to a single-point LOD: it shows how
#' detection performance varies across concentrations, rather than reducing
#' the whole curve to one number.
#'
#' @param calibration an `lm` object from `lm(y ~ x)`, OR a list with
#'   elements `x`, `y`, and `coef`.
#' @param blank_data optional numeric vector of blank measurements. If
#'   supplied, FPR is estimated empirically; otherwise, it is computed
#'   from the calibration residual SD assuming Gaussian noise.
#' @param conc_grid numeric vector of concentrations at which to
#'   evaluate the reliability curve. If NULL, defaults to 50 points
#'   spanning 0 to the upper calibration concentration.
#' @param decision_threshold the response threshold above which a sample
#'   is called "detected". If NULL, a one-sided prediction limit at
#'   `conf_level` is used for supplied blanks; otherwise the intercept plus
#'   `qnorm(conf_level)` times the effective noise SD is used.
#' @param conf_level confidence level for the reliability curves.
#' @param n_boot number of bootstrap replicates for uncertainty bands
#'   (default 300).
#' @param seed optional integer seed.
#' @param extra_sigma2 numeric scalar; additional variance to add to the
#'   calibration noise variance, e.g. from matrix effects, batch effects,
#'   or drift (see [batch_effect_quantify()], [drift_uncertainty()]).
#'   Defaults to 0.
#'
#' @return A `crossstatr_detection_reliability` object.
#'
#' @section Methodology:
#'
#'   Under the calibration model \eqn{y = a + b x + \epsilon} with
#'   \eqn{\epsilon \sim N(0, \sigma^2)}, the distribution of a measured
#'   response at concentration \eqn{x_0} is
#'   \deqn{y(x_0) \sim N(a + b x_0, \sigma^2).}
#'
#'   Given a decision threshold \eqn{L_C} (the critical level), the
#'   false-positive probability is
#'   \deqn{\mathrm{FPR} = P(y > L_C \mid x = 0)
#'     = 1 - \Phi\!\left(\frac{L_C - a}{\sigma}\right),}
#'   and the false-negative probability at concentration \eqn{x_0} is
#'   \deqn{\mathrm{FNR}(x_0) = P(y \leq L_C \mid x = x_0)
#'     = \Phi\!\left(\frac{L_C - a - b x_0}{\sigma}\right).}
#'
#'   The detection probability is
#'   \deqn{\mathrm{PD}(x_0) = 1 - \mathrm{FNR}(x_0).}
#'
#'   The classical LOD is the concentration where PD = 0.95
#'   (equivalently, FNR = 0.05), which under Gaussian assumptions is
#'   \deqn{x_{\mathrm{LOD}} = (L_C + z_{0.95}\,\sigma - a) / b.}
#'
#'   Matrix effects, batch effects, and drift can be incorporated by
#'   inflating \eqn{\sigma} additively:
#'   \deqn{\sigma^2_{\mathrm{eff}} = \sigma^2_{\mathrm{noise}} +
#'     \sigma^2_{\mathrm{matrix}} + \sigma^2_{\mathrm{batch}} +
#'     \sigma^2_{\mathrm{drift}},}
#'   provided those variance components have been estimated separately
#'   (e.g. via [batch_effect_quantify()], [drift_uncertainty()]). Pass
#'   these as additional variance contributions via `extra_sigma2`.
#'
#'   Bootstrap CIs on PD/FNR are obtained by sampling calibration
#'   coefficients from their MVN distribution, varying residual noise, and
#'   resampling supplied blank measurements when the threshold is estimated.
#'   Percentile bands and `lod_ci` are returned. Additional variance supplied
#'   through `extra_sigma2` is held fixed.
#'
#' @section Interpretation:
#'
#'   A single LOD summarises one decision threshold. The probability
#'   curve describes detection across concentrations under the fitted
#'   noise model; Currie (1999) supplies the analytical detection framework.
#'
#' @section Limitations:
#'
#'   - Only estimable positive-slope linear calibrations are supported.
#'   - Bands are percentile approximations; no BCa correction is applied.
#'     Small bootstrap counts give unstable tail quantiles. A fixed supplied
#'     threshold is held constant across bootstrap draws.
#'
#'   - The Gaussian noise assumption must hold. Heavy-tailed noise (e.g.
#'     from sporadic contamination) will inflate FPR beyond the nominal
#'     rate.
#'   - The default decision threshold is computed from the calibration
#'     data; if the user supplies a different one, it must be reported
#'     alongside the curves.
#'   - Reliability curves assume the calibration is valid across the
#'     entire `conc_grid` range. Extrapolating below the lowest
#'     calibration concentration is risky.
#'   - The bootstrap CI assumes the calibration model is correctly
#'     specified.
#'
#' @references
#'   - Currie, L. A. (1999). Detection and quantification limits: origins
#'     and historical overview. Analytica Chimica Acta, 391, 127-134.
#'     \doi{10.1016/S0003-2670(99)00105-1}.
#'
#'   - Currie, L. A. (1995). Nomenclature in evaluation of analytical
#'     methods including detection and quantification capabilities (IUPAC
#'     Recommendations 1995). Pure and Applied Chemistry, 67, 1699-1723.
#'     \doi{10.1351/pac199567101699}.
#'
#'   - Hubaux, A. and Vos, G. (1970). Decision and detection limits for
#'     calibration curves. Analytical Chemistry, 42, 849-855.
#'     \doi{10.1021/ac60290a013}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(0.1, 10, length.out = 25)
#' y <- 0.5 + 0.9 * x + rnorm(length(x), sd = 0.3)
#' cal <- lm(y ~ x)
#' res <- detection_reliability(cal, blank_data = rnorm(20, 0.5, 0.3))
#' print(res)
#'
#' @export
detection_reliability <- function(calibration, blank_data = NULL,
                                   conc_grid = NULL, decision_threshold = NULL,
                                   conf_level = 0.95, n_boot = 300L,
                                   seed = NULL, extra_sigma2 = 0) {
  call <- match.call()
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(n_boot, "n_boot")

  # Extract calibration parameters
  if (inherits(calibration, "lm")) {
    beta <- stats::coef(calibration)
    Vbeta <- stats::vcov(calibration)
    sigma <- summary(calibration)$sigma
    df <- calibration$df.residual
    x_cal <- calibration$model[, 2L]
  } else if (is.list(calibration) && all(c("x","y","coef") %in% names(calibration))) {
    fit <- stats::lm(calibration$y ~ calibration$x)
    beta <- stats::coef(fit)
    Vbeta <- stats::vcov(fit)
    sigma <- summary(fit)$sigma
    df <- fit$df.residual
    x_cal <- calibration$x
  } else {
    stop("'calibration' must be an lm object or a list with x, y, coef.",
         call. = FALSE)
  }

  if (length(beta) != 2L || any(!is.finite(beta)) || beta[2L] <= 0 || df <= 0) {
    stop("Detection reliability requires an estimable positive-slope linear calibration.", call. = FALSE)
  }
  if (!is.numeric(extra_sigma2) || length(extra_sigma2) != 1L ||
      !is.finite(extra_sigma2) || extra_sigma2 < 0) {
    stop("'extra_sigma2' must be a finite non-negative variance.", call. = FALSE)
  }
  if (!is.null(blank_data)) .assert_numeric_vector(blank_data, "blank_data", min_length = 2L)
  sigma_eff <- sqrt(sigma^2 + extra_sigma2)

  # Decision threshold: upper 95% prediction bound of the blank
  if (is.null(decision_threshold)) {
    if (!is.null(blank_data)) {
      m_blank <- mean(blank_data)
      s_blank <- stats::sd(blank_data)
      n_b <- length(blank_data)
      tcrit <- stats::qt(1 - (1 - conf_level), df = n_b - 1)
      L_C <- m_blank + tcrit * s_blank * sqrt(1 + 1 / n_b)
    } else {
      # Calibration-based: y at x=0 is just the intercept
      zc <- stats::qnorm(1 - (1 - conf_level))
      L_C <- beta[1L] + zc * sigma_eff
    }
  } else {
    L_C <- decision_threshold
  }

  # Grid
  if (is.null(conc_grid)) {
    conc_grid <- seq(0, max(x_cal), length.out = 50L)
  }

  # Reliability curves
  mu_at_x <- beta[1L] + beta[2L] * conc_grid
  pd <- 1 - stats::pnorm((L_C - mu_at_x) / sigma_eff)
  fnr <- 1 - pd
  fpr <- 1 - stats::pnorm((L_C - beta[1L]) / sigma_eff)

  # LOD: smallest x where PD >= 0.95
  lod <- (L_C + stats::qnorm(0.95) * sigma_eff - beta[1L]) / beta[2L]

  # Bootstrap CIs
  .set_seed(seed)
  boot_pd <- matrix(NA_real_, nrow = n_boot, ncol = length(conc_grid))
  boot_fpr <- numeric(n_boot)
  boot_lod <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    b_beta <- MASS::mvrnorm(1L, mu = beta, Sigma = Vbeta)
    # Optionally vary sigma_eff via chi-square draw
    b_sigma <- sqrt(sigma^2 * stats::rchisq(1L, df = df) / df + extra_sigma2)
    if (!is.null(blank_data) && is.null(decision_threshold)) {
      b_blank <- sample(blank_data, replace = TRUE)
      b_m_b <- mean(b_blank)
      b_s_b <- stats::sd(b_blank)
      b_tcrit <- stats::qt(conf_level, df = n_b - 1L)
      b_LC <- b_m_b + b_tcrit * b_s_b * sqrt(1 + 1 / n_b)
    } else {
      b_LC <- L_C  # use supplied/decision threshold; bootstrap calibration only
    }
    b_mu <- b_beta[1L] + b_beta[2L] * conc_grid
    boot_pd[b, ] <- 1 - stats::pnorm((b_LC - b_mu) / b_sigma)
    boot_fpr[b]  <- 1 - stats::pnorm((b_LC - b_beta[1L]) / b_sigma)
    boot_lod[b]  <- (b_LC + stats::qnorm(0.95) * b_sigma - b_beta[1L]) / b_beta[2L]
  }
  alpha <- 1 - conf_level
  pd_ci_low  <- apply(boot_pd, 2L, stats::quantile, probs = alpha / 2,
                     na.rm = TRUE)
  pd_ci_high <- apply(boot_pd, 2L, stats::quantile, probs = 1 - alpha / 2,
                     na.rm = TRUE)
  fpr_ci <- stats::quantile(boot_fpr, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  lod_ci <- stats::quantile(boot_lod, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  # Empirical FPR if blank data supplied
  fpr_emp <- if (!is.null(blank_data)) mean(blank_data > L_C) else NA_real_

  curve_df <- data.frame(
    concentration = conc_grid,
    mean_response = mu_at_x,
    pd = pd,
    pd_ci_low = pd_ci_low,
    pd_ci_high = pd_ci_high,
    fnr = fnr
  )

  .crossstatr_result(
    list(curve = curve_df,
         decision_threshold = L_C,
         sigma_eff = sigma_eff,
         fpr = fpr,
         fpr_ci = fpr_ci,
         fpr_empirical = fpr_emp,
         lod = lod,
         lod_ci = lod_ci,
         conf_level = conf_level,
         n_boot = n_boot,
         method = "Detection reliability curve",
         call = call),
    "crossstatr_detection_reliability")
}

#' @export
print.crossstatr_detection_reliability <- function(x, ...) {
  cat(sprintf("<crossStatR::detection_reliability> %s\n", x$method))
  cat(sprintf("Decision threshold L_C = %.4g, sigma_eff = %.4g\n",
              x$decision_threshold, x$sigma_eff))
  cat(sprintf("LOD (PD=0.95): %.4g\n", x$lod))
  cat(sprintf("Theoretical FPR: %.4f (95%% CI [%.4f, %.4f])\n",
              x$fpr, x$fpr_ci[1L], x$fpr_ci[2L]))
  if (!is.na(x$fpr_empirical)) {
    cat(sprintf("Empirical FPR (from blank data): %.4f\n", x$fpr_empirical))
  }
  cat("\nReliability curve (first 6 + last 6 rows):\n")
  curve <- x$curve
  cat("Low concentration:\n"); print(utils::head(curve, 6L), row.names = FALSE)
  cat("High concentration:\n"); print(utils::tail(curve, 6L), row.names = FALSE)
  invisible(x)
}
