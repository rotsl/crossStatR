#' Limits of detection and quantification
#'
#' Computes the limit of detection (LOD) and limit of quantification (LOQ)
#' for a calibration curve, using one of three classical methods: the
#' Currie (1968) critical-value approach, the simple blank-standard-deviation
#' approach, or the Hubaux-Vos (1970) prediction-interval approach. For
#' heteroscedastic calibration data the residual SD can be modelled as a
#' function of concentration via lowess, in which case the LOD is computed
#' at the local residual SD evaluated at the LOD itself.
#'
#' The calibration can be supplied either as an `lm` object fit by
#' `stats::lm(y ~ x)`, or as a list with components `x`, `y`, and
#' (optionally) `coef`. If `response` and `conc` are supplied, they
#' override the calibration curve as the source of the blank SD.
#'
#' @param calibration an `lm` object (from `stats::lm(y ~ x)`) or a list
#'   with components `x` (concentration), `y` (response), and optionally
#'   `coef` (length-2: intercept, slope). If `coef` is missing it is fit
#'   by OLS.
#' @param response optional; a numeric vector of blank-sample responses.
#'   Used to compute the blank SD for the Currie and blank-SD methods.
#'   If NULL, the calibration residual SD is used instead.
#' @param conc optional; a numeric vector of concentrations corresponding
#'   to `response`. If supplied with `response`, used for heteroscedastic
#'   noise modelling.
#' @param conf_level confidence level for the LOD CI.
#' @param false_positive_rate target false-positive rate (alpha). Default 0.05.
#' @param false_negative_rate target false-negative rate (beta). Default 0.05.
#'   The classical Currie 3*sigma LOD corresponds to alpha = beta ~ 0.00135
#'   (one-sided z = 3). The default 0.05 gives a more conservative LOD.
#' @param method character; `"currie"` (default), `"blank_sd"`, or
#'   `"hubaux_vos"`.
#'
#' @return A `crossstatr_lod` object with elements:
#'   \item{lod}{limit of detection.}
#'   \item{loq}{limit of quantification.}
#'   \item{lod_ci}{Monte Carlo CI on the LOD (length-2).}
#'   \item{method}{method used.}
#'   \item{blank_sd}{blank SD used in the calculation.}
#'   \item{slope, intercept}{calibration slope and intercept.}
#'   \item{residual_se}{residual standard error of the calibration fit.}
#'   \item{df_resid}{residual degrees of freedom.}
#'   \item{false_positive_rate, false_negative_rate}{as input.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   The calibration is a linear model \eqn{y = \alpha + \beta x + \epsilon}
#'   with \eqn{\epsilon \sim N(0, \sigma^2)}. The LOD is the lowest
#'   concentration that can be reliably distinguished from blank.
#'
#'   Currie (1968) critical value: the LOD is the concentration at which
#'   the probability of a false positive (signal exceeds blank threshold
#'   when true concentration is zero) and the probability of a false
#'   negative (signal below threshold when true concentration is LOD) are
#'   both controlled. The classical 3.3*sigma rule corresponds to
#'   \eqn{z_\alpha + z_\beta \approx 3.3} (with \eqn{\alpha = \beta \approx
#'   0.05}). The general form is
#'   \deqn{\mathrm{LOD} = \frac{(z_{1-\alpha} + z_{1-\beta})\, \sigma_b}{\beta},}
#'   where \eqn{\sigma_b} is the blank SD (or the calibration residual SD
#'   when no blanks are available) and \eqn{\beta} is the slope. The LOQ
#'   is conventionally
#'   \deqn{\mathrm{LOQ} = \frac{10\, \sigma_b}{\beta},}
#'   corresponding to a relative SD of 10% at the LOQ (Currie, 1968;
#'   IUPAC, 1995). The IUPAC simple form replaces \eqn{z_{1-\alpha} +
#'   z_{1-\beta}} with the constant 3.
#'
#'   Blank-SD method: same form as Currie but using the empirical SD of
#'   the blank sample supplied via `response`, with no normal-theory
#'   calibration-based SE. The classical 3*sigma rule is used.
#'
#'   Hubaux-Vos (1970): the LOD is the concentration at which the lower
#'   one-sided prediction bound on the response equals the upper one-sided
#'   prediction bound at zero concentration. The closed-form fixed-point
#'   equation is
#'   \deqn{\mathrm{LOD} = \frac{(t_{1-\alpha,\,\nu} + t_{1-\beta,\,\nu})\,
#'         \sigma}{\beta}\,
#'         \sqrt{1 + \frac{1}{n} + \frac{\mathrm{LOD}^2}{S_{xx}}},}
#'   where \eqn{\sigma} is the calibration residual SE, \eqn{\nu} is the
#'   residual df, \eqn{n} is the number of calibration points, and
#'   \eqn{S_{xx}} is the corrected sum of squares of the calibration
#'   concentrations. The equation is solved iteratively to convergence.
#'
#'   Heteroscedastic noise: if the residual SD varies appreciably with
#'   concentration (the ratio of high- to low-conc residual SD exceeds 2,
#'   assessed via lowess smoothing of \eqn{|r_i|} on \eqn{x_i}), the LOD
#'   is recomputed using the local \eqn{\sigma(x)} evaluated at the LOD
#'   estimate. This is iterated to convergence.
#'
#'   Monte Carlo CI on the LOD: simulate from the asymptotic distribution
#'   of \eqn{(\hat\alpha, \hat\beta, \hat\sigma)}, recompute the LOD for
#'   each simulation, and take the percentile interval.
#'
#' @section Limitations:
#'
#'   - All three methods assume the calibration is linear in the relevant
#'     concentration range. If the calibration is nonlinear near zero, the
#'     LOD is unreliable.
#'   - The Hubaux-Vos method assumes homoscedastic noise unless the
#'     heteroscedastic correction is invoked; even then the lowess fit is
#'     a smoother, not a model, and the resulting LOD is approximate.
#'   - The Monte Carlo CI assumes the asymptotic normal distribution of
#'     the calibration parameters. For small calibration samples, this
#'     underestimates the true uncertainty.
#'   - The 3*sigma convention is a rule of thumb, not a guarantee of any
#'     specific false-positive rate. Setting `false_positive_rate` and
#'     `false_negative_rate` explicitly is more honest.
#'
#' @references
#'   - Currie, L. A. (1968). Limits for qualitative detection and
#'     quantitative determination. Application to radiochemistry.
#'     Analytical Chemistry, 40, 586-593.
#'     \doi{10.1021/ac60259a007}.
#'
#'   - Hubaux, A. and Vos, G. (1970). Decision and detection limits for
#'     calibration curves. Analytical Chemistry, 42, 849-855.
#'     \doi{10.1021/ac60290a013}.
#'
#'   - Currie, L. A. (1995). Nomenclature in evaluation of analytical
#'     methods including detection and quantification capabilities (IUPAC
#'     Recommendations 1995). Pure and Applied Chemistry, 67, 1699-1723.
#'     \doi{10.1351/pac199567101699}.
#'
#' @examples
#' set.seed(1)
#' conc <- seq(0, 10, length.out = 20)
#' resp <- 0.5 + 2 * conc + rnorm(length(conc), sd = 0.3)
#' cal <- lm(resp ~ conc)
#' lod <- robust_lod(cal, method = "currie")
#' print(lod)
#'
#' # With explicit blank measurements
#' blank_responses <- rnorm(10, mean = 0.5, sd = 0.3)
#' lod2 <- robust_lod(cal, response = blank_responses, method = "currie")
#' print(lod2)
#'
#' # Hubaux-Vos prediction-interval approach
#' lod3 <- robust_lod(cal, method = "hubaux_vos")
#' print(lod3)
#'
#' @export
robust_lod <- function(calibration, response = NULL, conc = NULL,
                      conf_level = 0.95,
                      false_positive_rate = 0.05,
                      false_negative_rate = 0.05,
                      method = c("currie", "blank_sd", "hubaux_vos")) {
  call <- match.call()
  method <- match.arg(method)
  .assert_probability(conf_level, "conf_level")
  .assert_probability(false_positive_rate, "false_positive_rate")
  .assert_probability(false_negative_rate, "false_negative_rate")

  cal <- .coerce_calibration(calibration)
  slope <- cal$coef[2L]
  intercept <- cal$coef[1L]
  residual_se <- cal$residual_se
  df_resid <- cal$df_resid
  x_cal <- cal$x
  y_cal <- cal$y
  n <- length(x_cal)

  if (!is.null(response)) {
    .assert_numeric_vector(response, "response", min_length = 2L)
  }
  if (!is.null(conc)) {
    .assert_numeric_vector(conc, "conc", min_length = 2L)
    if (!is.null(response) && length(conc) != length(response)) {
      stop("'conc' and 'response' must have the same length.", call. = FALSE)
    }
  }

  # Blank SD: from explicit blanks, else from calibration residual SE.
  if (!is.null(response)) {
    blank_sd <- .safe_sd(response)
  } else {
    blank_sd <- residual_se
  }
  if (!is.finite(blank_sd) || blank_sd <= 0) {
    stop("Blank SD is zero or not finite; cannot compute LOD.",
         call. = FALSE)
  }
  if (!is.finite(slope) || abs(slope) < .Machine$double.eps) {
    stop("Calibration slope is zero or not finite; cannot compute LOD.",
         call. = FALSE)
  }

  # Detect heteroscedasticity from calibration residuals vs concentration.
  # The method smooths |r| on x to estimate E[|r| | x], then convert to a local
  # sigma estimate via sigma(x) = E[|r| | x] / sqrt(2/pi) (Gaussian case).
  # The ratio of max to min sigma(x) flags heteroscedasticity.
  cal_resid <- y_cal - (intercept + slope * x_cal)
  sigma_fun <- NULL
  hetero <- FALSE
  if (length(x_cal) >= 6L && .safe_sd(cal_resid) > 0) {
    abs_r <- abs(cal_resid)
    lo <- stats::lowess(x_cal, abs_r)
    sigma_lo <- lo$y * sqrt(pi / 2)
    # Floor the local sigma at the smallest positive observed |r| to avoid
    # driving the LOD to numerical zero during iteration.
    floor_val <- min(abs_r[abs_r > 0], na.rm = TRUE) * sqrt(pi / 2)
    if (!is.finite(floor_val) || floor_val <= 0) {
      floor_val <- .Machine$double.eps
    }
    sigma_lo[sigma_lo < floor_val] <- floor_val
    ratio <- max(sigma_lo, na.rm = TRUE) /
      max(min(sigma_lo, na.rm = TRUE), .Machine$double.eps)
    if (is.finite(ratio) && ratio > 2) {
      hetero <- TRUE
      sigma_fun <- stats::approxfun(lo$x, sigma_lo, rule = 2)
    }
  }

  z_alpha <- stats::qnorm(1 - false_positive_rate)
  z_beta  <- stats::qnorm(1 - false_negative_rate)

  # Initial LOD/LOQ with the global blank SD.
  lod_out <- .compute_lod(method, slope, intercept, residual_se,
                          blank_sd, x_cal, n, df_resid,
                          false_positive_rate, false_negative_rate,
                          z_alpha, z_beta)
  lod <- lod_out$lod
  loq <- 10 * blank_sd / abs(slope)

  # Heteroscedastic refinement: substitute local sigma at the LOD estimate
  # and iterate. The method only adjusts when the LOD falls within the calibration
  # range so the lowess extrapolation is not relied on.
  if (hetero && is.finite(lod)) {
    x_range <- range(x_cal, na.rm = TRUE)
    for (iter in seq_len(50L)) {
      if (lod < x_range[1L] || lod > x_range[2L]) break
      sigma_local <- sigma_fun(lod)
      if (!is.finite(sigma_local) || sigma_local <= 0) break
      lod_new <- .compute_lod(method, slope, intercept, residual_se,
                              sigma_local, x_cal, n, df_resid,
                              false_positive_rate, false_negative_rate,
                              z_alpha, z_beta)$lod
      if (!is.finite(lod_new)) break
      converged <- abs(lod_new - lod) < 1e-6 * max(1, abs(lod))
      lod <- lod_new
      if (converged) break
    }
    # Use sigma at the (final) LOD for LOQ and for the Monte Carlo CI.
    sigma_local <- sigma_fun(lod)
    if (is.finite(sigma_local) && sigma_local > 0) {
      blank_sd <- sigma_local
      loq <- 10 * sigma_local / abs(slope)
    }
  }

  # Monte Carlo CI on the LOD via the calibration parameter distribution.
  lod_ci <- .lod_montecarlo(cal, blank_sd, method, false_positive_rate,
                            false_negative_rate, conf_level, n,
                            hetero, sigma_fun)

  .crossstatr_result(
    list(
      lod = lod,
      loq = loq,
      lod_ci = lod_ci,
      method = method,
      blank_sd = blank_sd,
      slope = slope,
      intercept = intercept,
      residual_se = residual_se,
      df_resid = df_resid,
      heteroscedastic = hetero,
      false_positive_rate = false_positive_rate,
      false_negative_rate = false_negative_rate,
      conf_level = conf_level,
      call = call
    ),
    "crossstatr_lod"
  )
}

# Compute the LOD given the chosen method and a working sigma.
# `blank_sd` carries the working sigma to use for the LOD computation
# (the calibration residual SE in the homoscedastic case, or the local
# lowess-fitted sigma at the LOD in the heteroscedastic case).
# @noRd
.compute_lod <- function(method, slope, intercept, residual_se,
                         blank_sd, x_cal, n, df_resid,
                         alpha, beta, z_alpha, z_beta) {
  sigma <- if (is.finite(blank_sd) && blank_sd > 0) blank_sd else residual_se
  if (method == "currie") {
    # Currie (1968) critical value, general form: (z_alpha + z_beta)*sigma/slope.
    lod <- (z_alpha + z_beta) * sigma / abs(slope)
  } else if (method == "blank_sd") {
    # IUPAC simple form: 3*sigma/slope.
    lod <- 3 * sigma / abs(slope)
  } else {
    # Hubaux-Vos (1970) closed-form fixed-point iteration:
    #   LOD = (t_alpha + t_beta)*sigma/slope * sqrt(1 + 1/n + LOD^2/Sxx)
    lod <- .hubaux_vos_lod(x_cal, slope, sigma, df_resid,
                          alpha, beta)
  }
  list(lod = lod)
}

# Coerce `calibration` (lm or list) to a uniform structure.
# @noRd
.coerce_calibration <- function(calibration) {
  if (inherits(calibration, "lm")) {
    smry <- summary(calibration)
    coef <- stats::coef(calibration)
    if (length(coef) < 2L) {
      stop("Calibration 'lm' must have at least an intercept and slope.",
           call. = FALSE)
    }
    mf <- stats::model.frame(calibration)
    if (ncol(mf) < 2L) {
      stop("Calibration 'lm' needs both y and x.", call. = FALSE)
    }
    list(
      coef = coef[1L:2L],
      residual_se = smry$sigma,
      df_resid = calibration$df.residual,
      x = mf[[2L]],
      y = mf[[1L]]
    )
  } else if (is.list(calibration)) {
    if (is.null(calibration$x) || is.null(calibration$y)) {
      stop("'calibration' list must have 'x' and 'y' components.",
           call. = FALSE)
    }
    x <- calibration$x
    y <- calibration$y
    if (length(x) != length(y)) {
      stop("'calibration' x and y must have equal length.",
           call. = FALSE)
    }
    if (is.null(calibration$coef)) {
      fit <- stats::lm(y ~ x)
      coef <- stats::coef(fit)
      residual_se <- summary(fit)$sigma
      df_resid <- fit$df.residual
    } else {
      coef <- calibration$coef
      res <- y - (coef[1L] + coef[2L] * x)
      residual_se <- sqrt(sum(res^2) / (length(x) - 2L))
      df_resid <- length(x) - 2L
    }
    list(coef = coef[1L:2L], residual_se = residual_se,
         df_resid = df_resid, x = x, y = y)
  } else {
    stop("'calibration' must be an 'lm' object or a list.",
         call. = FALSE)
  }
}

# Hubaux-Vos (1970) prediction-interval LOD via fixed-point iteration.
# Solve  LOD = c * sqrt(1 + 1/n + LOD^2 / Sxx)
# where c = (t_alpha + t_beta) * sigma / |slope|.
# @noRd
.hubaux_vos_lod <- function(x_cal, slope, residual_se, df_resid,
                            alpha, beta, max_iter = 100L, tol = 1e-10) {
  n <- length(x_cal)
  ssx <- sum((x_cal - mean(x_cal))^2)
  if (ssx <= 0) {
    stop("Calibration x has zero variance; cannot compute Hubaux-Vos LOD.",
         call. = FALSE)
  }
  if (df_resid <= 0 || !is.finite(df_resid)) {
    # Fall back to z-based form when residual df unavailable.
    z_alpha <- stats::qnorm(1 - alpha)
    z_beta  <- stats::qnorm(1 - beta)
    return((z_alpha + z_beta) * residual_se / abs(slope))
  }
  t_alpha <- stats::qt(1 - alpha, df = df_resid)
  t_beta  <- stats::qt(1 - beta,  df = df_resid)
  c_const <- (t_alpha + t_beta) * residual_se / abs(slope)
  a_term  <- 1 + 1 / n
  # Initial guess: simple Currie form (ignores LOD^2/Sxx term).
  lod <- c_const
  for (iter in seq_len(max_iter)) {
    lod_new <- c_const * sqrt(a_term + lod^2 / ssx)
    if (abs(lod_new - lod) < tol * max(1, abs(lod_new))) {
      lod <- lod_new
      break
    }
    lod <- lod_new
  }
  lod
}

# Monte Carlo CI on the LOD via simulation from the calibration parameter
# distribution.
# @noRd
.lod_montecarlo <- function(cal, blank_sd, method, alpha, beta,
                            conf_level, n_cal, hetero = FALSE,
                            sigma_fun = NULL, n_sim = 1000L) {
  slope <- cal$coef[2L]
  intercept <- cal$coef[1L]
  # When explicit blank SD is supplied (or heteroscedastic local sigma was
  # used for the point estimate), draw sigma around that central value
  # rather than the calibration residual SE so that the CI reflects the
  # noise actually used to compute the LOD.
  s <- if (is.finite(blank_sd) && blank_sd > 0) blank_sd else cal$residual_se
  df_resid <- cal$df_resid
  V <- .calib_vcov(cal, n_cal)
  z_alpha <- stats::qnorm(1 - alpha)
  z_beta  <- stats::qnorm(1 - beta)
  sim_lod <- numeric(n_sim)
  for (b in seq_len(n_sim)) {
    params <- MASS::mvrnorm(1, mu = c(intercept, slope), Sigma = V)
    s_b <- s * sqrt(stats::rchisq(1, df = df_resid) / df_resid)
    a_b <- params[1L]
    b_b <- params[2L]
    if (abs(b_b) < .Machine$double.eps) {
      sim_lod[b] <- NA_real_
      next
    }
    if (method == "hubaux_vos") {
      sim_lod[b] <- tryCatch(
        .hubaux_vos_lod(cal$x, b_b, s_b, df_resid, alpha, beta),
        error = function(e) NA_real_
      )
    } else if (method == "blank_sd") {
      sim_lod[b] <- 3 * s_b / abs(b_b)
    } else {
      sim_lod[b] <- (z_alpha + z_beta) * s_b / abs(b_b)
    }
  }
  alpha_ci <- 1 - conf_level
  qs <- stats::quantile(sim_lod, probs = c(alpha_ci / 2, 1 - alpha_ci / 2),
                        na.rm = TRUE)
  names(qs) <- NULL
  qs
}

# Approximate variance-covariance matrix of calibration coefficients.
# @noRd
.calib_vcov <- function(cal, n_cal) {
  x <- cal$x
  y <- cal$y
  s <- cal$residual_se
  ssx <- sum((x - mean(x))^2)
  xm <- mean(x)
  V <- matrix(NA_real_, 2L, 2L)
  V[1L, 1L] <- s^2 * (1 / n_cal + xm^2 / ssx)
  V[2L, 2L] <- s^2 / ssx
  V[1L, 2L] <- V[2L, 1L] <- -s^2 * xm / ssx
  V
}

#' @export
print.crossstatr_lod <- function(x, ...) {
  cat(sprintf("<crossStatR::robust_lod> method = %s\n", x$method))
  cat(sprintf("LOD = %.4g  (95%% CI: [%.4g, %.4g])\n",
              x$lod, x$lod_ci[1L], x$lod_ci[2L]))
  cat(sprintf("LOQ = %.4g\n", x$loq))
  cat(sprintf("Slope = %.4g, intercept = %.4g, residual SE = %.4g (df = %d)\n",
              x$slope, x$intercept, x$residual_se, x$df_resid))
  cat(sprintf("Blank SD = %.4g\n", x$blank_sd))
  if (isTRUE(x$heteroscedastic)) {
    cat("Heteroscedastic noise detected: LOD recomputed with local sigma.\n")
  }
  cat(sprintf("False positive rate = %.3f, false negative rate = %.3f\n",
              x$false_positive_rate, x$false_negative_rate))
  invisible(x)
}
