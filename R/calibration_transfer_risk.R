#' Statistical framework for calibration transfer risk
#'
#' Quantifies the risk of transferring a calibration model from a source
#' instrument to a target instrument without full recalibration. The
#' function reports the transfer prediction error, an out-of-distribution
#' flag for unknown samples, and an uncertainty interval that accounts
#' for the additional uncertainty introduced by the transfer.
#'
#' @param source_cal an `lm` object from `lm(y ~ x)` on the source
#'   instrument. Required.
#' @param source_data,target_data data frames each with columns `x`
#'   (reference concentration) and `y` (instrument response). Required.
#' @param transfer_samples integer; number of samples shared between
#'   source and target used to estimate the transfer correction
#'   (default 6). If fewer common samples are present in `target_data`,
#'   this is reduced.
#' @param transfer_model `"offset"` (constant shift, default), `"slope_offset"`
#'   (linear correction), or `"none"` (no transfer correction; baseline
#'   for comparison).
#' @param conf_level confidence level for transfer prediction intervals.
#' @param n_boot number of bootstrap replicates for transfer uncertainty
#'   (default 500).
#' @param seed optional integer seed.
#'
#' @return A `crossstatr_cal_transfer_risk` object.
#'
#' @section Methodology:
#'
#'   Let \eqn{f_s(x) = a_s + b_s x} be the source calibration. For a
#'   target instrument the method models
#'   \deqn{y_t = f_t(x) + \epsilon = g(f_s(x)) + \delta + \epsilon,}
#'   where \eqn{g(\cdot)} is a transformation (identity for offset-only,
#'   affine for slope-offset) and \eqn{\delta} captures systematic bias.
#'
#'   The transfer-corrected prediction for a new target response
#'   \eqn{y_t^*} at unknown concentration \eqn{x^*} is obtained by
#'   inverting the source calibration:
#'   \deqn{\hat x^* = f_s^{-1}(g^{-1}(y_t^* - \hat\delta)),}
#'   where \eqn{\hat\delta} (and the slope correction, if any) are
#'   estimated from the shared transfer samples via OLS.
#'
#'   **Transfer prediction error (TPE)** is the standard deviation of
#'   \eqn{\hat x^* - x_{\mathrm{true}}} over the target evaluation set.
#'
#'   **Out-of-distribution detection**: for a new target sample with
#'   response \eqn{y_t^*}, the package computes the standardized transfer residual
#'   \deqn{r^* = (y_t^* - \hat\delta - g(\hat x^*))/\hat\sigma_t.}
#'   If \eqn{|r^*| > z_{1 - \alpha/2}}, the sample is flagged as
#'   out-of-distribution (matrix effect, sensor drift, or different
#'   analyte).
#'
#'   **Transfer uncertainty** combines the source calibration
#'   uncertainty, the transfer-correction estimator uncertainty, and the
#'   target measurement noise via a parametric bootstrap: simulate
#'   source calibration parameters from their MVN distribution, simulate
#'   transfer-correction parameters from their MVN, and recompute the
#'   inverse prediction; the resulting distribution of \eqn{\hat x^*}
#'   captures the joint uncertainty.
#'
#'   **Reliability summary**: a transfer is deemed acceptable if TPE is
#'   below a user-supplied tolerance (default: 2x the source calibration
#'   standard error) and the OOD rate on the evaluation set is below 5%.
#'
#' @section Limitations:
#'
#'   - Assumes the source calibration is correctly specified. If the
#'     source calibration itself is biased, the transfer correction
#'     cannot fix it.
#'   - The OOD detection is univariate; multivariate OOD (e.g. on
#'     full spectra) requires a separate approach.
#'   - The bootstrap is parametric; nonparametric bootstrap is an
#'     alternative when transfer samples are plentiful (n > 30).
#'   - For drift between source and target acquisition times, use
#'     [drift_uncertainty()] jointly with this function.
#'
#' @references
#'   - Wang, Y., Veltkamp, D. J. and Kowalski, B. R. (1991). Multivariate
#'     instrument standardization. Analytical Chemistry, 63, 2750-2756.
#'     \doi{10.1021/ac00023a016}.
#'
#'   - Bouveresse, E. and Massart, D. L. (1996). Standardisation of near-
#'     infrared spectrometric instruments: A review. Vibrational
#'     Spectroscopy, 11, 3-15.
#'     \doi{10.1016/0924-2031(95)00055-0}.
#'
#'   - Feudale, R. N., Woody, N. A., Tan, H., Myles, A. J., Brown, S. D.
#'     and Ferré, J. (2002). Transfer of multivariate calibration models: a
#'     review. Chemometrics and Intelligent Laboratory Systems, 64,
#'     181-192.
#'     \doi{10.1016/S0169-7439(02)00085-0}.
#'
#' @examples
#' set.seed(1)
#' # Source calibration data
#' x <- seq(0.1, 10, length.out = 30)
#' y_s <- 2 + 3 * x + rnorm(length(x), sd = 0.2)
#' source_cal <- lm(y_s ~ x)
#'
#' # Target: same samples with offset and slope change
#' y_t <- 2.5 + 2.7 * x + rnorm(length(x), sd = 0.25)
#' source_data <- data.frame(x = x, y = y_s)
#' target_data <- data.frame(x = x, y = y_t)
#'
#' res <- calibration_transfer_risk(source_cal, source_data, target_data)
#' print(res)
#'
#' @export
calibration_transfer_risk <- function(source_cal, source_data, target_data,
                                       transfer_samples = 6L,
                                       transfer_model = c("offset", "slope_offset", "none"),
                                       conf_level = 0.95,
                                       n_boot = 500L, seed = NULL) {
  call <- match.call()
  transfer_model <- match.arg(transfer_model)
  .assert_probability(conf_level, "conf_level")

  if (!inherits(source_cal, "lm")) {
    stop("'source_cal' must be an lm object from lm(y ~ x).", call. = FALSE)
  }
  .assert_data_frame(source_data, "source_data", required_cols = c("x", "y"))
  .assert_data_frame(target_data, "target_data", required_cols = c("x", "y"))

  # Source calibration parameters and covariance
  beta_s <- stats::coef(source_cal)
  Vbeta_s <- stats::vcov(source_cal)
  sigma_s <- summary(source_cal)$sigma
  df_s <- source_cal$df.residual

  # Align source and target on common x
  common_x <- intersect(source_data$x, target_data$x)
  if (length(common_x) < 3L) {
    stop("Source and target must share at least 3 common concentrations.",
         call. = FALSE)
  }
  n_transfer <- min(transfer_samples, length(common_x))
  # Pick transfer samples
  set.seed(seed %||% 1L)
  transfer_x <- sample(common_x, n_transfer)
  s_row <- match(transfer_x, source_data$x)
  t_row <- match(transfer_x, target_data$x)
  y_s_t <- source_data$y[s_row]
  y_t_t <- target_data$y[t_row]

  # Estimate transfer correction: y_t = g(y_s) + delta
  # offset only: delta = mean(y_t_t - y_s_t)
  # slope_offset: y_t = alpha + beta * y_s
  if (transfer_model == "offset") {
    delta <- mean(y_t_t - y_s_t)
    se_delta <- stats::sd(y_t_t - y_s_t) / sqrt(n_transfer)
    alpha_t <- delta; beta_t <- 1
    V_alpha_beta <- matrix(c(se_delta^2, 0, 0, 0), 2L, 2L)
  } else if (transfer_model == "slope_offset") {
    fit_t <- stats::lm(y_t_t ~ y_s_t)
    alpha_t <- stats::coef(fit_t)[1L]
    beta_t  <- stats::coef(fit_t)[2L]
    V_alpha_beta <- stats::vcov(fit_t)
    se_delta <- NA_real_
  } else {
    alpha_t <- 0; beta_t <- 1
    V_alpha_beta <- matrix(0, 2L, 2L)
    se_delta <- 0
  }

  # Predict x from target responses on evaluation set (all target points)
  eval_y_t <- target_data$y
  # Invert: y_s_corrected = (y_t - alpha) / beta; then x = (y_s_corr - a_s) / b_s
  y_s_corr <- (eval_y_t - alpha_t) / beta_t
  x_pred <- (y_s_corr - beta_s[1L]) / beta_s[2L]
  true_x <- target_data$x
  tpe <- sqrt(mean((x_pred - true_x)^2))

  # Out-of-distribution: standardized transfer residual
  pred_y_from_x <- alpha_t + beta_t * (beta_s[1L] + beta_s[2L] * true_x)
  resid_t <- eval_y_t - pred_y_from_x
  sigma_t <- stats::sd(resid_t)
  if (!is.finite(sigma_t) || sigma_t == 0) sigma_t <- sigma_s
  std_resid <- resid_t / sigma_t
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ood_flag <- abs(std_resid) > z_crit
  ood_rate <- mean(ood_flag)

  # Bootstrap transfer uncertainty
  .set_seed(seed)
  boot_x_pred <- replicate(n_boot, {
    # Sample source calibration parameters
    bs_beta <- MASS::mvrnorm(1L, mu = beta_s, Sigma = Vbeta_s)
    # Sample transfer-correction parameters
    if (transfer_model == "offset") {
      b_alpha <- rnorm(1L, mean = alpha_t, sd = se_delta)
      b_beta <- 1
    } else if (transfer_model == "slope_offset") {
      ab <- MASS::mvrnorm(1L, mu = c(alpha_t, beta_t), Sigma = V_alpha_beta)
      b_alpha <- ab[1L]; b_beta <- ab[2L]
    } else {
      b_alpha <- 0; b_beta <- 1
    }
    y_s_c <- (eval_y_t - b_alpha) / b_beta
    (y_s_c - bs_beta[1L]) / bs_beta[2L]
  })
  boot_x_pred <- as.vector(boot_x_pred)
  # Prediction interval per sample: take percentiles per row
  pi_mat <- matrix(boot_x_pred, nrow = length(eval_y_t), byrow = FALSE)
  alpha <- 1 - conf_level
  ci_low  <- apply(pi_mat, 1L, stats::quantile, probs = alpha / 2)
  ci_high <- apply(pi_mat, 1L, stats::quantile, probs = 1 - alpha / 2)

  # Reliability summary
  source_cal_se <- summary(source_cal)$sigma / abs(beta_s[2L])
  acceptable_tpe_threshold <- 2 * source_cal_se
  acceptable_ood_threshold <- 0.05
  acceptable <- tpe <= acceptable_tpe_threshold && ood_rate <= acceptable_ood_threshold

  reliability <- list(
    tpe = tpe,
    tpe_threshold = acceptable_tpe_threshold,
    ood_rate = ood_rate,
    ood_threshold = acceptable_ood_threshold,
    acceptable = acceptable,
    n_transfer = n_transfer,
    n_eval = nrow(target_data)
  )

  .crossstatr_result(
    list(predictions = data.frame(
            true_x = true_x, predicted_x = x_pred,
            ci_low = ci_low, ci_high = ci_high,
            resid = resid_t, std_resid = std_resid,
            ood = ood_flag),
         transfer_correction = list(
           model = transfer_model,
           alpha = alpha_t, beta = beta_t,
           se_alpha = if (transfer_model != "none") sqrt(V_alpha_beta[1L, 1L]) else 0,
           se_beta  = if (transfer_model == "slope_offset") sqrt(V_alpha_beta[2L, 2L]) else NA_real_
         ),
         reliability = reliability,
         source_cal = list(coef = beta_s, vcov = Vbeta_s, sigma = sigma_s,
                           df = df_s),
         n_boot = n_boot,
         conf_level = conf_level,
         method = sprintf("Calibration transfer (%s correction, %d transfer samples)",
                          transfer_model, n_transfer),
         call = call),
    "crossstatr_cal_transfer_risk")
}

#' @export
print.crossstatr_cal_transfer_risk <- function(x, ...) {
  cat(sprintf("<crossStatR::calibration_transfer_risk> %s\n", x$method))
  cat(sprintf("TPE = %.4g (acceptable: < %.4g): %s\n",
              x$reliability$tpe, x$reliability$tpe_threshold,
              if (x$reliability$tpe <= x$reliability$tpe_threshold) "PASS" else "FAIL"))
  cat(sprintf("OOD rate = %.0f%% (acceptable: < %.0f%%): %s\n",
              100 * x$reliability$ood_rate,
              100 * x$reliability$ood_threshold,
              if (x$reliability$ood_rate <= x$reliability$ood_threshold) "PASS" else "FAIL"))
  cat(sprintf("Overall transfer acceptable: %s\n",
              if (x$reliability$acceptable) "YES" else "NO"))
  invisible(x)
}
