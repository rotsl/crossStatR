#' Drift correction with uncertainty quantification
#'
#' Corrects a signal for instrumental drift by fitting a smooth drift
#' function of time to reference standards and subtracting the fitted
#' drift from the signal, anchored at a reference time. The drift model
#' is one of: linear, polynomial, or first-order exponential approach
#' to asymptote.
#'
#' When `reference = NULL`, the drift model is fit to a running median
#' of the signal itself. This is a less rigorous fallback: the running
#' median conflates drift with genuine slow signal variation, so the
#' correction will be biased if the signal is not approximately
#' constant in the absence of drift. The documentation recommends supplying real
#' reference-standard measurements whenever the experimental design
#' allows it.
#'
#' The correction uncertainty is propagated from the drift model via the
#' delta method. For the linear and polynomial models this is the
#' standard formula for the standard error of a fitted value, and for
#' the exponential model it uses the Jacobian of the prediction
#' function evaluated at the parameter estimates. The resulting SE is
#' on the corrected signal at each time point.
#'
#' @param y numeric vector of signal intensities.
#' @param time numeric vector of measurement times, same length as `y`.
#' @param reference optional data frame with columns `time` and `y`
#'   containing reference-standard measurements taken at known times. If
#'   NULL, the drift model is fit to a running median of `y` (see
#'   Description).
#' @param drift_model character; one of `"linear"`, `"polynomial"`, or
#'   `"exponential"`.
#' @param degree integer; polynomial degree for `drift_model =
#'   "polynomial"`. Default 2.
#' @param conf_level confidence level for the correction uncertainty band.
#'
#' @return A `crossstatr_drift` object with elements:
#'   \item{time, y_original, y_corrected}{vectors after alignment to
#'     `time`.}
#'   \item{drift_fit}{fitted drift values at each time point.}
#'   \item{drift_se}{delta-method standard error of the drift fit at
#'     each time point.}
#'   \item{correction_ci}{n x 2 approximate pointwise uncertainty interval
#'     combining drift-estimation error and measurement noise.}
#'   \item{se_corr, sigma_y}{drift-correction SE and assumed measurement-noise SD.}
#'   \item{model_summary}{list of model-specific summary statistics
#'     (coefficients, R^2 where applicable, residual SE).}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   The drift model is fit to the pairs \eqn{(t_j, r_j)} where \eqn{r_j}
#'   is either the reference-standard reading at time \eqn{t_j} (when
#'   `reference` is supplied) or the running median of `y`. For the
#'   linear model the drift is
#'   \deqn{\hat d(t) = \hat\alpha + \hat\beta\, t,}
#'   fit by OLS. For the polynomial model,
#'   \deqn{\hat d(t) = \sum_{k=0}^{m} \hat\beta_k\, t^k,}
#'   with \eqn{m} = `degree`. For the exponential model,
#'   \deqn{\hat d(t) = \hat a\, (1 - \exp(-\hat b\, t)),}
#'   fit by `stats::nls`.
#'
#'   The corrected signal is
#'   \deqn{\hat y_i^{\mathrm{corr}} = y_i - \hat d(t_i) + \hat d(t_{\mathrm{ref}}),}
#'   where \eqn{t_{\mathrm{ref}}} is the median of the measurement
#'   times. This anchors the corrected signal at \eqn{t_{\mathrm{ref}}}:
#'   the mean of the corrected signal equals the mean of the original
#'   signal at that time.
#'
#'   The delta-method variance of the drift fit at time \eqn{t} is
#'   \deqn{\mathrm{Var}(\hat d(t)) = J(t)^\top\, V_\theta\, J(t),}
#'   where \eqn{V_\theta} is the parameter covariance matrix and
#'   \eqn{J(t)} is the Jacobian of \eqn{\hat d(t)} with respect to the
#'   parameters. The SE on the correction
#'   \eqn{\hat d(t_i) - \hat d(t_{\mathrm{ref}})} additionally accounts
#'   for the covariance between the two fitted values:
#'   \deqn{\mathrm{Var}(\hat d(t_i) - \hat d(t_{\mathrm{ref}})) =
#'         \mathrm{Var}(\hat d(t_i)) + \mathrm{Var}(\hat d(t_{\mathrm{ref}}))
#'         - 2\, \mathrm{Cov}(\hat d(t_i), \hat d(t_{\mathrm{ref}})).}
#'
#' @section Limitations:
#'
#'   - The interval adds measurement-noise variance to drift-correction
#'     variance and uses a Student-t critical value with drift-fit residual df.
#'     Reference-standard residual noise is assumed to represent sample noise.
#'     Without standards, corrected-signal SD supplies a rough noise estimate.
#'     Dependence between a fitted signal and its correction is not propagated;
#'     these intervals are approximate and are not confidence bands for a mean.
#'
#'   - The drift model assumes the drift is a smooth function of time.
#'     If the instrument suffers step changes (e.g. a calibration is
#'     re-applied mid-run), no smooth model will fit well. Inspect the
#'     residuals before trusting the correction.
#'   - The exponential model assumes a monotonic approach to an
#'     asymptote. It will not fit a drift that reverses direction.
#'   - When `reference = NULL`, the drift model is fit to a running
#'     median of the signal, which conflates drift with genuine signal
#'     variation. This is only defensible when the signal is expected
#'     to be approximately constant in the absence of drift.
#'   - The delta-method SE is a local linearisation and may under-cover
#'     for highly nonlinear drift models. For critical applications,
#'     use a parametric bootstrap.
#'
#' @references
#'   - Bates, D. M. and Watts, D. G. (1988). Nonlinear Regression Analysis
#'     and Its Applications. Wiley.
#'     \doi{10.1002/9780470316757}.
#'
#'   - Oehlert, G. W. (1992). A Note on the Delta Method. The American
#'     Statistician, 46, 27-29.
#'     \doi{10.1080/00031305.1992.10475842}.
#'
#' @examples
#' set.seed(1)
#' time <- seq(0, 10, length.out = 100)
#' drift <- 0.5 + 0.1 * time
#' signal <- 5 * exp(-(time - 5)^2 / 2)
#' y <- signal + drift + rnorm(length(time), sd = 0.1)
#' # Reference standards at known times
#' ref_df <- data.frame(
#'   time = seq(0.5, 9.5, by = 1),
#'   y = 0.5 + 0.1 * seq(0.5, 9.5, by = 1) + rnorm(10, sd = 0.05))
#' fit <- drift_correct_stat(y, time, reference = ref_df,
#'                           drift_model = "linear")
#' print(fit)
#'
#' # Fallback to running median when no reference is supplied
#' fit2 <- drift_correct_stat(y, time, drift_model = "exponential")
#' print(fit2)
#'
#' @export
drift_correct_stat <- function(y, time, reference = NULL,
                              drift_model = c("linear", "polynomial",
                                              "exponential"),
                              degree = 2L, conf_level = 0.95) {
  call <- match.call()
  drift_model <- match.arg(drift_model)
  .assert_numeric_vector(y, "y", min_length = 3L)
  .assert_numeric_vector(time, "time", min_length = 3L)
  if (length(y) != length(time)) {
    stop("'y' and 'time' must have the same length.", call. = FALSE)
  }
  .assert_probability(conf_level, "conf_level")
  if (drift_model == "polynomial") {
    .assert_positive_integer(degree, "degree")
  }

  # Build the drift-fit dataset.
  if (!is.null(reference)) {
    .assert_data_frame(reference, "reference", required_cols = c("time", "y"))
    drift_time <- reference$time
    drift_y    <- reference$y
  } else {
    # Running median of y used as a pseudo-reference. This is less
    # rigorous than real reference standards because the running
    # median conflates drift with slow signal variation. Choose an odd
    # window of at most 7 observations.
    k <- max(3L, min(length(y), 7L))
    if (k %% 2L == 0L) k <- k + 1L
    drift_time <- time
    drift_y    <- stats::runmed(y, k)
  }

  n <- length(time)
  t_ref <- stats::median(time, na.rm = TRUE)

  drift_fit_obj <- .fit_drift_model(drift_time, drift_y, drift_model, degree)

  pred_at_time <- .predict_drift(drift_fit_obj, time, drift_model)
  pred_at_ref  <- .predict_drift(drift_fit_obj, t_ref, drift_model)

  drift_at_time <- pred_at_time$fit
  drift_at_ref  <- pred_at_ref$fit

  y_corrected <- y - drift_at_time + as.numeric(drift_at_ref)

  # SE of the drift fit at each time: Var(d(t)) = J' V J.
  drift_se <- .drift_fit_se(drift_fit_obj, time, drift_model)
  # SE of the correction (d(t) - d(t_ref)) including the covariance.
  se_corr <- .drift_correction_se(drift_fit_obj, time, t_ref, drift_model)

  inner_fit <- if (inherits(drift_fit_obj, "lm") || inherits(drift_fit_obj, "nls")) {
    drift_fit_obj
  } else if (is.list(drift_fit_obj) && !is.null(drift_fit_obj$model)) {
    drift_fit_obj$model
  } else {
    NULL
  }
  if (inherits(inner_fit, "lm")) {
    fit_sigma <- summary(inner_fit)$sigma
    if (is.null(reference)) {
      resid_y <- y - pred_at_time$fit + as.numeric(pred_at_ref$fit)
      sigma_y <- if (length(resid_y) > 1L) stats::sd(resid_y) else fit_sigma
    } else {
      sigma_y <- fit_sigma
    }
  } else if (inherits(inner_fit, "nls")) {
    sigma_y <- summary(inner_fit)$sigma
  } else {
    sigma_y <- NA_real_
  }
  se_total <- sqrt(se_corr^2 + if (is.finite(sigma_y)) sigma_y^2 else 0)
  alpha <- 1 - conf_level
  df_resid <- max(1L, stats::df.residual(inner_fit))
  tcrit <- stats::qt(1 - alpha / 2, df = df_resid)
  correction_ci <- cbind(y_corrected - tcrit * se_total,
                         y_corrected + tcrit * se_total)
  colnames(correction_ci) <- paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")

  model_summary <- .drift_model_summary(drift_fit_obj, drift_model)

  .crossstatr_result(
    list(
      time = time,
      y_original = y,
      y_corrected = y_corrected,
      drift_fit = drift_at_time,
      drift_se = drift_se,
      se_corr = se_corr,
      sigma_y = sigma_y,
      correction_ci = correction_ci,
      model_summary = model_summary,
      drift_model = drift_model,
      reference_time = t_ref,
      model = drift_fit_obj,
      conf_level = conf_level,
      call = call
    ),
    "crossstatr_drift"
  )
}

# Fit the drift model.
# @noRd
.fit_drift_model <- function(time, y, drift_model, degree) {
  if (drift_model == "linear") {
    fit <- stats::lm(y ~ time)
    list(model = fit, type = "linear", degree = NA_integer_)
  } else if (drift_model == "polynomial") {
    # Use raw polynomials so the package can construct the design matrix manually
    # for delta-method variance computation.
    fit <- stats::lm(y ~ stats::poly(time, degree = degree, raw = TRUE))
    list(model = fit, type = "polynomial", degree = degree)
  } else {
    # Exponential: y ~ a * (1 - exp(-b * time))
    # Robust starting values: a = dynamic range of y, b = 1 / time_range.
    a0 <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
    if (!is.finite(a0) || a0 <= 0) a0 <- 1
    time_range <- diff(range(time))
    if (!is.finite(time_range) || time_range <= 0) time_range <- 1
    b0 <- 1 / time_range
    form <- stats::as.formula("y ~ a * (1 - exp(-b * time))")
    fit <- tryCatch(
      stats::nls(form, data = data.frame(time = time, y = y),
                 start = list(a = a0, b = b0),
                 control = stats::nls.control(maxiter = 200L)),
      error = function(e) {
        stop("Exponential drift model failed to fit: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    list(model = fit, type = "exponential", degree = NA_integer_)
  }
}

# Predict drift at new times, returning fit and (where possible) SE.
# @noRd
.predict_drift <- function(obj, new_time, drift_model) {
  if (drift_model == "linear") {
    pred <- stats::predict(obj$model, newdata = data.frame(time = new_time),
                          se.fit = TRUE)
    list(fit = as.numeric(pred$fit), se = as.numeric(pred$se.fit))
  } else if (drift_model == "polynomial") {
    # For raw=TRUE polynomials, predict.lm with newdata correctly
    # projects the new points onto the original raw basis.
    pred <- stats::predict(obj$model,
                          newdata = data.frame(time = new_time),
                          se.fit = TRUE)
    list(fit = as.numeric(pred$fit), se = as.numeric(pred$se.fit))
  } else {
    # Exponential: use predict.nls (no se.fit available).
    pred <- stats::predict(obj$model,
                          newdata = data.frame(time = new_time))
    list(fit = as.numeric(pred), se = rep(NA_real_, length(new_time)))
  }
}

# Delta-method SE of the drift fit at each time: Var(d(t)) = J' V J.
# @noRd
.drift_fit_se <- function(obj, time, drift_model) {
  if (drift_model == "linear" || drift_model == "polynomial") {
    X_t <- .design_matrix(obj$model, time, drift_model)
    V <- stats::vcov(obj$model)
    sqrt(pmax(0, rowSums((X_t %*% V) * X_t)))
  } else {
    # Exponential: delta method via analytic Jacobian.
    theta <- stats::coef(obj$model)
    V <- stats::vcov(obj$model)
    J_t <- .exp_jac(theta, time)
    sqrt(pmax(0, rowSums((J_t %*% V) * J_t)))
  }
}

# Delta-method SE on the correction y - d(t) + d(t_ref).
# Var(correction) = Var(d(t)) + Var(d(t_ref)) - 2 Cov(d(t), d(t_ref))
# @noRd
.drift_correction_se <- function(obj, time, t_ref, drift_model) {
  if (drift_model == "linear" || drift_model == "polynomial") {
    X_t <- .design_matrix(obj$model, time, drift_model)
    X_r <- .design_matrix(obj$model, t_ref, drift_model)
    V <- stats::vcov(obj$model)
    var_t <- rowSums((X_t %*% V) * X_t)
    var_r <- as.numeric(X_r %*% V %*% t(X_r))
    cov_tr <- as.numeric(X_t %*% V %*% t(X_r))
    var_corr <- pmax(0, var_t + var_r - 2 * cov_tr)
    sqrt(var_corr)
  } else {
    .exponential_correction_se(obj$model, time, t_ref)
  }
}

# Design matrix for linear / polynomial drift model. Matches the columns
# of stats::poly(time, degree, raw=TRUE) used in the lm fit.
# @noRd
.design_matrix <- function(model, time, drift_model) {
  if (drift_model == "linear") {
    cbind(1, time)
  } else {
    deg <- length(stats::coef(model)) - 1L
    X <- matrix(NA_real_, nrow = length(time), ncol = deg + 1L)
    X[, 1L] <- 1
    for (k in seq_len(deg)) X[, k + 1L] <- time^k
    X
  }
}

# Delta-method SE for exponential drift correction.
# Var(d(t) - d(t_ref)) = J(t)' V J(t) + J(t_ref)' V J(t_ref) - 2 J(t)' V J(t_ref)
# @noRd
.exponential_correction_se <- function(model, time, t_ref) {
  theta <- stats::coef(model)
  V <- stats::vcov(model)
  J_t <- .exp_jac(theta, time)        # n x 2
  J_r <- .exp_jac(theta, t_ref)       # 1 x 2
  var_t <- rowSums((J_t %*% V) * J_t)
  var_r <- as.numeric(J_r %*% V %*% t(J_r))
  cov_tr <- as.numeric(J_t %*% V %*% t(J_r))
  var_corr <- pmax(0, var_t + var_r - 2 * cov_tr)
  sqrt(var_corr)
}

# Jacobian of d(t) = a * (1 - exp(-b*t)) with respect to (a, b).
# Returns n x 2 matrix.
# @noRd
.exp_jac <- function(theta, time) {
  a <- theta[1L]; b <- theta[2L]
  exp_term <- exp(-b * time)
  d_da <- 1 - exp_term
  d_db <- a * time * exp_term
  cbind(d_da, d_db)
}

# Build a small summary list of the drift model.
# @noRd
.drift_model_summary <- function(obj, drift_model) {
  if (drift_model == "linear" || drift_model == "polynomial") {
    smry <- summary(obj$model)
    list(
      coefficients = smry$coefficients,
      r_squared = smry$r.squared,
      adj_r_squared = smry$adj.r.squared,
      residual_se = smry$sigma,
      df_resid = obj$model$df.residual
    )
  } else {
    smry <- summary(obj$model)
    list(
      coefficients = smry$coefficients,
      residual_se = smry$sigma,
      df_resid = smry$df[2L],
      converged = obj$model$convInfo$isConv
    )
  }
}

#' @export
print.crossstatr_drift <- function(x, ...) {
  cat(sprintf("<crossStatR::drift_correct_stat> model = %s\n", x$drift_model))
  cat(sprintf("Reference time: %.4g\n", x$reference_time))
  cat(sprintf("Drift fit range: [%.4g, %.4g]\n",
              min(x$drift_fit, na.rm = TRUE),
              max(x$drift_fit, na.rm = TRUE)))
  cat(sprintf("Corrected signal range: [%.4g, %.4g]\n",
              min(x$y_corrected, na.rm = TRUE),
              max(x$y_corrected, na.rm = TRUE)))
  cat(sprintf("Mean drift SE: %.4g\n",
              mean(x$drift_se, na.rm = TRUE)))
  invisible(x)
}
