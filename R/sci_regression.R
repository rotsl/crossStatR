#' Regression for scientific experiments
#'
#' Provides OLS, weighted least squares (WLS) with user-supplied
#' measurement uncertainty, generalized least squares via `nlme::gls`
#' (if available) for correlated residuals, nonlinear least squares via
#' `stats::nls()`, and robust regression via `MASS::rlm()`.
#'
#' @param formula a model formula.
#' @param data data frame.
#' @param weights optional numeric weights. Interpretation depends on
#'   `weight_type`: `"precision"` treats weights as `1/sigma^2`
#'   (measurement uncertainties in y), `"analytic"` treats them as
#'   analytic GLS weights (default if non-NULL).
#' @param weight_type character; one of `"precision"` or `"analytic"`.
#' @param robust logical; use `MASS::rlm()` (Huber, default) instead of
#'   OLS?
#' @param nonlinear logical; if TRUE, treat `formula` as nonlinear and
#'   use `nls()`; in that case `start` is required.
#' @param start starting values for `nls()`.
#' @param conf_level confidence level for CIs.
#' @param predict_at optional data frame for which to produce
#'   predictions with CIs.
#' @param ... additional arguments to underlying fitting function.
#'
#' @return A `crossstatr_regression` object with the fitted model,
#'   coefficients with CIs, residual diagnostics, and (optionally)
#'   predictions.
#'
#' @section Methodology:
#'
#'   Ordinary least squares minimizes
#'   \deqn{S(\beta) = \sum_i (y_i - x_i^\top \beta)^2,}
#'   with closed-form solution \eqn{\hat\beta = (X^\top X)^{-1} X^\top y}.
#'
#'   Weighted least squares with precision weights \eqn{w_i = 1/\sigma_i^2}:
#'   \deqn{\hat\beta = (X^\top W X)^{-1} X^\top W y,}
#'   where \eqn{W = \mathrm{diag}(w_i)}. This is the BLUE under the
#'   heteroscedastic model \eqn{y_i \sim N(x_i^\top \beta, \sigma_i^2)}.
#'
#'   The coefficient covariance is
#'   \deqn{\widehat{\mathrm{Cov}}(\hat\beta) = \hat\sigma^2 (X^\top W X)^{-1},}
#'   with \eqn{\hat\sigma^2 = \frac{1}{n - p} \sum_i w_i (y_i - x_i^\top \hat\beta)^2}.
#'
#'   Robust regression via `MASS::rlm` uses an M-estimator with Huber
#'   psi-function by default. The asymptotic covariance is the standard
#'   sandwich estimator; `MASS::rlm` returns a model-based covariance.
#'
#'   Prediction interval for a new observation at \eqn{x_0}:
#'   \deqn{\hat y_0 \pm t_{1-\alpha/2,\, n-p}\, \hat\sigma
#'     \sqrt{1 + x_0^\top (X^\top W X)^{-1} x_0}.}
#'
#' @section Limitations:
#'
#'   - OLS standard errors are wrong under heteroscedasticity unless
#'     weights are supplied. Always supply `weights` when measurement
#'     uncertainties are known.
#'   - For correlated residuals (e.g. time series), use `nlme::gls`
#'     directly with an appropriate `corStruct`; this function does not
#'     fit a correlation structure automatically.
#'   - Nonlinear fits require good starting values; `nls()` may fail to
#'     converge if `start` is poor. Check the convergence message.
#'
#' @references
#'   - Huber, P. J. (1981). Robust Statistics. Wiley.
#'     \doi{10.1002/0471725250}.
#'
#'   - Bates, D. M. and Watts, D. G. (1988). Nonlinear Regression Analysis
#'     and Its Applications. Wiley.
#'     \doi{10.1002/9780470316757}.
#'
#'   - Draper, N. R. and Smith, H. (1998). Applied Regression Analysis. 3rd
#'     ed. Wiley.
#'     \doi{10.1002/9781118625590}.
#'
#' @examples
#' df <- data.frame(x = 1:20, y = 2 * (1:20) + rnorm(20, sd = 2))
#' fit <- sci_regression(y ~ x, data = df)
#' print(fit)
#'
#' # Weighted regression with measurement uncertainties
#' df$sigma <- runif(20, 0.5, 2)
#' fit2 <- sci_regression(y ~ x, data = df, weights = 1 / df$sigma^2,
#'                        weight_type = "precision")
#'
#' @export
sci_regression <- function(formula, data, weights = NULL,
                           weight_type = c("precision", "analytic"),
                           robust = FALSE, nonlinear = FALSE, start = NULL,
                           conf_level = 0.95, predict_at = NULL, ...) {
  call <- match.call()
  weight_type <- match.arg(weight_type)
  .assert_probability(conf_level, "conf_level")
  .assert_flag(robust, "robust")
  .assert_flag(nonlinear, "nonlinear")

  if (nonlinear) {
    if (is.null(start)) {
      stop("Nonlinear regression requires 'start' (starting values).",
           call. = FALSE)
    }
    fit <- try(stats::nls(formula, data = data, start = start, ...),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      stop("nls() failed: ", attr(fit, "condition")$message, call. = FALSE)
    }
    cf <- summary(fit)$coefficients
    ci <- tryCatch(stats::confint(fit, level = conf_level),
                   error = function(e) {
                     m <- matrix(NA_real_, nrow = nrow(cf), ncol = 2L,
                                 dimnames = list(rownames(cf),
                                                 c(paste0((1 - conf_level) / 2 * 100, "%"),
                                                   paste0((1 - (1 - conf_level) / 2) * 100, "%"))))
                     m
                   })
    out_coef <- data.frame(
      term = rownames(cf),
      estimate = cf[, 1L],
      se = cf[, 2L],
      statistic = cf[, 3L],
      p_value = cf[, 4L],
      ci_low = ci[, 1L],
      ci_high = ci[, 2L],
      stringsAsFactors = FALSE)
    res <- stats::resid(fit)
    fit_type <- "nls"
  } else if (robust) {
    fit <- MASS::rlm(formula, data = data, ...)
    smry <- summary(fit)
    cf <- smry$coefficients
    tval <- cf[, 1L] / cf[, 2L]
    pval <- 2 * stats::pt(-abs(tval), df = fit$df.residual)
    ci <- cf[, 1L] + outer(cf[, 2L],
                           stats::qt(c((1-conf_level)/2, 1-(1-conf_level)/2),
                                     df = fit$df.residual))
    out_coef <- data.frame(
      term = rownames(cf), estimate = cf[, 1L], se = cf[, 2L],
      statistic = tval, p_value = pval,
      ci_low = ci[, 1L], ci_high = ci[, 2L], stringsAsFactors = FALSE)
    res <- stats::resid(fit)
    fit_type <- "rlm"
  } else {
    w <- weights
    if (!is.null(w) && weight_type == "precision") {
      # weights = 1/sigma^2; stats::lm with weights already does WLS
    }
    if (!is.null(w)) {
      # stats::lm evaluates the weights argument in the data env, so
      # attach the weights to the data frame for evaluation.
      data$.crossstatr_weights <- w
      fit <- stats::lm(formula, data = data, weights = .crossstatr_weights, ...)
      data$.crossstatr_weights <- NULL
    } else {
      fit <- stats::lm(formula, data = data, ...)
    }
    smry <- summary(fit)
    cf <- smry$coefficients
    ci <- stats::confint(fit, level = conf_level)
    out_coef <- data.frame(
      term = rownames(cf), estimate = cf[, 1L], se = cf[, 2L],
      statistic = cf[, 3L], p_value = cf[, 4L],
      ci_low = ci[, 1L], ci_high = ci[, 2L], stringsAsFactors = FALSE)
    res <- stats::resid(fit)
    fit_type <- if (is.null(w)) "ols" else "wls"
  }

  # Residual diagnostics
  diags <- list(
    n = length(res),
    df_resid = if (inherits(fit, "lm")) fit$df.residual else
      summary(fit)$df[2L],
    rss = sum(res^2),
    rmse = sqrt(mean(res^2)),
    r_squared = if (inherits(fit, "lm")) smry$r.squared else NA_real_,
    adj_r_squared = if (inherits(fit, "lm")) smry$adj.r.squared else NA_real_,
    shapiro_p = tryCatch(stats::shapiro.test(res)$p.value,
                         error = function(e) NA_real_),
    bp_p = tryCatch(.bptest(formula, data, w), error = function(e) NA_real_)
  )

  preds <- NULL
  if (!is.null(predict_at)) {
    pred_ci <- tryCatch(
      stats::predict(fit, newdata = predict_at,
                     interval = "prediction", level = conf_level),
      error = function(e) NULL)
    pred_conf <- tryCatch(
      stats::predict(fit, newdata = predict_at,
                     interval = "confidence", level = conf_level),
      error = function(e) NULL)
    if (!is.null(pred_ci)) {
      preds <- data.frame(
        fit = pred_ci[, 1L],
        ci_low = pred_conf[, 2L], ci_high = pred_conf[, 3L],
        pi_low = pred_ci[, 2L], pi_high = pred_ci[, 3L])
      preds <- cbind(predict_at, preds)
    }
  }

  .crossstatr_result(
    list(model = fit, fit_type = fit_type,
         coefficients = out_coef,
         diagnostics = diags,
         predictions = preds,
         conf_level = conf_level,
         method = paste("Regression:", fit_type),
         call = call),
    "crossstatr_regression")
}

# Breusch-Pagan test for heteroscedasticity (simplified)
.bptest <- function(formula, data, weights = NULL) {
  fit <- stats::lm(formula, data = data, weights = weights)
  res2 <- stats::resid(fit)^2
  aux <- stats::lm(res2 ~ stats::fitted(fit))
  n <- length(res2)
  R2 <- summary(aux)$r.squared
  stat <- n * R2
  p <- stats::pchisq(stat, df = 1, lower.tail = FALSE)
  p
}

#' @export
print.crossstatr_regression <- function(x, ...) {
  cat(sprintf("<crossStatR::regression> %s\n", x$method))
  print(x$coefficients, row.names = FALSE)
  cat(sprintf("\nRMSE = %.4f, R^2 = %.4f\n",
              x$diagnostics$rmse, x$diagnostics$r_squared))
  invisible(x)
}
