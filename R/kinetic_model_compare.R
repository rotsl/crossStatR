#' Compare candidate kinetic models
#'
#' Fits a set of candidate kinetic models to time-concentration data and
#' reports side-by-side AIC, BIC, residual sum of squares, parameter
#' estimates with confidence intervals, and identifiability diagnostics. A
#' best model is selected by AIC, and AIC deltas are reported so the analyst
#' can see whether one model is clearly superior or whether several models
#' are essentially indistinguishable.
#'
#' @param time numeric vector of sampling times.
#' @param conc numeric vector of measured concentrations (same length as
#'   `time`).
#' @param models named list of functions, each taking `(params, t)` and
#'   returning predicted concentrations. If `NULL`, default models are used:
#'   first-order decay, second-order decay, biexponential decay, and
#'   Michaelis-Menten saturation.
#' @param start named list of starting values, one named numeric vector per
#'   model. If `NULL`, data-driven defaults are computed from `time` and
#'   `conc`.
#' @param conf_level confidence level for parameter intervals.
#'
#' @return A `crossstatr_kinetic_compare` object with elements:
#'   \item{comparison}{data frame with one row per model: `model`, `AIC`,
#'     `BIC`, `RSS`, `n_params`, `condition_number`, `identifiable`,
#'     `delta_AIC`.}
#'   \item{best_model}{name of the model with smallest AIC.}
#'   \item{parameter_estimates}{named list of data frames (one per model)
#'     with columns `parameter`, `estimate`, `se`, `ci_low`, `ci_high`.}
#'   \item{models}{the model functions used.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Each model is fit by minimising the residual sum of squares
#'   \deqn{\mathrm{RSS}(\theta) = \sum_{i=1}^{n}(y_i - f(\theta, t_i))^2}
#'   with `stats::nlminb` (box-constrained to non-negative parameters).
#'
#'   AIC and BIC use the Gaussian likelihood form
#'   \deqn{\mathrm{AIC} = n\log(\mathrm{RSS}/n) + 2k, \qquad
#'         \mathrm{BIC} = n\log(\mathrm{RSS}/n) + k\log(n),}
#'   where `k` is the number of parameters. AIC delta is
#'   \eqn{\Delta_i = \mathrm{AIC}_i - \min_j \mathrm{AIC}_j}.
#'
#'   Standard errors for parameters are computed via linearisation. With
#'   sensitivity matrix \eqn{J_{ij} = \partial f_i / \partial \theta_j}
#'   (computed by central finite differences), the parameter
#'   covariance is \eqn{\mathrm{Cov}(\hat{\theta}) \approx \hat{\sigma}^2
#'   (J'J)^{-1}}, with \eqn{\hat{\sigma}^2 = \mathrm{RSS} / (n - k)}.
#'   Confidence intervals use a Student-t critical value with `n - k` df.
#'
#'   Identifiability is diagnosed by the condition number
#'   \eqn{\kappa(J) = \sigma_{\max}(J) / \sigma_{\min}(J)}. When
#'   \eqn{\kappa > 10^8} the parameters are flagged as poorly
#'   identifiable, following the rule of thumb discussed in Seber and Wild
#'   (2003).
#'
#' @section Limitations:
#'
#'   - The default optimiser is `nlminb` with non-negativity bounds. Models
#'     that need signed parameters will require a user-supplied wrapper.
#'   - Standard errors derived from `J'J` are asymptotic and can be
#'     unreliable for small samples or near-flat likelihoods. Profile-likelihood
#'     intervals (`stats::confint` on an `nls` object) are more accurate but
#'     require refitting with `nls`; this function does not do that.
#'   - Likelihood-ratio tests for nested models are not computed
#'     automatically because detecting nesting from a black-box function is
#'     not reliable. The user should inspect AIC deltas and, where nesting
#'     is known, perform an LRT manually using the RSS values.
#'   - Identical or near-identical fits across models are common when the
#'     data are sparse; a small AIC delta (< 2) should be interpreted as
#'     evidence of model non-identifiability, not as evidence for the simpler
#'     model.
#'
#' @references
#'   - Seber, G. A. F. and Wild, C. J. (1989). Nonlinear Regression. Wiley.
#'     \doi{10.1002/0471725315}.
#'
#'   - Akaike, H. (1974). A new look at the statistical model
#'     identification. IEEE Transactions on Automatic Control, 19, 716-723.
#'     \doi{10.1109/TAC.1974.1100705}.
#'
#'   - Bates, D. M. and Watts, D. G. (1988). Nonlinear Regression Analysis
#'     and Its Applications. Wiley.
#'     \doi{10.1002/9780470316757}.
#'
#' @examples
#' set.seed(1)
#' t <- seq(0, 10, length.out = 40)
#' conc <- 5 * exp(-0.3 * t) + rnorm(t, sd = 0.1)
#' result <- kinetic_model_compare(t, conc)
#' print(result)
#'
#' # Parameter estimates for the best model
#' result$parameter_estimates[[result$best_model]]
#'
#' @export
kinetic_model_compare <- function(time, conc, models = NULL, start = NULL,
                                   conf_level = 0.95) {
  call <- match.call()
  .assert_numeric_vector(time, "time", min_length = 4L)
  .assert_numeric_vector(conc, "conc", min_length = 4L)
  if (length(time) != length(conc)) {
    stop("'time' and 'conc' must have equal length.", call. = FALSE)
  }
  .assert_probability(conf_level, "conf_level")

  ok <- is.finite(time) & is.finite(conc)
  if (sum(ok) < 4L) {
    stop("Need at least 4 finite (time, conc) pairs.", call. = FALSE)
  }
  time <- time[ok]; conc <- conc[ok]

  if (is.null(models)) {
    models <- list(
      first_order = function(p, t) p[1L] * exp(-p[2L] * t),
      second_order = function(p, t) 1 / (1 / p[1L] + p[2L] * t),
      biexp = function(p, t) p[1L] * exp(-p[3L] * t) + p[2L] * exp(-p[4L] * t),
      saturation = function(p, t) p[1L] * t / (p[2L] + t)
    )
  }
  if (!is.list(models) || is.null(names(models)) ||
      any(names(models) == "")) {
    stop("'models' must be a named list of functions.", call. = FALSE)
  }
  for (nm in names(models)) {
    if (!is.function(models[[nm]])) {
      stop(sprintf("model '%s' is not a function.", nm), call. = FALSE)
    }
  }

  if (is.null(start)) {
    A0 <- max(conc, na.rm = TRUE)
    t_med <- max(stats::median(time[time > 0]), 1)
    start <- list(
      first_order  = c(A = A0, k = 1 / t_med),
      second_order = c(conc0 = A0, k = 1 / t_med),
      biexp        = c(A = A0 * 0.6, B = A0 * 0.4,
                       k1 = 1 / t_med, k2 = 0.5 / t_med),
      saturation   = c(Vmax = A0 * 1.1, Km = t_med)
    )
  }
  if (!is.list(start)) {
    stop("'start' must be a list of named numeric vectors.", call. = FALSE)
  }

  fit_results <- lapply(names(models), function(nm) {
    st <- start[[nm]]
    if (is.null(st)) {
      stop(sprintf("No starting values for model '%s'.", nm), call. = FALSE)
    }
    .fit_kinetic_model(time, conc, models[[nm]], st, conf_level, nm)
  })
  names(fit_results) <- names(models)

  comp <- do.call(rbind, lapply(names(models), function(nm) {
    fr <- fit_results[[nm]]
    data.frame(
      model = nm,
      AIC = fr$AIC,
      BIC = fr$BIC,
      RSS = fr$RSS,
      n_params = fr$n_params,
      condition_number = fr$condition_number,
      identifiable = fr$identifiable,
      stringsAsFactors = FALSE
    )
  }))
  comp$delta_AIC <- comp$AIC - min(comp$AIC)
  rownames(comp) <- NULL

  best_model <- as.character(comp$model[which.min(comp$AIC)])

  param_est <- lapply(fit_results, function(fr) {
    data.frame(
      parameter = names(fr$params),
      estimate = as.numeric(fr$params),
      se = fr$se,
      ci_low = fr$ci_low,
      ci_high = fr$ci_high,
      stringsAsFactors = FALSE
    )
  })

  .crossstatr_result(
    list(comparison = comp,
         best_model = best_model,
         parameter_estimates = param_est,
         models = models,
         tidy = comp,
         call = call),
    "crossstatr_kinetic_compare"
  )
}

# Fit a single kinetic model and return diagnostics.
# @noRd
.fit_kinetic_model <- function(time, conc, model_func, start, conf_level,
                               model_name) {
  pnames <- names(start)
  start_vec <- unname(start)
  k <- length(start_vec)
  n <- length(time)
  lower <- rep(0, k)

  obj <- function(p) {
    pred <- tryCatch(model_func(p, time), error = function(e) NULL)
    if (is.null(pred) || !is.numeric(pred) || length(pred) != n ||
        any(!is.finite(pred))) {
      return(1e10)
    }
    sum((conc - pred)^2)
  }

  fit <- tryCatch(
    stats::nlminb(start_vec, obj, lower = lower),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$objective) || fit$objective >= 1e9) {
    fit <- tryCatch(
      stats::optim(start_vec, obj, method = "L-BFGS-B", lower = lower),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      stop(sprintf("Failed to fit model '%s'.", model_name), call. = FALSE)
    }
    params <- fit$par
    rss <- fit$value
  } else {
    params <- fit$par
    rss <- fit$objective
  }
  names(params) <- pnames

  sigma2 <- rss / max(1, n - k)
  AIC <- n * log(rss / n) + 2 * k
  BIC <- n * log(rss / n) + k * log(n)

  J <- tryCatch(
    .num_jacobian(params, function(p) model_func(p, time)),
    error = function(e) matrix(NA_real_, nrow = n, ncol = k)
  )
  sv <- tryCatch(svd(J)$d, error = function(e) NA_real_)
  sv <- sv[is.finite(sv) & sv > 0]
  cond <- if (length(sv) >= 2L) max(sv) / min(sv) else NA_real_
  identifiable <- is.finite(cond) && cond < 1e8

  JtJ <- tryCatch(t(J) %*% J, error = function(e) NULL)
  JtJ_inv <- if (!is.null(JtJ)) tryCatch(solve(JtJ), error = function(e) NULL) else NULL
  if (is.null(JtJ_inv)) {
    se <- rep(NA_real_, k)
    ci_lo <- rep(NA_real_, k)
    ci_hi <- rep(NA_real_, k)
  } else {
    cov_par <- JtJ_inv * sigma2
    se <- sqrt(pmax(diag(cov_par), 0))
    df_resid <- max(1, n - k)
    tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = df_resid)
    ci_lo <- as.numeric(params) - tcrit * se
    ci_hi <- as.numeric(params) + tcrit * se
  }

  list(params = params, se = se, ci_low = ci_lo, ci_high = ci_hi,
       RSS = rss, AIC = AIC, BIC = BIC,
       n_params = k, condition_number = cond,
       identifiable = identifiable,
       jacobian = J)
}

#' @export
print.crossstatr_kinetic_compare <- function(x, ...) {
  cat(sprintf("<crossStatR::kinetic_model_compare> best model: %s\n",
              x$best_model))
  cat("\nModel comparison:\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}
