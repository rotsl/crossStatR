#' Compare two calibration curves via ANCOVA
#'
#' Compares two linear calibration models (typically linear regressions of
#' response on concentration) to test whether they share the same slope
#' and the same intercept. The comparison uses a single ANCOVA fit on
#' the stacked data with an interaction between concentration and a
#' model indicator, then tests the interaction (slopes equal?) and the
#' main effect of the indicator (intercepts equal given equal slopes).
#'
#' Both models can be supplied either as `lm` objects (from
#' `stats::lm(y ~ x)`) or as data frames with columns `x` and `y`. If a
#' data frame is supplied, the linear model is fit internally by OLS.
#'
#' @param model1,model2 `lm` objects (from `stats::lm(y ~ x)`) or data
#'   frames with columns `x` and `y`.
#' @param conf_level confidence level for the slope-difference and
#'   intercept-difference CIs.
#' @param method character; currently only `"ancova"` is implemented.
#'
#' @return A `crossstatr_calibration_compare` object with elements:
#'   \item{summary}{data frame with one row per comparison (slope
#'     difference, intercept difference), each with estimate, SE, CI, and
#'     p-value.}
#'   \item{ancova_fit}{the combined `lm` object from the ANCOVA fit.}
#'   \item{slope_equal}{logical; TRUE if slopes are not significantly
#'     different at the chosen `conf_level`.}
#'   \item{slope_equal_p}{p-value for the slope-equality test.}
#'   \item{intercept_equal}{logical; TRUE if intercepts are not
#'     significantly different.}
#'   \item{intercept_equal_p}{p-value for the intercept-equality test.}
#'   \item{slope1, slope2, intercept1, intercept2}{the fitted slopes and
#'     intercepts of each calibration.}
#'   \item{residual_se1, residual_se2}{residual SE of each model.}
#'   \item{n1, n2}{number of observations in each calibration.}
#'   \item{f_statistic, f_p_value}{F statistic and p-value for the
#'     overall model-comparison test.}
#'   \item{method}{backend used.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Let the two calibrations be
#'   \deqn{y_{1i} = \alpha_1 + \beta_1 x_{1i} + \epsilon_{1i}, \qquad
#'         y_{2j} = \alpha_2 + \beta_2 x_{2j} + \epsilon_{2j},}
#'   with \eqn{\epsilon \sim N(0, \sigma^2)} assumed common. Define an
#'   indicator \eqn{g \in \{0, 1\}} for the second calibration. The
#'   combined ANCOVA model is
#'   \deqn{y = \alpha + \beta\, x + \delta\, g + \gamma\, (x \cdot g) +
#'         \epsilon,}
#'   where \eqn{\delta = \alpha_2 - \alpha_1} and \eqn{\gamma = \beta_2 -
#'   \beta_1}. The slope-difference test is \eqn{H_0: \gamma = 0}, and the
#'   intercept-difference test (given equal slopes) is \eqn{H_0: \delta =
#'   0} in the reduced model without interaction.
#'
#'   The CIs on the slope and intercept differences are read directly off
#'   the ANCOVA fit coefficients. The overall F statistic tests the joint
#'   null \eqn{H_0: \delta = \gamma = 0} via
#'   \deqn{F = \frac{(\mathrm{RSS}_{\mathrm{reduced}} -
#'         \mathrm{RSS}_{\mathrm{full}}) / 2}
#'        {\mathrm{RSS}_{\mathrm{full}} / (n_1 + n_2 - 4)},}
#'   with \eqn{2} and \eqn{n_1 + n_2 - 4} degrees of freedom under the
#'   null.
#'
#' @section Limitations:
#'
#'   - The ANCOVA assumes equal residual variance across the two
#'     calibrations. If the residual SDs differ by more than a factor of
#'     2, the F test is anti-conservative. Inspect `residual_se1` and
#'     `residual_se2`.
#'   - The test assumes linear calibrations. If either calibration is
#'     visibly nonlinear (e.g. a quadratic term is significant), the
#'     slope-comparison test is not meaningful.
#'   - The CIs on slope and intercept differences assume the two
#'     calibrations are independent. If they share samples, instruments,
#'     or operators, this is wrong; the CIs will be too narrow.
#'   - The combined ANCOVA fit uses the same residual variance for both
#'     calibrations. If the data clearly show heteroscedasticity between
#'     groups, use `nlme::gls` with a `varIdent` variance structure
#'     instead.
#'
#' @references
#'   - Draper, N. R. and Smith, H. (1998). Applied Regression Analysis. 3rd
#'     ed. Wiley.
#'     \doi{10.1002/9781118625590}.
#'
#' @examples
#' set.seed(1)
#' x1 <- seq(0, 10, length.out = 15)
#' x2 <- seq(0, 10, length.out = 12)
#' y1 <- 0.5 + 2.0 * x1 + rnorm(length(x1), sd = 0.3)
#' y2 <- 0.7 + 1.9 * x2 + rnorm(length(x2), sd = 0.3)
#' m1 <- lm(y1 ~ x1)
#' m2 <- lm(y2 ~ x2)
#' cmp <- calibration_compare(m1, m2)
#' print(cmp)
#'
#' # Different data frame inputs
#' df1 <- data.frame(x = x1, y = y1)
#' df2 <- data.frame(x = x2, y = y2)
#' cmp2 <- calibration_compare(df1, df2)
#' print(cmp2)
#'
#' @export
calibration_compare <- function(model1, model2,
                              conf_level = 0.95,
                              method = c("ancova")) {
  call <- match.call()
  method <- match.arg(method)
  .assert_probability(conf_level, "conf_level")

  fit1 <- .coerce_to_lm(model1, name = "model1")
  fit2 <- .coerce_to_lm(model2, name = "model2")

  coef1 <- stats::coef(fit1)
  coef2 <- stats::coef(fit2)
  slope1 <- coef1[2L]; intercept1 <- coef1[1L]
  slope2 <- coef2[2L]; intercept2 <- coef2[1L]
  smry1 <- summary(fit1); smry2 <- summary(fit2)
  resid_se1 <- smry1$sigma
  resid_se2 <- smry2$sigma
  n1 <- length(stats::fitted(fit1))
  n2 <- length(stats::fitted(fit2))
  alpha <- 1 - conf_level

  # Combined ANCOVA: stack the data, fit y ~ x * group.
  d1 <- data.frame(x = stats::model.frame(fit1)[[2L]],
                   y = stats::model.frame(fit1)[[1L]],
                   group = 0)
  d2 <- data.frame(x = stats::model.frame(fit2)[[2L]],
                   y = stats::model.frame(fit2)[[1L]],
                   group = 1)
  combined <- rbind(d1, d2)

  ancova_fit <- stats::lm(y ~ x * group, data = combined)
  ancova_smry <- summary(ancova_fit)
  cf <- ancova_smry$coefficients

  # Row order in cf: (Intercept), x, group, x:group.
  # The "group" coefficient is the intercept difference (intercept2 - intercept1).
  # The "x:group" coefficient is the slope difference (slope2 - slope1).
  get_row <- function(nm, default = NA_real_) {
    if (nm %in% rownames(cf)) cf[nm, ] else rep(default, 4L)
  }
  slope_row <- get_row("x:group")
  int_row   <- get_row("group")
  slope_diff     <- slope_row[1L]
  slope_diff_se  <- slope_row[2L]
  slope_diff_p   <- slope_row[4L]
  int_diff       <- int_row[1L]
  int_diff_se    <- int_row[2L]
  int_diff_p     <- int_row[4L]

  # Overall F test: compare the full ANCOVA (separate slope and intercept)
  # against the reduced common-line model y ~ x.
  reduced_fit <- stats::lm(y ~ x, data = combined)
  an <- stats::anova(reduced_fit, ancova_fit)
  f_stat <- an$F[2L]
  f_p    <- an$`Pr(>F)`[2L]

  z_or_t <- stats::qt(1 - alpha / 2, df = n1 + n2 - 4L)
  slope_diff_ci <- slope_diff + c(-1, 1) * z_or_t * slope_diff_se
  int_diff_ci  <- int_diff  + c(-1, 1) * z_or_t * int_diff_se

  summary_df <- data.frame(
    comparison = c("slope_diff", "intercept_diff"),
    estimate = c(slope_diff, int_diff),
    se = c(slope_diff_se, int_diff_se),
    ci_low = c(slope_diff_ci[1L], int_diff_ci[1L]),
    ci_high = c(slope_diff_ci[2L], int_diff_ci[2L]),
    p_value = c(slope_diff_p, int_diff_p),
    stringsAsFactors = FALSE
  )

  slope_equal <- !is.na(slope_diff_p) && slope_diff_p > alpha
  intercept_equal <- !is.na(int_diff_p) && int_diff_p > alpha

  .crossstatr_result(
    list(
      summary = summary_df,
      ancova_fit = ancova_fit,
      slope_equal = slope_equal,
      slope_equal_p = slope_diff_p,
      intercept_equal = intercept_equal,
      intercept_equal_p = int_diff_p,
      slope1 = slope1,
      slope2 = slope2,
      intercept1 = intercept1,
      intercept2 = intercept2,
      slope_diff = slope_diff,
      intercept_diff = int_diff,
      residual_se1 = resid_se1,
      residual_se2 = resid_se2,
      n1 = n1,
      n2 = n2,
      f_statistic = f_stat,
      f_p_value = f_p,
      method = method,
      conf_level = conf_level,
      call = call
    ),
    "crossstatr_calibration_compare"
  )
}

# Coerce a model1/model2 argument to an lm fit.
# Accepts: an lm object or a data.frame with columns x and y.
# @noRd
.coerce_to_lm <- function(model, name = "model") {
  if (inherits(model, "lm")) {
    return(model)
  }
  if (is.data.frame(model)) {
    if (!("x" %in% names(model)) || !("y" %in% names(model))) {
      stop(sprintf("'%s' data.frame must have columns 'x' and 'y'.", name),
           call. = FALSE)
    }
    return(stats::lm(y ~ x, data = model))
  }
  stop(sprintf("'%s' must be an lm object or a data.frame.", name),
       call. = FALSE)
}

#' @export
print.crossstatr_calibration_compare <- function(x, ...) {
  cat(sprintf("<crossStatR::calibration_compare> method = %s\n", x$method))
  cat(sprintf("n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("Slope1 = %.4g (SE_resid = %.4g), Slope2 = %.4g (SE_resid = %.4g)\n",
              x$slope1, x$residual_se1, x$slope2, x$residual_se2))
  cat(sprintf("Intercept1 = %.4g, Intercept2 = %.4g\n",
              x$intercept1, x$intercept2))
  cat(sprintf("Slopes equal: %s (p = %.4g), Intercepts equal: %s (p = %.4g)\n",
              x$slope_equal, x$slope_equal_p,
              x$intercept_equal, x$intercept_equal_p))
  cat(sprintf("Overall F = %.4g (p = %.4g)\n", x$f_statistic, x$f_p_value))
  cat("\nComparison summary:\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
