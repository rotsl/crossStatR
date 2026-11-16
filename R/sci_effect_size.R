#' Standardized effect sizes with confidence intervals
#'
#' Computes Cohen's d, Hedges' g, paired-d (Cohen's \eqn{d_z}), Pearson r as
#' an effect size, and partial eta-squared, each with a confidence interval
#' and an interpretation against Cohen's (1988) conventions. The function is
#' a single entry point so that downstream code can request an effect size
#' without having to remember which helper computes which.
#'
#' The documentation does not recommend Cohen's benchmarks as anything more than a rough
#' sanity check. They were derived for the behavioral sciences and may be
#' inappropriate in domains where the meaningful size of an effect depends
#' on measurement scale or cost (e.g. electrochemical measurements, where a
#' "small" shift in half-wave potential can be diagnostic). The
#' `interpretation` field is included as a convenience, not as a substitute
#' for domain judgement.
#'
#' @param x numeric vector (first sample).
#' @param y numeric vector (second sample), or `NULL` for one-sample d.
#' @param paired logical; compute a paired effect size? Required for
#'   `type = "cohens_d_z"` and `type = "paired_d"`.
#' @param type character; `"auto"` (default), `"cohens_d"`, `"hedges_g"`,
#'   `"cohens_d_z"`, `"paired_d"`, `"r"`, or `"eta2_partial"`.
#'   `"auto"` dispatches as follows: if `ss_effect` is supplied, partial
#'   eta-squared; if `y` is `NULL`, one-sample Cohen's d; if `paired` is
#'   `TRUE`, paired d (`d_z`); otherwise Cohen's d for two independent
#'   samples.
#' @param conf_level confidence level (default 0.95).
#' @param n1,n2 optional explicit sample sizes used to compute the SE. If
#'   `NULL`, the lengths of `x` and `y` are used.
#' @param ss_effect,ss_error sums of squares for the effect and the error
#'   term. Required when `type = "eta2_partial"`.
#'
#' @return A `crossstatr_effect_size` object with elements:
#'   \item{type}{effect size label.}
#'   \item{estimate}{numeric point estimate.}
#'   \item{se}{standard error (delta method for d-family; Fisher-z for r).}
#'   \item{ci_low, ci_high}{CI endpoints.}
#'   \item{interpretation}{text label: negligible / small / medium / large.}
#'   \item{method}{short method label.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   **Cohen's d** (independent samples, pooled SD):
#'   \deqn{d = \frac{\bar{x}_1 - \bar{x}_2}{s_p}, \qquad
#'         s_p^2 = \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1 + n_2 - 2}.}
#'
#'   **Hedges' g** applies the small-sample bias correction
#'   \deqn{J = 1 - \frac{3}{4(n_1 + n_2 - 2) - 1}, \qquad g = J \cdot d.}
#'
#'   **Paired d (\eqn{d_z})** is the mean of the paired differences divided
#'   by the SD of the differences:
#'   \deqn{d_z = \frac{\bar{d}}{s_d}, \qquad d_i = x_i - y_i.}
#'   This is larger in magnitude than the corresponding independent-samples
#'   `d` when the paired correlation is positive; do not confuse the two
#'   (Cohen, 1988; Morris & DeShon, 2002).
#'
#'   The standard error for the d-family uses the delta-method formula of
#'   Hedges & Olkin (1985, eq. 8.21):
#'   \deqn{\mathrm{SE}(d)^2 \approx \frac{n_1 + n_2}{n_1 n_2} +
#'         \frac{d^2}{2(n_1 + n_2)}}
#'   for independent samples (with a small-sample correction). The CI is then
#'   \eqn{d \pm z_{1-\alpha/2}\, \mathrm{SE}(d)}.
#'
#'   **Pearson r as effect size**: \eqn{r = \mathrm{cor}(x, y)}, with the
#'   CI from the Fisher-z transform, as in [sci_correlation()].
#'
#'   **Partial eta-squared**:
#'   \deqn{\eta_p^2 = \frac{SS_{\mathrm{effect}}}{SS_{\mathrm{effect}} +
#'         SS_{\mathrm{error}}}.}
#'   A CI based on the noncentral F distribution requires the F statistic
#'   and its degrees of freedom. The package does not estimate that here; the package reports a
#'   point estimate and a caveat in `interpretation`.
#'
#'   Cohen's (1988) conventions:
#'   - d: small 0.2, medium 0.5, large 0.8
#'   - r: small 0.1, medium 0.3, large 0.5
#'   - eta²: small 0.01, medium 0.06, large 0.14
#'
#' @section Limitations:
#'
#'   - The d-family SE and CI rely on a large-sample normal approximation
#'     and can have poor coverage for `n < 20` per group or for very large
#'     effects. For small samples, use bootstrap (e.g. [sci_bootstrap()]).
#'   - Partial eta-squared is returned without a CI; do not interpret the
#'     point estimate as if it had one.
#'   - Hedges' g and Cohen's d are computed under an equal-variance
#'     assumption. For unequal variances consider Glass's \eqn{\Delta}
#'     (not implemented here) or a Welch-style standardizer.
#'   - The auto dispatch is intentionally conservative; if the analyst wants a
#'     specific effect size, name it explicitly via `type`.
#'
#' @references
#'   - Cohen, J. (1988). Statistical Power Analysis for the Behavioral
#'     Sciences. 2nd ed. Lawrence Erlbaum Associates.
#'     \doi{10.4324/9780203771587}.
#'
#'   - Hedges, L. V. and Olkin, I. (1985). Statistical Methods for Meta-
#'     Analysis. Academic Press.
#'     \doi{10.1016/C2009-0-03396-0}.
#'
#'   - Morris, S. B. and DeShon, R. P. (2002). Combining effect size
#'     estimates in meta-analysis with repeated measures and independent-
#'     groups designs. Psychological Methods, 7, 105-125.
#'     \doi{10.1037/1082-989X.7.1.105}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(30, mean = 1)
#' y <- rnorm(30, mean = 0)
#' sci_effect_size(x, y, type = "cohens_d")
#'
#' # Paired
#' sci_effect_size(x, y, paired = TRUE, type = "paired_d")
#'
#' # r
#' sci_effect_size(x, y, type = "r")
#'
#' # Partial eta-squared
#' sci_effect_size(ss_effect = 12.5, ss_error = 80, type = "eta2_partial")
#'
#' @export
sci_effect_size <- function(x, y = NULL, paired = FALSE,
                            type = c("auto", "cohens_d", "hedges_g",
                                     "cohens_d_z", "paired_d", "r",
                                     "eta2_partial"),
                            conf_level = 0.95, n1 = NULL, n2 = NULL,
                            ss_effect = NULL, ss_error = NULL) {
  call <- match.call()
  type <- match.arg(type)
  .assert_probability(conf_level, "conf_level")
  .assert_flag(paired, "paired")

  if (type == "auto") {
    if (!is.null(ss_effect)) {
      type <- "eta2_partial"
    } else if (is.null(y)) {
      type <- "cohens_d"
    } else if (paired) {
      type <- "paired_d"
    } else {
      type <- "cohens_d"
    }
  }

  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (type == "eta2_partial") {
    if (is.null(ss_effect) || is.null(ss_error)) {
      stop("type = 'eta2_partial' requires 'ss_effect' and 'ss_error'.",
           call. = FALSE)
    }
    if (!is.numeric(ss_effect) || !is.numeric(ss_error) ||
        ss_effect < 0 || ss_error <= 0) {
      stop("'ss_effect' (>=0) and 'ss_error' (>0) must be numeric.",
           call. = FALSE)
    }
    est <- ss_effect / (ss_effect + ss_error)
    interp <- .interpret_effect(est, "eta2")
    out <- .crossstatr_result(
      list(type = "eta2_partial",
           estimate = est,
           se = NA_real_,
           ci_low = NA_real_,
           ci_high = NA_real_,
           interpretation = interp,
           method = "Partial eta-squared (no CI; noncentral F not computed)",
           call = call),
      "crossstatr_effect_size")
    return(out)
  }

  if (type == "r") {
    if (is.null(y)) stop("type = 'r' requires both 'x' and 'y'.",
                         call. = FALSE)
    .assert_numeric_vector(x, "x", min_length = 3L, finite = FALSE)
    .assert_numeric_vector(y, "y", min_length = 3L, finite = FALSE)
    if (length(x) != length(y)) {
      stop("'x' and 'y' must have the same length for type = 'r'.",
           call. = FALSE)
    }
    ok <- is.finite(x) & is.finite(y)
    xx <- x[ok]; yy <- y[ok]
    n <- length(xx)
    if (n < 3L) stop("Need at least 3 finite pairs for r.", call. = FALSE)
    r_est <- stats::cor(xx, yy)
    se_r <- (1 - r_est^2) / sqrt(n - 1)  # approximate SE on r scale
    if (abs(r_est) < 1) {
      z <- atanh(r_est)
      se_z <- 1 / sqrt(n - 3)
      ci_z <- c(z - zcrit * se_z, z + zcrit * se_z)
      ci <- tanh(ci_z)
    } else {
      ci <- c(r_est, r_est)
    }
    interp <- .interpret_effect(abs(r_est), "r")
    out <- .crossstatr_result(
      list(type = "r",
           estimate = r_est,
           se = se_r,
           ci_low = ci[1L],
           ci_high = ci[2L],
           interpretation = interp,
           method = "Pearson r (Fisher-z CI)",
           call = call),
      "crossstatr_effect_size")
    return(out)
  }

  # d-family (cohens_d, hedges_g, cohens_d_z / paired_d)
  if (type %in% c("cohens_d_z", "paired_d")) {
    paired <- TRUE
    type <- "paired_d"
  }

  if (is.null(y)) {
    .assert_numeric_vector(x, "x", min_length = 2L, finite = FALSE)
    if (is.null(n1)) n1 <- sum(is.finite(x))
    if (paired) {
      stop("paired = TRUE requires 'y'.", call. = FALSE)
    }
    d_val <- .cohens_d(x, hedges = (type == "hedges_g"))
    se_val <- .se_cohens_d(d_val, n1, paired = FALSE)
    method_str <- sprintf("%s (one-sample)",
                          if (type == "hedges_g") "Hedges' g" else "Cohen's d")
  } else {
    .assert_numeric_vector(x, "x", min_length = 2L, finite = FALSE)
    .assert_numeric_vector(y, "y", min_length = 2L, finite = FALSE)
    if (paired) {
      if (length(x) != length(y)) {
        stop("Paired effect size requires equal-length 'x' and 'y'.",
             call. = FALSE)
      }
      if (is.null(n1)) n1 <- sum(is.finite(x) & is.finite(y))
      d_val <- .cohens_d(x, y, paired = TRUE,
                        hedges = (type == "hedges_g"))
      se_val <- .se_cohens_d(d_val, n1, paired = TRUE)
      method_str <- if (type == "hedges_g")
        "Hedges' g (paired)" else "Cohen's d_z (paired)"
    } else {
      if (is.null(n1)) n1 <- sum(is.finite(x))
      if (is.null(n2)) n2 <- sum(is.finite(y))
      d_val <- .cohens_d(x, y, paired = FALSE,
                        hedges = (type == "hedges_g"))
      se_val <- .se_cohens_d(d_val, n1, n2, paired = FALSE)
      method_str <- if (type == "hedges_g")
        "Hedges' g (independent)" else "Cohen's d (independent)"
    }
  }

  ci_lo <- d_val - zcrit * se_val
  ci_hi <- d_val + zcrit * se_val
  interp <- .interpret_effect(abs(d_val), "d")

  out <- .crossstatr_result(
    list(type = type,
         estimate = d_val,
         se = se_val,
         ci_low = ci_lo,
         ci_high = ci_hi,
         interpretation = interp,
         method = method_str,
         n1 = n1,
         n2 = if (is.null(n2)) NA_integer_ else n2,
         paired = paired,
         call = call),
    "crossstatr_effect_size")
  out
}

# Interpret an effect size magnitude against Cohen (1988) conventions.
# @param val magnitude of effect (absolute value).
# @param family one of "d", "r", "eta2".
# @noRd
.interpret_effect <- function(val, family) {
  if (is.na(val)) return(NA_character_)
  cutoffs <- switch(family,
                    d = c(0.2, 0.5, 0.8),
                    r = c(0.1, 0.3, 0.5),
                    eta2 = c(0.01, 0.06, 0.14),
                    c(0.2, 0.5, 0.8))
  labels <- c("negligible", "small", "medium", "large")
  if (val < cutoffs[1L]) labels[1L]
  else if (val < cutoffs[2L]) labels[2L]
  else if (val < cutoffs[3L]) labels[3L]
  else labels[4L]
}

#' @export
print.crossstatr_effect_size <- function(x, ...) {
  cat(sprintf("<crossStatR::effect_size> %s\n", x$method))
  cat(sprintf("estimate = %.4f,  se = %s\n",
              x$estimate,
              if (is.na(x$se)) "n/a" else sprintf("%.4f", x$se)))
  if (!is.na(x$ci_low) && !is.na(x$ci_high)) {
    cat(sprintf("CI: [%.4f, %.4f]\n", x$ci_low, x$ci_high))
  } else {
    cat("CI: not computed\n")
  }
  cat(sprintf("interpretation: %s\n", x$interpretation))
  invisible(x)
}
