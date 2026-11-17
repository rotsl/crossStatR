#' Group comparisons for experimental data
#'
#' Provides the standard battery of two-group comparisons (Student's t,
#' Welch's t, paired t, Wilcoxon rank-sum, Wilcoxon signed-rank) plus a
#' permutation-based alternative for unequal-variance or non-Gaussian
#' settings. Effect sizes (Cohen's d / Hedges' g, rank-biserial for the
#' Wilcoxon tests) and confidence intervals are returned alongside the
#' test statistic.
#'
#' **Test selection is left to the analyst.** This function does not
#' pre-test for normality and pick a test on that basis; doing so inflates
#' Type I error (Rochon et al., 2012).
#'
#' @param x numeric vector (first sample) or a formula referencing a
#'   two-column data.frame.
#' @param y numeric vector (second sample) for two-sample tests; NULL for
#'   paired tests where `y` is implicit, or for formula interface.
#' @param paired logical; use the paired test?
#' @param test character; one of `"auto"` (Welch if n>=30 per group or
#'   balanced, otherwise Student or Wilcoxon — see Details),
#'   `"student"`, `"welch"`, `"perm"` (permutation), `"wilcoxon"`.
#' @param alternative `"two.sided"`, `"less"`, `"greater"`.
#' @param conf_level confidence level for the difference in means and for
#'   the effect size.
#' @param n_perm number of permutations for `test = "perm"` (default 9999).
#' @param seed optional integer seed for permutation reproducibility.
#' @param data optional data frame when `x` is a formula.
#' @param ... additional arguments passed to underlying tests.
#'
#' @return A `crossstatr_compare` object.
#'
#' @section Methodology:
#'
#'   Student's t-statistic (equal variance):
#'   \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{s_p\sqrt{1/n_1 + 1/n_2}}, \qquad
#'         s_p^2 = \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1 + n_2 - 2},}
#'   with \eqn{df = n_1 + n_2 - 2}.
#'
#'   Welch's t-statistic:
#'   \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{s_1^2/n_1 + s_2^2/n_2}}, \qquad
#'         df \approx \frac{(s_1^2/n_1 + s_2^2/n_2)^2}
#'           {(s_1^2/n_1)^2/(n_1-1) + (s_2^2/n_2)^2/(n_2-1)}.}
#'
#'   Paired t-statistic: applied to differences \eqn{d_i = x_i - y_i},
#'   reducing to a one-sample t-test on `d`.
#'
#'   Permutation test: under the null hypothesis of exchangeability,
#'   group labels are randomly permuted `n_perm` times; the two-sided
#'   p-value is
#'   \deqn{p = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\{|T_b| \geq |T_{\mathrm{obs}}|\}}
#'                {B + 1}.}
#'
#'   Wilcoxon rank-sum (Mann-Whitney U) and Wilcoxon signed-rank are
#'   computed via `stats::wilcox.test` with continuity correction.
#'
#'   Effect size Cohen's d / Hedges' g is computed as in
#'   [sci_effect_size()].
#'
#' @section Limitations:
#'
#'   - Standardised-effect intervals are approximate: small samples use a
#'     Student-t critical value when effect-size df is below 30. Hedges g
#'     and its SE use the same small-sample correction factor.
#'
#'   - The permutation test is exact only when the exchangeability
#'     assumption holds (typically true for random assignment, not for
#'     observational data).
#'   - Effect size CIs use the large-sample normal approximation and may
#'     be inaccurate for small samples or extreme effect sizes.
#'   - Paired tests require `length(x) == length(y)`.
#'
#' @references
#'   - Welch, B. L. (1947). The generalization of Student's problem when
#'     several different population variances are involved. Biometrika, 34,
#'     28-35.
#'     \doi{10.1093/biomet/34.1-2.28}.
#'
#'   - Rochon, J., Gondan, M. and Kieser, M. (2012). To test or not to
#'     test: Preliminary assessment of normality when comparing two
#'     independent samples. BMC Medical Research Methodology, 12, 81.
#'     \doi{10.1186/1471-2288-12-81}.
#'
#'   - Davison, A. C. and Hinkley, D. V. (1997). Bootstrap Methods and
#'     their Application. Cambridge University Press.
#'     \doi{10.1017/CBO9780511802843}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(30, mean = 10)
#' y <- rnorm(30, mean = 11)
#' sci_compare(x, y, test = "welch")
#'
#' @export
sci_compare <- function(x, y = NULL, paired = FALSE,
                       test = c("auto", "student", "welch", "perm",
                                "wilcoxon"),
                       alternative = c("two.sided", "less", "greater"),
                       conf_level = 0.95, n_perm = 9999L, seed = NULL,
                       data = NULL, ...) {
  call <- match.call()
  test <- match.arg(test)
  alternative <- match.arg(alternative)
  .assert_probability(conf_level, "conf_level")

  # Formula interface
  if (inherits(x, "formula")) {
    if (is.null(data)) stop("Formula interface requires 'data'.",
                            call. = FALSE)
    mf <- stats::model.frame(x, data = data)
    if (ncol(mf) != 2L) stop("Formula must be of form 'y ~ group'.",
                              call. = FALSE)
    vals <- mf[[1L]]; grp <- mf[[2L]]
    if (length(unique(grp)) != 2L)
      stop("Two-sample comparison requires exactly 2 groups.", call. = FALSE)
    glev <- as.character(sort(unique(grp)))
    xx <- vals[grp == glev[1L]]; yy <- vals[grp == glev[2L]]
  } else {
    .assert_numeric_vector(x, "x", min_length = 1L, finite = FALSE)
    if (is.null(y)) stop("Two-sample tests require 'y'.", call. = FALSE)
    .assert_numeric_vector(y, "y", min_length = 1L, finite = FALSE)
    xx <- x; yy <- y
  }

  if (paired && length(xx) != length(yy)) {
    stop("Paired test requires equal-length samples.", call. = FALSE)
  }

  # Drop NA / Inf
  ok_x <- is.finite(xx); ok_y <- is.finite(yy)
  if (paired) {
    ok <- ok_x & ok_y
    xx <- xx[ok]; yy <- yy[ok]
  } else {
    xx <- xx[ok_x]; yy <- yy[ok_y]
  }
  nx <- length(xx); ny <- length(yy)
  if (nx < 2L || ny < 2L) {
    stop("Each sample needs at least 2 finite values.", call. = FALSE)
  }

  if (test == "auto") {
    # Choose Welch by default (more conservative than Student, robust to unequal var)
    test <- "welch"
  }

  # Compute the test
  if (paired) {
    d <- xx - yy
    if (test == "wilcoxon") {
      tt <- stats::wilcox.test(d, alternative = alternative,
                              conf.level = conf_level,
                              conf.int = TRUE)
      stat <- as.numeric(tt$statistic)
      p <- tt$p.value
      eff <- as.numeric(tt$estimate)  # pseudo-median of differences
      eff_name <- "pseudomedian(d)"
    } else if (test == "perm") {
      obs_t <- mean(d) / (stats::sd(d) / sqrt(length(d)))
      .set_seed(seed)
      perm_t <- replicate(n_perm, {
        ds <- d * sample(c(-1, 1), length(d), replace = TRUE)
        mean(ds) / (stats::sd(ds) / sqrt(length(ds)))
      })
      p <- (1 + sum(abs(perm_t) >= abs(obs_t))) / (n_perm + 1)
      stat <- obs_t
      eff <- mean(d)
      eff_name <- "mean(d)"
    } else {
      tt <- stats::t.test(d, alternative = alternative,
                         conf.level = conf_level)
      stat <- as.numeric(tt$statistic)
      p <- tt$p.value
      ci <- tt$conf.int
      eff <- mean(d)
      eff_name <- "mean(d)"
    }
  } else {
    if (test == "wilcoxon") {
      tt <- stats::wilcox.test(xx, yy, alternative = alternative,
                              conf.level = conf_level,
                              conf.int = TRUE)
      stat <- as.numeric(tt$statistic)
      p <- tt$p.value
      # Rank-biserial correlation (signed-rank summary)
      eff <- as.numeric(tt$estimate)  # Hodges-Lehmann
      eff_name <- "Hodges-Lehmann(d)"
    } else if (test == "perm") {
      obs_diff <- mean(xx) - mean(yy)
      pooled <- c(xx, yy)
      labels <- c(rep(1L, nx), rep(2L, ny))
      .set_seed(seed)
      perm_diff <- replicate(n_perm, {
        s <- sample(labels)
        mean(pooled[s == 1L]) - mean(pooled[s == 2L])
      })
      p <- (1 + sum(abs(perm_diff) >= abs(obs_diff))) / (n_perm + 1)
      stat <- obs_diff
      eff <- obs_diff
      eff_name <- "mean(x) - mean(y)"
    } else if (test == "student") {
      tt <- stats::t.test(xx, yy, alternative = alternative,
                         var.equal = TRUE, conf.level = conf_level)
      stat <- as.numeric(tt$statistic)
      p <- tt$p.value
      eff <- mean(xx) - mean(yy)
      eff_name <- "mean(x) - mean(y)"
    } else { # welch
      tt <- stats::t.test(xx, yy, alternative = alternative,
                         var.equal = FALSE, conf.level = conf_level)
      stat <- as.numeric(tt$statistic)
      p <- tt$p.value
      eff <- mean(xx) - mean(yy)
      eff_name <- "mean(x) - mean(y)"
    }
  }

  # Effect size (Cohen's d, Hedges' g)
  d_val <- .cohens_d(xx, yy, paired = paired, hedges = FALSE)
  g_val <- .cohens_d(xx, yy, paired = paired, hedges = TRUE)
  d_se <- .se_cohens_d(d_val, nx, if (paired) nx else ny, paired = paired)
  df_d <- if (paired) nx - 1L else nx + ny - 2L
  crit <- if (df_d >= 30L) stats::qnorm(1 - (1 - conf_level) / 2)
          else stats::qt(1 - (1 - conf_level) / 2, df = df_d)
  d_ci <- c(d_val - crit * d_se, d_val + crit * d_se)
  g_se <- (1 - 3 / (4 * df_d - 1)) * d_se
  g_ci <- c(g_val - crit * g_se, g_val + crit * g_se)

  out <- .crossstatr_result(
    list(test = test,
         statistic = unname(stat),
         p_value = p,
         effect = setNames(eff, eff_name),
         effect_size = data.frame(
           name = c("Cohen_d", "Hedges_g"),
           estimate = c(d_val, g_val),
           se = c(d_se, g_se),
           ci_low = c(d_ci[1L], g_ci[1L]),
           ci_high = c(d_ci[2L], g_ci[2L])),
         n1 = nx, n2 = if (paired) NA_integer_ else ny,
         alternative = alternative,
         conf_level = conf_level,
         method = sprintf("Comparison: %s%s", test,
                          if (paired) " (paired)" else ""),
         call = call),
    "crossstatr_compare")
  out
}

#' @export
print.crossstatr_compare <- function(x, ...) {
  cat(sprintf("<crossStatR::compare> %s\n", x$method))
  cat(sprintf("statistic = %.4f,  p = %g\n", x$statistic, x$p_value))
  cat(sprintf("effect (%s) = %.4f\n", names(x$effect), as.numeric(x$effect)))
  cat("Effect sizes:\n"); print(x$effect_size, row.names = FALSE)
  invisible(x)
}
