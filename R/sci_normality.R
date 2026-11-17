#' Normality assessment for experimental data
#'
#' Provides a multi-pronged assessment of distributional assumptions.
#' **A non-significant normality test is not evidence of normality**: with
#' small samples the test has very little power to detect deviations, and
#' with large samples it will reject for practically irrelevant
#' deviations. Both regimes are reported explicitly.
#'
#' @param x numeric vector.
#' @param group optional grouping vector; assessments are done per group.
#' @param test character; one of `"shapiro"` (Shapiro-Wilk, default for
#'   `n <= 5000`), `"ad"` (Anderson-Darling, requires the `nortest` package
#'   at runtime via `requireNamespace`), or `"auto"` (picks Shapiro-Wilk
#'   for `n <= 5000`, otherwise suggests Anderson-Darling).
#' @param alpha significance threshold for flagging non-normality
#'   (default 0.05). Used only for the qualitative verdict.
#' @param min_n minimum sample size below which normality testing is
#'   not attempted (default 3 for Shapiro-Wilk).
#' @param plot logical; produce a Q-Q plot when `ggplot2` is available?
#'
#' @return A `crossstatr_normality` object with per-group test results,
#'   sample sizes, distributional moments, skewness and kurtosis, and an
#'   explicit qualitative verdict that distinguishes "not rejected"
#'   from "evidence of normality".
#'
#' @section Methodology:
#'
#'   The Shapiro-Wilk statistic is
#'   \deqn{W = \frac{\bigl(\sum_i a_i x_{(i)}\bigr)^2}
#'                {\sum_i (x_i - \bar{x})^2},}
#'   where `a_i` are tabulated coefficients that depend only on `n`
#'   (Royston, 1992). It is the most powerful general-purpose test for
#'   normality for `n` up to a few thousand.
#'
#'   Anderson-Darling is more sensitive in the tails. Its statistic is
#'   \deqn{A^2 = -n - \frac{1}{n}\sum_{i=1}^{n}(2i-1)
#'     \bigl[\ln z_{(i)} + \ln(1 - z_{(n+1-i)})\bigr],}
#'   with \eqn{z_{(i)} = \Phi((x_{(i)} - \bar{x})/s)}.
#'
#'   Skewness uses the bias-corrected Fisher-Pearson sample estimator
#'   (Doane & Seward, 2011):
#'   \deqn{g_1 = \frac{\sqrt{n(n-1)}}{n-2}\,g_1^\ast,}
#'   with \eqn{g_1^\ast} the biased sample skewness.
#'
#' @section Limitations:
#'
#'   - Normality testing is most useful as a *diagnostic*, not a decision
#'     rule. Pre-test selection of `t` vs. Wilcoxon on the basis of
#'     a Shapiro-Wilk p-value inflates Type I error of the downstream
#'     test (Rochon et al., 2012). Prefer robust methods or specify the
#'     test a priori.
#'   - For very large `n` (>5000), Shapiro-Wilk will reject almost any
#'     real dataset; use the Anderson-Darling test only as a sanity check
#'     and rely on residual diagnostics instead.
#'
#' @references
#'   - Royston, P. (1992). Approximating the Shapiro-Wilk W-test for non-
#'     normality. Statistics and Computing, 2, 117-119.
#'     \doi{10.1007/BF01891203}.
#'
#'   - Stephens, M. A. (1974). EDF Statistics for Goodness of Fit and Some
#'     Comparisons. Journal of the American Statistical Association, 69,
#'     730-737.
#'     \doi{10.1080/01621459.1974.10480196}.
#'
#'   - Doane, D. P. and Seward, L. E. (2011). Measuring Skewness: A
#'     Forgotten Statistic? Journal of Statistics Education, 19.
#'     \doi{10.1080/10691898.2011.11889611}.
#'
#'   - Rochon, J., Gondan, M. and Kieser, M. (2012). To test or not to
#'     test: Preliminary assessment of normality when comparing two
#'     independent samples. BMC Medical Research Methodology, 12, 81.
#'     \doi{10.1186/1471-2288-12-81}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(50)
#' sci_normality(x)
#'
#' @export
sci_normality <- function(x, group = NULL, test = c("shapiro", "ad", "auto"),
                         alpha = 0.05, min_n = 3L, plot = FALSE) {
  call <- match.call()
  test <- match.arg(test)
  .assert_numeric_vector(x, "x", min_length = 1L, finite = FALSE)
  .assert_probability(alpha, "alpha")
  .assert_positive_integer(min_n, "min_n")
  if (length(x) != length(group) && !is.null(group)) {
    stop("Length of 'group' must match 'x'.", call. = FALSE)
  }
  if (is.null(group)) group <- factor(rep(1L, length(x)))
  groups <- as.character(sort(unique(group)))

  res <- lapply(groups, function(g) {
    xi <- x[group == g]
    xi <- xi[is.finite(xi)]
    n <- length(xi)
    if (n < min_n) {
      return(list(group = g, n = n, test = NA_character_,
                 statistic = NA_real_, p_value = NA_real_,
                 skewness = NA_real_, kurtosis = NA_real_,
                 verdict = "too few observations"))
    }
    # Skewness (bias-corrected Fisher-Pearson) and excess kurtosis
    m <- mean(xi); s <- stats::sd(xi)
    if (is.na(s) || s == 0) {
      sk <- NA_real_; ku <- NA_real_
    } else {
      n <- length(xi)
      m2 <- mean((xi - m)^2)
      m3 <- mean((xi - m)^3)
      m4 <- mean((xi - m)^4)
      sk <- (sqrt(n * (n - 1)) / (n - 2)) * m3 / m2^(3 / 2)
      ku <- ((n - 1) * (n + 1) / ((n - 2) * (n - 3))) *
        (m4 / m2^2 - 3 * (n - 1) / (n + 1))
    }
    # Pick the test
    chosen <- if (test == "auto") {
      if (n <= 5000L) "shapiro" else "ad"
    } else test
    if (chosen == "shapiro") {
      if (n > 5000L) {
        warning("Shapiro-Wilk is unreliable for n > 5000; consider Anderson-Darling.",
                call. = FALSE)
      }
      st <- try(stats::shapiro.test(xi), silent = TRUE)
      if (inherits(st, "try-error")) {
        stat <- NA_real_; p <- NA_real_
      } else {
        stat <- as.numeric(st$statistic)
        p <- st$p.value
      }
    } else if (chosen == "ad") {
      if (!requireNamespace("nortest", quietly = TRUE)) {
        warning("Anderson-Darling requires the 'nortest' package; falling back to Shapiro-Wilk.",
                call. = FALSE)
        st <- stats::shapiro.test(xi)
        stat <- as.numeric(st$statistic); p <- st$p.value
        chosen <- "shapiro"
      } else {
        st <- nortest::ad.test(xi)
        stat <- as.numeric(st$statistic); p <- st$p.value
      }
    }
    verdict <- if (is.na(p)) "indeterminate"
    else if (p < alpha) "rejected (non-normal)"
    else if (n < 30) "not rejected (low power)"
    else if (n > 5000) "not rejected (large-sample caveat)"
    else "not rejected"
    list(group = g, n = n, test = chosen,
         statistic = stat, p_value = p,
         skewness = sk, kurtosis = ku,
         verdict = verdict)
  })

  out_df <- do.call(rbind, lapply(res, function(r) {
    data.frame(group = r$group, n = r$n, test = as.character(r$test),
               statistic = r$statistic, p_value = r$p_value,
               skewness = r$skewness, kurtosis = r$kurtosis,
               verdict = r$verdict,
               stringsAsFactors = FALSE)
  }))

  if (plot) {
    .ensure_ggplot()
    qq_df <- do.call(rbind, lapply(groups, function(g) {
      xi <- x[group == g]
      xi <- xi[is.finite(xi)]
      if (length(xi) < 2L) return(NULL)
      qq <- stats::qqnorm(xi, plot.it = FALSE)
      data.frame(group = g, theoretical = qq$x, sample = qq$y)
    }))
    if (!is.null(qq_df)) {
      p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50",
                             linetype = "dashed") +
        ggplot2::facet_wrap(~ group) +
        ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                      title = "Q-Q plot by group") +
        ggplot2::theme_minimal()
      print(p)
    }
  }

  .crossstatr_result(
    list(summary = out_df,
         alpha = alpha,
         method = "Normality assessment",
         call = call),
    "crossstatr_normality")
}

#' @export
print.crossstatr_normality <- function(x, ...) {
  cat(sprintf("<crossStatR::normality> alpha = %.3f\n", x$alpha))
  print(x$summary, row.names = FALSE)
  invisible(x)
}
