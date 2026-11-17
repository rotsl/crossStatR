#' Test for differences across complete spectra
#'
#' Tests whether complete spectra (rows of a matrix) differ across groups.
#' Three complementary methods are offered: a Hotelling T-squared test for
#' two groups (multivariate analogue of the two-sample t-test), a PERMANOVA
#' F-test based on a Euclidean distance matrix with permutation-derived
#' p-values, and a battery of pointwise t-tests with multiplicity
#' adjustment.
#'
#' The Hotelling test is the right tool when the dimensionality is small
#' relative to the sample size and the analyst wantss an omnibus statement.
#' PERMANOVA is preferred when the dimensionality exceeds the sample size
#' or when the test must generalise to more than two groups. Pointwise tests
#' answer a different question (where do the spectra differ?) and should be
#' treated as exploratory follow-up rather than as a primary omnibus test.
#'
#' @param spectra numeric matrix or data frame; rows are observations,
#'   columns are wavelengths (or m/z channels, time bins, etc.).
#' @param group factor or coercible vector of group memberships, length
#'   equal to `nrow(spectra)`. Must have at least 2 levels.
#' @param wavelengths optional numeric vector of wavelength values, length
#'   equal to `ncol(spectra)`. If `NULL`, column indices are used.
#' @param method `"t2"`, `"permanova"`, or `"pointwise_t"`.
#' @param conf_level confidence level (used only for some auxiliary intervals).
#' @param n_perm integer; number of permutations for PERMANOVA.
#' @param seed optional integer seed for permutation reproducibility.
#' @param adjust p-value adjustment method for the pointwise tests; any
#'   method accepted by [stats::p.adjust()].
#'
#' @return A `crossstatr_spectral_change` object with elements:
#'   \item{test_statistic}{numeric scalar; the F-statistic, T-squared, or
#'     max abs t-statistic (depending on `method`).}
#'   \item{p_value}{numeric scalar p-value.}
#'   \item{method}{method label.}
#'   \item{n_perm}{number of permutations used (NA for non-permutation methods).}
#'   \item{pointwise}{for `pointwise_t`, a data frame of `wavelength`,
#'     `t`, `p`, `p_adj`, `significant`; otherwise `NULL`.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Hotelling T-squared (two-group case):
#'   \deqn{T^2 = \frac{n_1 n_2}{n_1 + n_2}(\bar{y}_1 - \bar{y}_2)' S_p^{-1}(\bar{y}_1 - \bar{y}_2),}
#'   where \eqn{S_p} is the pooled within-group covariance. When
#'   \eqn{p < n_1 + n_2 - 1} and \eqn{S_p} is invertible, the F-transform
#'   \deqn{F = \frac{n_1 + n_2 - p - 1}{p(n_1 + n_2 - 2)}T^2}
#'   is referenced to \eqn{F_{p,\,n_1+n_2-p-1}} (Hotelling, 1931).
#'
#'   When the pooled covariance is singular or the F-transform has
#'   non-positive denominator df, the statistic is descriptive only.
#'   A warning and `p_value = NA` are returned; `method = "permanova"`
#'   supplies a separate permutation-based comparison.
#'
#'   PERMANOVA pseudo-F (Anderson, 2001). Let \eqn{d_{ij}} be the Euclidean
#'   distance between rows `i` and `j` of `spectra`. The total sum of
#'   squares is
#'   \deqn{\mathrm{SS}_T = \frac{1}{n}\sum_{i=1}^{n}\sum_{j=1}^{n}d_{ij}^2,}
#'   the within-group sum of squares is
#'   \deqn{\mathrm{SS}_W = \sum_g \frac{1}{n_g}\sum_{i,j \in g}d_{ij}^2,}
#'   and the between-group sum of squares is
#'   \eqn{\mathrm{SS}_A = \mathrm{SS}_T - \mathrm{SS}_W}. The pseudo-F is
#'   \deqn{F = \frac{\mathrm{SS}_A / (g - 1)}{\mathrm{SS}_W / (n - g)}.}
#'   Permutation of group labels gives the empirical p-value
#'   \deqn{p = \frac{1 + \#\{F_b \geq F_{\mathrm{obs}}\}}{B + 1}.}
#'
#'   Pointwise two-sample t-tests are run at each wavelength, with
#'   p-values adjusted using `stats::p.adjust` (default `"fdr"`, Benjamini
#'   and Hochberg, 1995). The reported `test_statistic` is the maximum
#'   absolute t-statistic; the reported `p_value` is the minimum adjusted
#'   p-value, which provides a global test with strong control under
#'   independence.
#'
#' @section Limitations:
#'
#'   - The Hotelling F-transform requires \eqn{p < n_1 + n_2 - 1}. For
#'     wider spectra or singular pooled covariance, no Hotelling p-value is
#'     returned. PERMANOVA requires its own exchangeability assumptions.
#'   - PERMANOVA assumes exchangeability of observations under the null.
#'     When groups have very different dispersions (not just means) the
#'     test can be significant even without location differences; see
#'     Anderson (2006) and Warton et al. (2012).
#'   - The pointwise tests ignore the multivariate dependence across
#'     wavelengths. The minimum adjusted p-value should be read as
#'     exploratory rather than as a definitive global test.
#'   - Rows with missing values are silently dropped before testing; this
#'     may bias the analysis if missingness is informative.
#'
#' @references
#'   - Anderson, M. J. (2001). A new method for non‐parametric multivariate
#'     analysis of variance. Austral Ecology, 26, 32-46.
#'     \doi{10.1111/j.1442-9993.2001.01070.pp.x}.
#'
#'   - Anderson, M. J. (2006). Distance‐Based Tests for Homogeneity of
#'     Multivariate Dispersions. Biometrics, 62, 245-253.
#'     \doi{10.1111/j.1541-0420.2005.00440.x}.
#'
#'   - Hotelling, H. (1931). The Generalization of Student's Ratio. The
#'     Annals of Mathematical Statistics, 2, 360-378.
#'     \doi{10.1214/aoms/1177732979}.
#'
#'   - Benjamini, Y. and Hochberg, Y. (1995). Controlling the False
#'     Discovery Rate: A Practical and Powerful Approach to Multiple
#'     Testing. Journal of the Royal Statistical Society Series B:
#'     Statistical Methodology, 57, 289-300.
#'     \doi{10.1111/j.2517-6161.1995.tb02031.x}.
#'
#'   - Warton, D. I., Wright, S. T. and Wang, Y. (2012). Distance‐based
#'     multivariate analyses confound location and dispersion effects.
#'     Methods in Ecology and Evolution, 3, 89-101.
#'     \doi{10.1111/j.2041-210X.2011.00127.x}.
#'
#' @examples
#' set.seed(1)
#' # 12 spectra x 50 wavelengths, 2 groups differing in 5 channels
#' wl <- seq(400, 798, length.out = 50)
#' Y <- matrix(rnorm(12 * 50, sd = 0.3), nrow = 12, ncol = 50)
#' Y[1:6, 21:25] <- Y[1:6, 21:25] + 0.8
#' g <- factor(rep(c("ctrl", "trt"), each = 6))
#'
#' # PERMANOVA on the full spectra
#' spectral_change_test(Y, g, wavelengths = wl,
#'                      method = "permanova", n_perm = 199)
#'
#' # Hotelling T^2
#' spectral_change_test(Y, g, method = "t2")
#'
#' # Pointwise t-tests with BH adjustment
#' r <- spectral_change_test(Y, g, wavelengths = wl, method = "pointwise_t")
#' head(r$pointwise)
#'
#' @export
spectral_change_test <- function(spectra, group, wavelengths = NULL,
                                  method = c("t2", "permanova",
                                             "pointwise_t"),
                                  conf_level = 0.95, n_perm = 999L,
                                  seed = NULL, adjust = "fdr") {
  call <- match.call()
  method <- match.arg(method)
  if (!is.matrix(spectra) && !is.data.frame(spectra)) {
    stop("'spectra' must be a matrix or data frame.", call. = FALSE)
  }
  spectra <- as.matrix(spectra)
  if (!is.numeric(spectra)) {
    stop("'spectra' must be numeric.", call. = FALSE)
  }
  if (is.null(wavelengths)) {
    wavelengths <- seq_len(ncol(spectra))
  } else {
    .assert_numeric_vector(wavelengths, "wavelengths", min_length = 1L)
    if (length(wavelengths) != ncol(spectra)) {
      stop("Length of 'wavelengths' must match ncol(spectra).", call. = FALSE)
    }
  }
  if (length(group) != nrow(spectra)) {
    stop("Length of 'group' must match nrow(spectra).", call. = FALSE)
  }
  group <- as.factor(group)
  if (nlevels(group) < 2L) {
    stop("'group' must have at least 2 levels.", call. = FALSE)
  }
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(n_perm, "n_perm")
  if (!is.character(adjust) || length(adjust) != 1L) {
    stop("'adjust' must be a single character string.", call. = FALSE)
  }

  ok <- complete.cases(spectra, group)
  if (!all(ok)) {
    spectra <- spectra[ok, , drop = FALSE]
    group <- droplevels(group[ok])
  }
  n <- nrow(spectra)
  p <- ncol(spectra)
  g <- nlevels(group)

  if (method == "t2") {
    if (g != 2L) {
      stop("Hotelling T^2 requires exactly 2 groups. Use 'permanova' for >2.",
           call. = FALSE)
    }
    lev <- levels(group)
    g1 <- spectra[group == lev[1L], , drop = FALSE]
    g2 <- spectra[group == lev[2L], , drop = FALSE]
    n1 <- nrow(g1); n2 <- nrow(g2)
    if (n1 < 2L || n2 < 2L) {
      stop("Each group must have at least 2 observations.", call. = FALSE)
    }
    ybar1 <- colMeans(g1)
    ybar2 <- colMeans(g2)
    diff <- ybar1 - ybar2
    S1 <- stats::cov(g1)
    S2 <- stats::cov(g2)
    Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    Sp_inv <- tryCatch(solve(Sp), error = function(e) NULL)
    if (is.null(Sp_inv)) {
      warning(sprintf(
        "Hotelling T^2 requires invertible pooled covariance; here p=%d, n1+n2-1=%d. ",
        p, n1 + n2 - 1L),
        "Use method='permanova' instead; returning p_value=NA.",
        call. = FALSE, immediate. = TRUE)
      Sp_inv <- MASS::ginv(Sp)
      T2 <- (n1 * n2 / (n1 + n2)) * sum(diff * (Sp_inv %*% diff))
      stat <- T2
      p_value <- NA_real_
      method_label <- sprintf("Hotelling T^2 (singular Sp, p=%d, n1+n2-1=%d; p_value=NA)",
                              p, n1 + n2 - 1L)
    } else {
      T2 <- (n1 * n2 / (n1 + n2)) * sum(diff * (Sp_inv %*% diff))
      df1 <- p
      df2 <- n1 + n2 - p - 1L
      if (df2 > 0L) {
        F_stat <- ((n1 + n2 - p - 1L) / (p * (n1 + n2 - 2L))) * T2
        p_value <- stats::pf(F_stat, df1, df2, lower.tail = FALSE)
        stat <- F_stat
        method_label <- sprintf("Hotelling T^2 (F-transform, df1=%d, df2=%d)",
                                df1, df2)
      } else {
        warning(sprintf(
          "Hotelling T^2: df2 = n1+n2-p-1 = %d <= 0 (p=%d, n1+n2-1=%d). ",
          df2, p, n1 + n2 - 1L),
          "Use method='permanova' instead; returning p_value=NA.",
          call. = FALSE, immediate. = TRUE)
        p_value <- NA_real_
        stat <- T2
        method_label <- "Hotelling T^2 (df2<=0; p_value=NA)"
      }
    }
    out <- list(test_statistic = unname(stat),
                p_value = p_value,
                method = method_label,
                n_perm = NA_integer_,
                pointwise = NULL,
                call = call)
  } else if (method == "permanova") {
    old <- .set_seed(seed)
    on.exit(.restore_seed(old), add = TRUE)
    D <- as.matrix(dist(spectra))^2
    F_obs <- .permanova_F(D, group)
    if (!is.finite(F_obs)) {
      stop("PERMANOVA F-statistic is degenerate (likely zero within-group SS).",
           call. = FALSE)
    }
    perm_Fs <- replicate(n_perm, {
      pg <- sample(group)
      .permanova_F(D, pg)
    })
    p_value <- (1 + sum(perm_Fs >= F_obs, na.rm = TRUE)) / (n_perm + 1)
    out <- list(test_statistic = unname(F_obs),
                p_value = p_value,
                method = sprintf("PERMANOVA (Euclidean, %d permutations)",
                                 n_perm),
                n_perm = as.integer(n_perm),
                pointwise = NULL,
                call = call)
  } else {
    if (g != 2L) {
      stop("pointwise_t currently supports only 2 groups.", call. = FALSE)
    }
    lev <- levels(group)
    g1 <- spectra[group == lev[1L], , drop = FALSE]
    g2 <- spectra[group == lev[2L], , drop = FALSE]
    tvals <- vapply(seq_len(p), function(j) {
      unname(stats::t.test(g1[, j], g2[, j])$statistic)
    }, numeric(1))
    pvals <- vapply(seq_len(p), function(j) {
      stats::t.test(g1[, j], g2[, j])$p.value
    }, numeric(1))
    p_adj <- stats::p.adjust(pvals, method = adjust)
    pw_df <- data.frame(
      wavelength = wavelengths,
      t = tvals,
      p = pvals,
      p_adj = p_adj,
      significant = p_adj < (1 - conf_level),
      stringsAsFactors = FALSE
    )
    out <- list(test_statistic = max(abs(tvals)),
                p_value = min(p_adj),
                method = sprintf("Pointwise t-tests (%s-adjusted)", adjust),
                n_perm = NA_integer_,
                pointwise = pw_df,
                n_significant = sum(pw_df$significant),
                call = call)
  }

  out$conf_level <- conf_level
  out$tidy <- if (!is.null(out$pointwise)) out$pointwise else
    data.frame(test_statistic = out$test_statistic,
               p_value = out$p_value,
               method = out$method,
               stringsAsFactors = FALSE)

  .crossstatr_result(out, "crossstatr_spectral_change")
}

# PERMANOVA pseudo-F statistic from squared distance matrix and grouping.
# @noRd
.permanova_F <- function(D, group) {
  n <- nrow(D)
  group <- as.factor(group)
  g <- nlevels(group)
  if (g < 2L || n < g + 1L) return(NA_real_)
  SSt <- sum(D) / n
  SSw <- 0
  for (lv in levels(group)) {
    idx <- which(group == lv)
    ng <- length(idx)
    if (ng > 0L) {
      SSw <- SSw + sum(D[idx, idx, drop = FALSE]) / ng
    }
  }
  SSa <- SSt - SSw
  df_between <- g - 1
  df_within <- n - g
  if (df_within <= 0L || SSw <= 0) return(NA_real_)
  (SSa / df_between) / (SSw / df_within)
}

#' @export
print.crossstatr_spectral_change <- function(x, ...) {
  cat(sprintf("<crossStatR::spectral_change_test> %s\n", x$method))
  cat(sprintf("test_statistic = %.4g,  p_value = %g\n",
              x$test_statistic, x$p_value))
  if (is.na(x$n_perm)) {
    cat("Permutations: (not used)\n")
  } else {
    cat(sprintf("Permutations: %d\n", x$n_perm))
  }
  if (!is.null(x$pointwise)) {
    cat(sprintf("\nPointwise tests: %d/%d significant at level %g\n",
                sum(x$pointwise$significant), nrow(x$pointwise),
                1 - x$conf_level))
  }
  invisible(x)
}
