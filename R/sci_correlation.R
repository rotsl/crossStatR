#' Correlation analysis with confidence intervals
#'
#' Computes Pearson, Spearman, or Kendall correlations between two numeric
#' vectors with confidence intervals derived from Fisher's variance-stabilising
#' z-transform. When a grouping variable is supplied, the correlation is
#' computed within each group and the resulting p-values are adjusted for
#' multiplicity across groups. This is the natural use case when the analyst has
#' repeated within-subject or within-condition measures and want to know
#' whether the association generalises across strata, rather than trusting a
#' single pooled estimate that may be driven by one subgroup.
#'
#' The package deliberately do not pre-test for normality and switch methods on that
#' basis. The choice between Pearson (linear, sensitive to outliers, exact CI
#' under bivariate normality) and Spearman/Kendall (rank-based, robust to
#' monotone non-linearities and outliers) belongs to the analyst and should be
#' motivated by the measurement scale and the underlying science.
#'
#' @param x numeric vector.
#' @param y numeric vector of the same length as `x`.
#' @param method character; `"pearson"` (default), `"spearman"`, or
#'   `"kendall"`.
#' @param alternative `"two.sided"` (default), `"less"`, or `"greater"`.
#' @param conf_level confidence level for the interval on the correlation
#'   (default 0.95).
#' @param group optional grouping variable (factor or coercible) of the same
#'   length as `x`. If supplied, per-group correlations are computed and the
#'   p-values in `groupwise` are adjusted using `p_adjust`.
#' @param p_adjust character; p-value adjustment method for the groupwise
#'   p-values (default `"holm"`). Any method accepted by [stats::p.adjust()]
#'   is allowed.
#'
#' @return A `crossstatr_correlation` object with elements:
#'   \item{method}{correlation method used.}
#'   \item{correlation}{numeric scalar point estimate (over all data).}
#'   \item{p_value}{p-value of the test of `H0: \rho = 0`.}
#'   \item{ci}{length-2 numeric confidence interval.}
#'   \item{n}{number of complete (finite) observations used.}
#'   \item{alternative}{alternative hypothesis.}
#'   \item{conf_level}{confidence level used.}
#'   \item{groupwise}{data frame of per-group estimates and adjusted p-values,
#'     or `NULL` if `group` was not supplied.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   For Pearson's `r` The package uses Fisher's variance-stabilising transform
#'   \deqn{z = \mathrm{atanh}(r) = \tfrac{1}{2}\ln\!\frac{1+r}{1-r},}
#'   which is approximately Gaussian with variance \eqn{1/(n-3)} under the
#'   null of no correlation (Fisher, 1915; 1921). The confidence interval on
#'   `z` is
#'   \deqn{z \pm t_{1-\alpha/2,\,n-3}\,\frac{1}{\sqrt{n-3}},}
#'   which is then back-transformed to the correlation scale via
#'   \eqn{r = \mathrm{tanh}(z)}. The package uses a Student-t critical value rather than
#'   the Gaussian one; this gives a slightly wider interval for small `n` and
#'   is what `stats::cor.test` uses internally for Pearson intervals.
#'
#'   For Spearman's \eqn{\rho}, the same Fisher-z approximation is applied to
#'   the point estimate. This works well for moderate `n` (Bonett & Wright,
#'   2000), although for very small samples or extreme correlations the
#'   coverage can be slightly off.
#'
#'   For Kendall's \eqn{\tau}, `stats::cor.test` is called with
#'   `exact = FALSE` to avoid the combinatorial cost of the exact variance
#'   computation when `n` is large. The asymptotic variance used is the one
#'   for Kendall rank association. Confidence intervals for Kendall's tau are also
#'   constructed via the Fisher-z approximation on `atanh(tau)`.
#'
#'   When a `group` is supplied, p-values are adjusted across groups using
#'   `stats::p.adjust` with the requested method.
#'
#' @section Limitations:
#'
#'   - The Fisher-z interval for Spearman and Kendall is approximate and can
#'     be inaccurate for very small samples or when the correlation is close
#'     to \eqn{\pm 1}. For `n < 10` consider exact permutation or bootstrap
#'     intervals instead.
#'   - Kendall's tau with ties uses the normal approximation built into
#'     `cor.test`; the exact variance is not computed.
#'   - Groupwise correlations discard groups with fewer than 3 finite
#'     observations. The reported adjusted p-values correspond only to the
#'     groups that were actually computed; do not interpret them as a global
#'     omnibus test.
#'   - Pooling across groups via a single correlation ignores possible
#'     between-group heterogeneity. If the analyst supplies a `group`, read the
#'     `groupwise` table; the top-level `correlation` is the pooled estimate
#'     and may be misleading if the groups differ in location or scale.
#'
#' @references
#'   - Fisher, R. A. (1915). Frequency Distribution of the Values of the
#'     Correlation Coefficient in Samples from an Indefinitely Large
#'     Population. Biometrika, 10, 507-521.
#'     \doi{10.2307/2331838}.
#'
#'   - Fisher, R. A. (1921). On the probable error of a coefficient of
#'     correlation deduced from a small sample. Metron, 1, 3-32.
#'     \url{https://digital.library.adelaide.edu.au/items/002ad8fb-c23c-407b-8a89-f036d8da6030}.
#'
#'   - Bonett, D. G. and Wright, T. A. (2000). Sample Size Requirements for
#'     Estimating Pearson, Kendall and Spearman Correlations.
#'     Psychometrika, 65, 23-28.
#'     \doi{10.1007/BF02294183}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(40)
#' y <- 0.6 * x + rnorm(40, sd = 0.8)
#' result <- sci_correlation(x, y, method = "pearson")
#' print(result)
#'
#' # Groupwise
#' g <- rep(c("A", "B"), each = 20)
#' x2 <- c(rnorm(20), rnorm(20) + 1)
#' y2 <- c(rnorm(20), 0.8 * rnorm(20) + 1)
#' sci_correlation(x2, y2, group = g, method = "spearman")
#'
#' @export
sci_correlation <- function(x, y,
                            method = c("pearson", "spearman", "kendall"),
                            alternative = c("two.sided", "less", "greater"),
                            conf_level = 0.95, group = NULL,
                            p_adjust = "holm") {
  call <- match.call()
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  .assert_probability(conf_level, "conf_level")
  .assert_numeric_vector(x, "x", min_length = 3L, finite = FALSE)
  .assert_numeric_vector(y, "y", min_length = 3L, finite = FALSE)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  if (!is.null(group) && length(group) != length(x)) {
    stop("'group' must have the same length as 'x'.", call. = FALSE)
  }
  if (!is.character(p_adjust) || length(p_adjust) != 1L) {
    stop("'p_adjust' must be a single character string.", call. = FALSE)
  }

  # Overall correlation on all data
  overall <- .compute_correlation(x, y, method, alternative, conf_level)

  groupwise_df <- NULL
  if (!is.null(group)) {
    groups <- if (is.factor(group)) levels(droplevels(factor(group))) else
      sort(unique(group))
    groups <- as.character(groups)
    rows <- lapply(groups, function(g) {
      idx <- (group == g) & is.finite(x) & is.finite(y)
      n_g <- sum(idx)
      if (n_g < 3L) {
        return(data.frame(group = g, method = method,
                          correlation = NA_real_, p_value = NA_real_,
                          ci_low = NA_real_, ci_high = NA_real_,
                          n = n_g, p_adj = NA_real_,
                          stringsAsFactors = FALSE))
      }
      rr <- .compute_correlation(x[idx], y[idx], method, alternative,
                                 conf_level)
      data.frame(group = g, method = method,
                 correlation = rr$correlation, p_value = rr$p_value,
                 ci_low = rr$ci[1L], ci_high = rr$ci[2L],
                 n = rr$n, p_adj = NA_real_,
                 stringsAsFactors = FALSE)
    })
    groupwise_df <- do.call(rbind, rows)
    # Adjust across the groups for which a calculation was possible for a p-value
    valid <- !is.na(groupwise_df$p_value)
    if (any(valid)) {
      groupwise_df$p_adj[valid] <- .adjust_p(groupwise_df$p_value[valid],
                                             method = p_adjust)
    }
  }

  out <- .crossstatr_result(
    list(method = method,
         correlation = overall$correlation,
         p_value = overall$p_value,
         ci = overall$ci,
         n = overall$n,
         alternative = alternative,
         conf_level = conf_level,
         groupwise = groupwise_df,
         p_adjust_method = if (is.null(group)) NULL else p_adjust,
         call = call),
    "crossstatr_correlation"
  )
  out
}

# Compute correlation, p-value, and Fisher-z CI for a single (x, y) pair.
# Returns a list with correlation, p_value, ci (length 2), n.
# @noRd
.compute_correlation <- function(x, y, method, alternative, conf_level) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3L) {
    return(list(correlation = NA_real_, p_value = NA_real_,
                 ci = c(NA_real_, NA_real_), n = n))
  }
  alpha <- 1 - conf_level
  # Kendall: disable exact computation (costly for large n)
  if (method == "kendall") {
    ct <- stats::cor.test(x, y, method = "kendall",
                         alternative = alternative, exact = FALSE)
  } else {
    ct <- stats::cor.test(x, y, method = method,
                         alternative = alternative, exact = FALSE)
  }
  r <- as.numeric(ct$estimate)
  p <- ct$p.value
  # Fisher-z CI (works for Pearson, approximate for Spearman/Kendall)
  if (is.finite(r) && abs(r) < 1) {
    z <- atanh(r)
    se_z <- 1 / sqrt(n - 3)
    if (alternative == "two.sided") {
      tcrit <- stats::qt(1 - alpha / 2, df = n - 3)
      ci_z_lo <- z - tcrit * se_z
      ci_z_hi <- z + tcrit * se_z
    } else if (alternative == "greater") {
      tcrit <- stats::qt(1 - alpha, df = n - 3)
      ci_z_lo <- z - tcrit * se_z
      ci_z_hi <- Inf
    } else {
      tcrit <- stats::qt(1 - alpha, df = n - 3)
      ci_z_lo <- -Inf
      ci_z_hi <- z + tcrit * se_z
    }
    ci <- c(tanh(ci_z_lo), tanh(ci_z_hi))
  } else {
    # |r| == 1 or NA: degenerate interval
    ci <- c(r, r)
  }
  list(correlation = r, p_value = p, ci = ci, n = n)
}

#' @export
print.crossstatr_correlation <- function(x, ...) {
  cat(sprintf("<crossStatR::correlation> method = %s\n", x$method))
  cat(sprintf("r = %.4f,  p = %g,  n = %d\n",
              x$correlation, x$p_value, x$n))
  cat(sprintf("%g%% CI: [%.4f, %.4f]  (alternative: %s)\n",
              100 * x$conf_level, x$ci[1L], x$ci[2L], x$alternative))
  if (!is.null(x$groupwise)) {
    cat("\nGroupwise:\n")
    print(x$groupwise, row.names = FALSE)
  }
  invisible(x)
}
