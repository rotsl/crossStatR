#' Statistical identification of regions where experimental curves differ
#'
#' Constructs simultaneous confidence bands for the difference between
#' two (or two groups of) experimental curves, and identifies contiguous
#' regions where the bands exclude zero. The method handles correlated
#' measurements, unequal sampling grids, and replicate variability. It
#' is genuinely useful for spectroscopic and electrochemical datasets
#' where naive pointwise testing inflates the family-wise error rate.
#'
#' @param x numeric vector of x-values.
#' @param y1 numeric matrix or data frame of curves in group 1. Each row
#'   is a replicate curve observed at `x`. If curves are observed on
#'   different grids, supply `x1`, `x2` instead.
#' @param y2 numeric matrix or data frame of curves in group 2.
#' @param x1,x2 optional x-vectors for each group; if supplied, `x` is
#'   ignored and both groups are interpolated onto `x_grid` (or a union
#'   of the two grids).
#' @param x_grid optional target grid; if NULL, the union of `x1` and
#'   `x2` (or `x`) is used, restricted to the common range.
#' @param conf_level nominal coverage of the simultaneous band
#'   (default 0.95).
#' @param method `"max-t"` (default) for the max-t simultaneous band
#'   (max-statistic resampling; see Westfall and Young, 1993, when
#'   `n_perm > 0`). The accepted `"boot"` value currently follows the same
#'   computation; a separate bootstrap band is not implemented.
#' @param n_perm number of permutations for the Westfall-Young
#'   adjustment. 0 disables permutation adjustment (faster, more
#'   conservative).
#' @param seed optional integer seed.
#' @param min_region_width minimum width (in x units) for a region to
#'   be reported. Shorter regions are merged into neighbours.
#'
#' @return A `crossstatr_sdr` object with:
#'   \item{x_grid}{common grid.}
#'   \item{diff}{pointwise mean difference.}
#'   \item{se}{pointwise standard error.}
#'   \item{band_low, band_high}{lower/upper simultaneous band.}
#'   \item{regions}{data frame of contiguous significant regions
#'     (start, end, max_abs_diff, sign).}
#'   \item{n_perm}{number of permutations used.}
#'   \item{method}{method label.}
#'
#' @section Methodology:
#'
#'   Let \eqn{\hat d(x) = \bar y_1(x) - \bar y_2(x)} be the pointwise
#'   mean difference at grid point \eqn{x}, and let
#'   \eqn{s(x) = \sqrt{s_1^2(x)/n_1 + s_2^2(x)/n_2}} be its standard
#'   error. The pointwise t-statistic is
#'   \deqn{T(x) = \hat d(x) / s(x).}
#'
#'   A naive pointwise 95% interval \eqn{\hat d(x) \pm t_{0.975}\,s(x)}
#'   has family-wise coverage far below 95% when the family is the
#'   whole grid. The max-t band instead uses the critical value
#'   \deqn{c_\alpha = Q_{1-\alpha}\left[\max_x |T(x)|\right],}
#'   where the quantile is taken over the distribution of the maximum
#'   absolute t-statistic under the null. The simultaneous band is
#'   then \eqn{\hat d(x) \pm c_\alpha\, s(x)}.
#'
#'   The distribution of the max-t statistic is estimated by a
#'   Westfall-Young permutation scheme: under the null hypothesis of
#'   no group difference, group labels are exchanged between replicate
#'   curves; the max-t statistic is recomputed for each of `n_perm`
#'   permutations; the empirical \eqn{1 - \alpha} quantile of these
#'   permuted maxima is \eqn{c_\alpha}.
#'
#'   When `n_perm = 0`, the Bonferroni-adjusted critical value is used:
#'   \eqn{c_\alpha = t_{1 - \alpha/(2m), \nu}} with \eqn{m} the number
#'   of grid points. This is conservative but does not require
#'   permutation.
#'
#'   Contiguous significant regions are extracted as runs of grid
#'   points where the band excludes zero. Region boundaries are reported
#'   with explicit x-values; the width uncertainty of each boundary is
#'   one grid step.
#'
#' @section Limitations:
#'
#'   - The Bonferroni fallback uses a common degrees-of-freedom approximation
#'     rather than propagating pointwise Welch df. Unequal variances and small
#'     replicate counts can affect simultaneous coverage.
#'
#'   - Assumes replicate curves are independent. For paired designs,
#'     permute within pairs (not implemented; use paired-test
#'     infrastructure instead).
#'   - The permutation exchangeability assumption is violated when
#'     groups have different noise structures. In that case, use
#'     `n_perm = 0` (Bonferroni) instead.
#'   - Region boundaries are limited by the grid resolution; sub-grid
#'     localization requires interpolation and is not done here.
#'   - For very large grids (>1000 points) the Bonferroni adjustment
#'     becomes extremely conservative; permutation is recommended.
#'
#' @references
#'   - Westfall, P. H. and Young, S. S. (1993). Resampling-Based Multiple
#'     Testing: Examples and Methods for p-Value Adjustment. Wiley.
#'     \url{https://www.wiley.com/en-gb/Resampling-Based+Multiple+Testing%3A+Examples+and+Methods+for+p-Value+Adjustment-p-9780471557616}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(400, 800, length.out = 80)
#' # Two groups of curves with a true difference in 500-560 region
#' make_curve <- function(mu, region) {
#'   base <- sin((x - 400) * pi / 200) * 0.2 + mu
#'   if (region) base <- base + dnorm(x, mean = 530, sd = 18) * 0.8
#'   base + rnorm(length(x), sd = 0.15)
#' }
#' y1 <- matrix(replicate(15, make_curve(1, FALSE)), nrow = 15, byrow = TRUE)
#' y2 <- matrix(replicate(15, make_curve(1, TRUE)),  nrow = 15, byrow = TRUE)
#' res <- signal_difference_regions(x, y1, y2, n_perm = 200)
#' print(res)
#'
#' @export
signal_difference_regions <- function(x, y1, y2, x1 = NULL, x2 = NULL,
                                       x_grid = NULL, conf_level = 0.95,
                                       method = c("max-t", "boot"),
                                       n_perm = 0L, seed = NULL,
                                       min_region_width = 0) {
  call <- match.call()
  method <- match.arg(method)
  .assert_probability(conf_level, "conf_level")

  # Coerce to matrices
  if (!is.matrix(y1)) y1 <- as.matrix(y1)
  if (!is.matrix(y2)) y2 <- as.matrix(y2)
  if (ncol(y1) != length(x) && !is.null(x1)) {
    # use x1 as the per-column x for group 1
  } else if (ncol(y1) != length(x)) {
    stop("ncol(y1) must match length(x).", call. = FALSE)
  }
  if (ncol(y2) != length(x) && !is.null(x2)) {
    # use x2 as the per-column x for group 2
  } else if (ncol(y2) != length(x)) {
    stop("ncol(y2) must match length(x).", call. = FALSE)
  }

  n1 <- nrow(y1); n2 <- nrow(y2)
  if (n1 < 2L || n2 < 2L) {
    stop("Each group needs at least 2 replicate curves.", call. = FALSE)
  }

  # Interpolate onto common grid if x1/x2 supplied
  if (!is.null(x1) && !is.null(x2)) {
    if (is.null(x_grid)) {
      xr1 <- range(x1); xr2 <- range(x2)
      lo <- max(xr1[1L], xr2[1L]); hi <- min(xr1[2L], xr2[2L])
      x_grid <- sort(unique(c(x1, x2)))
      x_grid <- x_grid[x_grid >= lo & x_grid <= hi]
    }
    y1g <- t(apply(y1, 1L, function(row) .interp_to_grid(x1, row, x_grid)))
    y2g <- t(apply(y2, 1L, function(row) .interp_to_grid(x2, row, x_grid)))
  } else {
    x_grid <- x
    y1g <- y1; y2g <- y2
  }

  m <- length(x_grid)
  if (m < 2L) stop("Common grid too small.", call. = FALSE)

  # Pointwise difference and SE
  d_hat <- colMeans(y1g) - colMeans(y2g)
  v1 <- apply(y1g, 2L, var) / n1
  v2 <- apply(y2g, 2L, var) / n2
  s_hat <- sqrt(v1 + v2)
  # Avoid division by zero
  s_hat[s_hat < .Machine$double.eps] <- NA_real_
  t_stat <- d_hat / s_hat

  # Effective df via Welch-Satterthwaite
  df_eff <- (v1 + v2)^2 / (v1^2 / pmax(n1 - 1, 1L) + v2^2 / pmax(n2 - 1, 1L))
  df_eff[is.na(df_eff) | !is.finite(df_eff)] <- n1 + n2 - 2

  # Critical value for the simultaneous band
  if (n_perm > 0L) {
    .set_seed(seed)
    # Westfall-Young max-t: permute group labels, recompute max|T|
    pooled <- rbind(y1g, y2g)
    labels <- c(rep(1L, n1), rep(2L, n2))
    max_T_perm <- replicate(n_perm, {
      sp <- sample(labels)
      g1 <- pooled[sp == 1L, , drop = FALSE]
      g2 <- pooled[sp == 2L, , drop = FALSE]
      d_p <- colMeans(g1) - colMeans(g2)
      v1_p <- apply(g1, 2L, var) / sum(sp == 1L)
      v2_p <- apply(g2, 2L, var) / sum(sp == 2L)
      s_p <- sqrt(v1_p + v2_p)
      s_p[s_p < .Machine$double.eps] <- NA
      T_p <- d_p / s_p
      max(abs(T_p), na.rm = TRUE)
    })
    c_alpha <- stats::quantile(max_T_perm, probs = conf_level, na.rm = TRUE)
  } else {
    # Bonferroni: t_{1 - alpha/(2m), df}
    alpha <- 1 - conf_level
    c_alpha <- stats::qt(1 - alpha / (2 * m), df = max(1L, min(df_eff, na.rm = TRUE)))
  }

  band_low  <- d_hat - c_alpha * s_hat
  band_high <- d_hat + c_alpha * s_hat

  # Extract contiguous significant regions
  significant <- (!is.na(band_low)) & (band_low > 0 | band_high < 0)
  regions <- .extract_regions(x_grid, significant, d_hat, min_region_width)

  .crossstatr_result(
    list(x_grid = x_grid,
         diff = d_hat,
         se = s_hat,
         t_stat = t_stat,
         band_low = band_low,
         band_high = band_high,
         critical_value = c_alpha,
         regions = regions,
         n_perm = n_perm,
         n1 = n1, n2 = n2,
         method = if (n_perm > 0L)
           sprintf("max-t with Westfall-Young permutation (%d perms)", n_perm)
         else sprintf("max-t (Bonferroni, m=%d)", m),
         conf_level = conf_level,
         call = call),
    "crossstatr_sdr")
}

.extract_regions <- function(x, sig, d, min_width) {
  if (length(sig) == 0L || !any(sig)) {
    return(data.frame(start = numeric(0), end = numeric(0),
                      n_points = integer(0),
                      max_abs_diff = numeric(0),
                      sign = character(0)))
  }
  rle_obj <- rle(sig)
  ends <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1L
  # Keep only runs where value is TRUE
  keep <- rle_obj$values
  starts <- starts[keep]; ends <- ends[keep]
  res <- data.frame(
    start = x[starts], end = x[ends],
    n_points = ends - starts + 1L,
    max_abs_diff = vapply(seq_along(starts), function(k)
      max(abs(d[starts[k]:ends[k]])), numeric(1L)),
    sign = vapply(seq_along(starts), function(k)
      ifelse(mean(d[starts[k]:ends[k]]) > 0, "+", "-"), character(1L))
  )
  # Merge regions narrower than min_width into neighbours (drop them)
  if (min_width > 0 && nrow(res) > 0L) {
    res <- res[res$end - res$start >= min_width, , drop = FALSE]
  }
  res
}

#' @export
print.crossstatr_sdr <- function(x, ...) {
  cat(sprintf("<crossStatR::signal_difference_regions> %s\n", x$method))
  cat(sprintf("Nominal coverage: %.0f%%\n", 100 * x$conf_level))
  cat(sprintf("Critical value c_alpha = %.4f\n", x$critical_value))
  cat(sprintf("Grid size: %d points; groups n1=%d, n2=%d\n",
              length(x$x_grid), x$n1, x$n2))
  if (nrow(x$regions) == 0L) {
    cat("No significant regions detected.\n")
  } else {
    cat("\nSignificant regions:\n")
    print(x$regions, row.names = FALSE)
  }
  invisible(x)
}
