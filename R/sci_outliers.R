#' Outlier diagnostics for experimental data
#'
#' Flags observations that are statistically extreme relative to a
#' distributional model. **Flagging does not remove observations.** The
#' function returns a per-observation diagnostic table; removal is left
#' to the analyst, who should record any decision with justification.
#'
#' @param x numeric vector.
#' @param group optional grouping vector.
#' @param method character vector; which diagnostics to compute. Subset of
#'   `c("iqr", "mad", "zscore", "robust_z")`.
#' @param iqr_k IQR multiplier (default 1.5; use 3.0 for "extreme" only).
#' @param z_threshold absolute z-score threshold (default 3.0).
#' @param robust_z_threshold absolute MAD-based z-score threshold
#'   (default 3.5, Iglewicz & Hoaglin, 1993).
#'
#' @return A `crossstatr_outliers` object with elements:
#'   \item{flags}{data frame with one row per observation.}
#'   \item{counts}{per-method flag counts per group.}
#'   \item{method}{method label.}
#'
#' @section Methodology:
#'
#'   IQR rule (Tukey, 1977):
#'   \deqn{[Q_1 - k(Q_3 - Q_1),\; Q_3 + k(Q_3 - Q_1)]}
#'
#'   MAD rule (Hampel, 1974):
#'   \deqn{[m - k\,\mathrm{MAD},\; m + k\,\mathrm{MAD}],}
#'   with `MAD` rescaled by 1.4826 to match the Gaussian scale. The package uses
#'   `k = 3.5` by default.
#'
#'   Standardized z-score:
#'   \deqn{z_i = (x_i - \bar{x}) / s,}
#'   flagged when \eqn{|z_i| > z_{\mathrm{threshold}}}.
#'
#'   Robust z-score (Iglewicz & Hoaglin, 1993):
#'   \deqn{z^R_i = 0.6745\,(x_i - m)/\mathrm{MAD},}
#'   flagged when \eqn{|z^R_i| > 3.5}.
#'
#' @section Limitations:
#'
#'   - Robust z scores divide deviations by the normal-consistent MAD.
#'     If the MAD is zero, median-valued observations receive zero and
#'     nonzero deviations receive signed infinity, preserving their flags.
#'
#'   - For very small samples (n < 5) the IQR/MAD rules have low power
#'     and the z-score rule is unstable. Results are returned with a
#'     `n < 5` warning.
#'   - Groupwise diagnostics assume the group itself is homogeneous; if
#'     within-group structure exists (e.g. time trends), use model-based
#'     influence diagnostics from [sci_regression()] instead.
#'   - Multiple-test inflation is *not* corrected here because the
#'     diagnostics are exploratory.
#'
#' @references
#'   - Tukey, J. W. (1977). Exploratory Data Analysis. Addison-Wesley. ISBN
#'     9780201076165.
#'     \url{https://www.pearson.com/en-us/subject-catalog/p/exploratory-data-analysis-classic-version-/P200000006400}.
#'
#'   - Hampel, F. R. (1974). The Influence Curve and its Role in Robust
#'     Estimation. Journal of the American Statistical Association, 69,
#'     383-393.
#'     \doi{10.1080/01621459.1974.10482962}.
#'
#'   - Iglewicz, B. and Hoaglin, D. C. (1993). How to Detect and Handle
#'     Outliers. ASQC Quality Press.
#'     \url{https://asq.org/quality-press/display-item?item=E0801}.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(40), 8, 12)
#' sci_outliers(x)
#'
#' @export
sci_outliers <- function(x, group = NULL,
                        method = c("iqr", "mad", "zscore", "robust_z"),
                        iqr_k = 1.5, z_threshold = 3.0,
                        robust_z_threshold = 3.5) {
  call <- match.call()
  .assert_numeric_vector(x, "x", min_length = 1L, finite = FALSE)
  if (is.null(group)) group <- factor(rep(1L, length(x)))
  if (length(x) != length(group)) {
    stop("Length of 'group' must match 'x'.", call. = FALSE)
  }
  method <- match.arg(method, several.ok = TRUE)

  out <- do.call(rbind, lapply(as.character(sort(unique(group))), function(g) {
    idx <- which(group == g)
    xi <- x[idx]
    keep <- is.finite(xi)
    xi_f <- xi[keep]
    n <- length(xi_f)
    if (n < 2L) {
      return(data.frame(
        index = idx[keep], group = g, value = xi_f,
        iqr_flag = NA, mad_flag = NA, zscore = NA, z_flag = NA,
        robust_z = NA, robust_z_flag = NA,
        stringsAsFactors = FALSE))
    }
    q <- stats::quantile(xi_f, c(0.25, 0.5, 0.75))
    m <- q[2L]; q1 <- q[1L]; q3 <- q[3L]
    iqr_f <- q3 - q1
    iqr_lo <- q1 - iqr_k * iqr_f; iqr_hi <- q3 + iqr_k * iqr_f
    s <- stats::sd(xi_f)
    mad_v <- stats::mad(xi_f, constant = 1.4826)

    flag_iqr <- if ("iqr" %in% method)
      (xi_f < iqr_lo) | (xi_f > iqr_hi) else NA
    flag_mad <- if ("mad" %in% method)
      abs(xi_f - m) > robust_z_threshold * mad_v else NA

    zsc <- if ("zscore" %in% method) {
      if (is.finite(s) && s > 0) (xi_f - mean(xi_f)) / s else rep(0, length(xi_f))
    } else NA
    rzsc <- if ("robust_z" %in% method) {
      if (is.finite(mad_v) && mad_v > 0) (xi_f - m) / mad_v
      else ifelse(xi_f == m, 0, sign(xi_f - m) * Inf)
    } else NA

    z_flag <- if ("zscore" %in% method) abs(zsc) > z_threshold else NA
    rz_flag <- if ("robust_z" %in% method)
      abs(rzsc) > robust_z_threshold else NA

    data.frame(
      index = idx[keep], group = g, value = xi_f,
      iqr_flag = flag_iqr, mad_flag = flag_mad,
      zscore = zsc, z_flag = z_flag,
      robust_z = rzsc, robust_z_flag = rz_flag,
      stringsAsFactors = FALSE)
  }))

  counts <- data.frame(
    group = unique(out$group),
    n_total = as.integer(table(out$group)[unique(out$group)]),
    iqr_flags = if ("iqr" %in% method)
      as.integer(tapply(out$iqr_flag, out$group, sum, na.rm = TRUE)) else NA,
    mad_flags = if ("mad" %in% method)
      as.integer(tapply(out$mad_flag, out$group, sum, na.rm = TRUE)) else NA,
    z_flags = if ("zscore" %in% method)
      as.integer(tapply(out$z_flag, out$group, sum, na.rm = TRUE)) else NA,
    robust_z_flags = if ("robust_z" %in% method)
      as.integer(tapply(out$robust_z_flag, out$group, sum, na.rm = TRUE)) else NA
  )

  .crossstatr_result(
    list(flags = out, counts = counts, method = "Outlier diagnostics",
         call = call),
    "crossstatr_outliers")
}

#' @export
print.crossstatr_outliers <- function(x, ...) {
  cat("<crossStatR::outliers> Flag counts per group:\n")
  print(x$counts, row.names = FALSE)
  cat("\nFlagged observations (any method):\n")
  flagged <- x$flags[
    rowSums(x$flags[, c("iqr_flag","mad_flag","z_flag","robust_z_flag")],
            na.rm = TRUE) > 0, ]
  if (nrow(flagged) > 0L) print(flagged, row.names = FALSE) else
    cat("  (none)\n")
  invisible(x)
}
