#' Descriptive statistics for scientific experiments
#'
#' Computes a comprehensive set of descriptive statistics tailored to
#' experimental data, with explicit support for technical and biological
#' replicates. Technical replicates are *not* treated as independent
#' biological observations: when a `bio_rep` grouping is supplied,
#' technical replicates are summarized to a per-biological-unit value
#' before the final statistics are computed.
#'
#' @param x numeric vector of observations.
#' @param group optional factor or coercible vector giving experimental
#'   groups (e.g. treatment levels). If supplied, statistics are computed
#'   per group.
#' @param bio_rep optional vector identifying biological replicates.
#'   When supplied together with `tech_rep`, the function first averages
#'   technical replicates within each `bio_rep` x `group` cell, then
#'   summarises biological replicates. This prevents pseudo-replication.
#' @param tech_rep optional vector identifying technical replicates.
#' @param units character scalar; measurement units for `x`. Stored on
#'   the result for downstream labelling; not used in computation.
#' @param conf_level confidence level for the mean and median intervals
#'   (default 0.95).
#' @param na.rm logical; remove non-finite observations before computation?
#' @param quantile_probs numeric vector of probabilities for quantile
#'   reporting (default `c(0, 0.25, 0.5, 0.75, 1)`).
#'
#' @return A `crossstatr_describe` object with elements:
#'   \item{summary}{data frame of per-group (or overall) statistics.}
#'   \item{n_obs, n_bio_rep, n_tech_rep}{observation and replicate counts.}
#'   \item{units}{passed-through measurement units.}
#'   \item{method}{method label.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   The mean and variance use the standard unbiased estimators with
#'   Bessel's correction:
#'
#'   \deqn{\bar{x} = \frac{1}{n}\sum_{i=1}^{n} x_i, \qquad
#'         s^2 = \frac{1}{n-1}\sum_{i=1}^{n}(x_i - \bar{x})^2.}
#'
#'   The standard error is \eqn{s / \sqrt{n}}. The confidence interval
#'   for the mean uses a Student-t critical value with `n - 1` df:
#'
#'   \deqn{\bar{x} \pm t_{1-\alpha/2, n-1}\, s / \sqrt{n}.}
#'
#'   The median interval uses order-statistic endpoints with a normal
#'   approximation to binomial bounds:
#'
#'   \deqn{\bigl[x_{(k)}, x_{(n-k+1)}\bigr], \qquad
#'         k = \frac{n-1}{2} - z_{1-\alpha/2}\,\frac{\sqrt{n}}{2}.}
#'
#'   The coefficient of variation is \eqn{s / \bar{x}} (reported as a
#'   proportion, not a percentage). The MAD uses the standard
#'   consistency constant 1.4826 so that it estimates the same scale
#'   parameter as `sd` for Gaussian data.
#'
#'   When `bio_rep` is supplied, per-biological-replicate means are
#'   computed first by averaging technical replicates; this is the
#'   recommended workflow for experimental data and matches the analysis
#'   of Vaux et al. (2012, *EMBO reports*) on avoiding
#'   pseudo-replication. The reported `n` then equals the number of
#'   biological replicates, not the total number of measurements.
#'
#' @section Limitations:
#'
#'   - Confidence intervals assume approximate independence of the
#'     biological replicates. If the biological units are themselves
#'     clustered (e.g. mice within a cage), use [replicate_decompose()]
#'     instead.
#'   - The median CI requires at least 5 observations; with fewer,
#'     it is reported as `NA`.
#'   - The CV is meaningless for signed data near zero; the function
#'     still returns it but emits a warning when the mean changes sign
#'     within a group.
#'
#' @references
#'   - Vaux, D. L., Fidler, F. and Cumming, G. (2012). Replicates and
#'     repeats—what is the difference and is it significant? The EMBO
#'     Reports, 13, 291-296.
#'     \doi{10.1038/embor.2012.36}.
#'
#'   - Hettmansperger, T. P. and McKean, J. W. (2010). Robust Nonparametric
#'     Statistical Methods. 2nd ed. CRC Press.
#'     \doi{10.1201/b10451}.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, mean = 10), rnorm(20, mean = 12))
#' g <- rep(c("A", "B"), each = 20)
#' sci_describe(x, group = g)
#'
#' # With biological/technical replicates
#' df <- data.frame(
#'   y = rnorm(60, mean = rep(c(8, 10, 12), each = 20)),
#'   grp = rep(c("A","B","C"), each = 20),
#'   bio = rep(1:6, each = 10),
#'   tech = rep(1:10, times = 6))
#' sci_describe(df$y, group = df$grp, bio_rep = df$bio, tech_rep = df$tech)
#'
#' @export
sci_describe <- function(x, group = NULL, bio_rep = NULL, tech_rep = NULL,
                        units = "", conf_level = 0.95, na.rm = TRUE,
                        quantile_probs = c(0, 0.25, 0.5, 0.75, 1)) {
  call <- match.call()
  .assert_numeric_vector(x, "x", min_length = 1L, finite = FALSE)
  .assert_probability(conf_level, "conf_level")
  .assert_flag(na.rm, "na.rm")

  if (length(x) != length(group) && !is.null(group)) {
    stop("Length of 'group' must match 'x'.", call. = FALSE)
  }
  if (length(x) != length(bio_rep) && !is.null(bio_rep)) {
    stop("Length of 'bio_rep' must match 'x'.", call. = FALSE)
  }
  if (length(x) != length(tech_rep) && !is.null(tech_rep)) {
    stop("Length of 'tech_rep' must match 'x'.", call. = FALSE)
  }

  n_total <- length(x)
  n_bio_rep_total <- if (is.null(bio_rep)) NA_integer_ else
    length(unique(bio_rep))
  n_tech_per_bio <- if (is.null(tech_rep) || is.null(bio_rep)) NA_integer_ else
    as.integer(round(n_total / n_bio_rep_total))

  # Reduce to biological-replicate level when tech_rep info is given
  if (!is.null(bio_rep) && !is.null(tech_rep)) {
    df0 <- data.frame(x = x, group = if (is.null(group)) 1L else group,
                      bio = bio_rep, tech = tech_rep)
    if (na.rm) df0 <- df0[is.finite(df0$x), , drop = FALSE]
    df0 <- stats::aggregate(x ~ group + bio, data = df0, FUN = mean,
                            na.action = na.pass)
    x <- df0$x
    group <- if (is.null(group)) NULL else df0$group
  } else if (!is.null(bio_rep)) {
    df0 <- data.frame(x = x, group = if (is.null(group)) 1L else group,
                      bio = bio_rep)
    if (na.rm) df0 <- df0[is.finite(df0$x), , drop = FALSE]
    df0 <- stats::aggregate(x ~ group + bio, data = df0, FUN = mean,
                            na.action = na.pass)
    x <- df0$x
    group <- if (is.null(group)) NULL else df0$group
  }
  if (is.null(group)) group <- factor(rep(1L, length(x)))

  groups <- if (is.factor(group)) levels(droplevels(group)) else
    sort(unique(group))
  groups <- as.character(groups)

  out <- do.call(rbind, lapply(groups, function(g) {
    xi <- x[group == g]
    if (na.rm) xi <- xi[is.finite(xi)]
    n <- length(xi)
    if (n == 0L) {
      return(data.frame(group = g, n = 0L, mean = NA_real_, sd = NA_real_,
                        var = NA_real_, se = NA_real_,
                        median = NA_real_, mad = NA_real_,
                        cv = NA_real_, iqr = NA_real_,
                        mean_ci_low = NA_real_, mean_ci_high = NA_real_,
                        median_ci_low = NA_real_, median_ci_high = NA_real_,
                        n_missing = 0L))
    }
    m <- mean(xi)
    s <- if (n >= 2L) stats::sd(xi) else NA_real_
    v <- if (n >= 2L) stats::var(xi) else NA_real_
    se <- if (n >= 2L) s / sqrt(n) else NA_real_
    med <- stats::median(xi)
    md <- .safe_mad(xi)
    cv <- if (is.finite(m) && m != 0) s / m else NA_real_
    if (!is.na(cv) && sign(m) != sign(m + s)) {
      warning(sprintf("Group '%s': CV near zero-mean; reported value is unstable.",
                      g), call. = FALSE)
    }
    qs <- stats::quantile(xi, probs = quantile_probs, na.rm = TRUE,
                          names = TRUE)
    iqr <- stats::IQR(xi, na.rm = TRUE)
    mean_ci <- .ci_t(xi, conf_level = conf_level, na_rm = na.rm)
    med_ci <- .ci_median(xi, conf_level = conf_level, na_rm = na.rm)
    n_missing <- length(xi) - sum(is.finite(xi)) + if (na.rm)
      sum(!is.finite(x[group == g])) else 0L
    data.frame(
      group = g, n = n, mean = m, sd = s, var = v, se = se,
      median = med, mad = md, cv = cv, iqr = iqr,
      mean_ci_low = mean_ci[1L], mean_ci_high = mean_ci[2L],
      median_ci_low = med_ci[1L], median_ci_high = med_ci[2L],
      n_missing = as.integer(n_missing),
      stringsAsFactors = FALSE)
  }))
  colnames(out)[1L] <- if (!is.null(bio_rep)) "group" else "group"

  .crossstatr_result(
    list(summary = out,
         n_obs = n_total,
         n_bio_rep = n_bio_rep_total,
         n_tech_per_bio = n_tech_per_bio,
         units = units,
         method = if (!is.null(bio_rep))
           "Descriptive statistics (bio-replicate level)" else
           "Descriptive statistics",
         call = call),
    "crossstatr_describe")
}

#' @export
print.crossstatr_describe <- function(x, ...) {
  cat(sprintf("<crossStatR::describe> %s\n", x$method))
  cat(sprintf("Total observations: %d  Biological replicates: %s\n",
              x$n_obs,
              if (is.na(x$n_bio_rep)) "n/a" else as.character(x$n_bio_rep)))
  if (x$units != "") cat(sprintf("Units: %s\n", x$units))
  print(x$summary, row.names = FALSE)
  invisible(x)
}
