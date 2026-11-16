#' Internal helpers and shared infrastructure for crossStatR
#'
#' @importFrom stats approx complete.cases dist model.matrix na.pass rnorm setNames var
#' @importFrom utils str
#' @keywords internal
#' @noRd
NULL

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

.assert_numeric_vector <- function(x, name = "x", min_length = 1L,
                                   finite = TRUE, allow_na = FALSE) {
  if (is.null(x)) {
    stop(sprintf("'%s' must not be NULL.", name), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be a numeric vector; got %s.",
                 name, class(x)[1L]), call. = FALSE)
  }
  if (length(x) < min_length) {
    stop(sprintf("'%s' must have length >= %d (got %d).",
                 name, min_length, length(x)), call. = FALSE)
  }
  if (finite && !allow_na) {
    if (any(!is.na(x) & !is.finite(x))) {
      bad <- which(!is.na(x) & !is.finite(x))[1L]
      stop(sprintf("'%s' contains non-finite values at position %d.",
                   name, bad), call. = FALSE)
    }
  }
  invisible(TRUE)
}

.assert_data_frame <- function(x, name = "data", required_cols = character(0)) {
  if (is.null(x)) {
    stop(sprintf("'%s' must not be NULL.", name), call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop(sprintf("'%s' must be a data.frame; got %s.", name, class(x)[1L]), call. = FALSE)
  }
  if (nrow(x) == 0L) {
    stop(sprintf("'%s' has zero rows.", name), call. = FALSE)
  }
  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(sprintf("'%s' must contain columns: %s (missing: %s).",
                 name, paste(required_cols, collapse = ", "),
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.assert_positive_integer <- function(x, name = "n") {
  if (length(x) != 1L || !is.numeric(x) || x <= 0 ||
      !is.finite(x) || x != round(x)) {
    stop(sprintf("'%s' must be a single positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.assert_probability <- function(x, name = "p") {
  if (length(x) != 1L || !is.numeric(x) ||
      !is.finite(x) || x <= 0 || x >= 1) {
    stop(sprintf("'%s' must be a single probability in (0, 1).", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_flag <- function(x, name = "flag") {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be a single TRUE/FALSE.", name), call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------

.safe_sd <- function(x, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  stats::sd(x)
}

.safe_var <- function(x, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  stats::var(x)
}

.safe_se <- function(x, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  stats::sd(x) / sqrt(n)
}

.safe_mad <- function(x, constant = 1.4826, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::mad(x, constant = constant)
}

.ci_t <- function(x, conf_level = 0.95, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(c(NA_real_, NA_real_))
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(c(m, m))
  }
  alpha <- 1 - conf_level
  tcrit <- stats::qt(1 - alpha / 2, df = n - 1)
  se <- s / sqrt(n)
  c(m - tcrit * se, m + tcrit * se)
}

.ci_median <- function(x, conf_level = 0.95, na_rm = TRUE) {
  if (na_rm) x <- x[is.finite(x)]
  n <- length(x)
  if (n < 5L) return(c(NA_real_, NA_real_))
  alpha <- 1 - conf_level
  # Order-statistic-based interval (normal approximation to binomial order-statistic bounds)
  z <- stats::qnorm(1 - alpha / 2)
  k <- (n - 1) / 2 - z * sqrt(n) / 2
  k <- max(1L, floor(k))
  sx <- sort(x)
  lo_idx <- k
  hi_idx <- n - k + 1L
  lo_idx <- min(max(lo_idx, 1L), n)
  hi_idx <- min(max(hi_idx, 1L), n)
  c(sx[lo_idx], sx[hi_idx])
}

# Cochran's bias correction for sample variance with replicates (N-1 in denom)
# Used by replicate-aware variance functions.
.cochran_correction <- function(n, sigma2) {
  if (n < 2L) return(NA_real_)
  sigma2 * (1 - 1 / (4 * (n - 1)) + 1 / (32 * (n - 1)^2))
}

# ---------------------------------------------------------------------------
# Standardized effect sizes (shared by A4, A5, A10)
# ---------------------------------------------------------------------------

.cohens_d <- function(x, y = NULL, paired = FALSE, na_rm = TRUE,
                     hedges = FALSE) {
  if (na_rm) {
    x <- x[is.finite(x)]
    if (!is.null(y)) y <- y[is.finite(y)]
  }
  if (is.null(y)) {
    if (length(x) < 2L) return(NA_real_)
    m <- mean(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(NA_real_)
    d <- m / s
    if (hedges) {
      n <- length(x)
      J <- 1 - 3 / (4 * (n - 1) - 1)
      d <- d * J
    }
    return(d)
  }
  nx <- length(x); ny <- length(y)
  if (nx < 2L || ny < 2L) return(NA_real_)
  if (paired) {
    if (nx != ny) {
      stop("paired effect size requires equal-length samples.",
           call. = FALSE)
    }
    d <- mean(x - y) / stats::sd(x - y)
  } else {
    # Pooled SD, equal-variance formulation (Cohen 1988, eq 2.42)
    sp <- sqrt(((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) /
                 (nx + ny - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    d <- (mean(x) - mean(y)) / sp
  }
  if (hedges) {
    J <- 1 - 3 / (4 * (nx + ny - 2) - 1)
    d <- d * J
  }
  d
}

# Approximate SE of Cohen's d using the delta method (Hedges & Olkin 1985)
.se_cohens_d <- function(d, n1, n2 = NULL, paired = FALSE) {
  if (is.null(n2) || paired) {
    if (is.null(n2)) n2 <- n1
    n <- n1
    if (n < 3L) return(NA_real_)
    sqrt((1 / n + d^2 / (2 * n)) * (1 + 1 / (n - 1)))
  } else {
    neff <- (n1 * n2) / (n1 + n2)
    sqrt(((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))) *
           (1 + 1 / (n1 + n2 - 3)))
  }
}

# ---------------------------------------------------------------------------
# Cluster bootstrap (used by sci_bootstrap)
# ---------------------------------------------------------------------------

.resample_clusters <- function(ids, rng = NULL) {
  if (!is.null(rng)) {
    old <- .Random.seed
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    assign(".Random.seed", rng, envir = globalenv())
  }
  levs <- unique(ids)
  sampled <- sample(levs, size = length(levs), replace = TRUE)
  # Expand back to observation level
  unlist(lapply(seq_along(sampled), function(k) {
    which(ids == sampled[k])
  }), use.names = FALSE)
}

# ---------------------------------------------------------------------------
# S3 result class
# ---------------------------------------------------------------------------

#' Construct a crossStatR result object
#'
#' @param x named list of components
#' @param class subclass name
#' @return crossStatR S3 object
#' @noRd
.crossstatr_result <- function(x, class) {
  structure(x, class = c(class, "crossstatr_result"))
}

#' Generic: tidy a crossStatR result
#' @param x result object
#' @param ... additional arguments
#' @export
tidy.crossstatr_result <- function(x, ...) {
  # Default: look for a 'tidy' element or coerce to data.frame
  if (!is.null(x$tidy)) {
    return(as.data.frame(x$tidy))
  }
  if (!is.null(x$coefficients)) {
    return(as.data.frame(x$coefficients))
  }
  if (!is.null(x$summary)) {
    return(as.data.frame(x$summary))
  }
  as.data.frame(x[!vapply(x, is.function, logical(1))], stringsAsFactors = FALSE)
}

#' @export
as.data.frame.crossstatr_result <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  tidy.crossstatr_result(x)
}

#' @export
print.crossstatr_result <- function(x, ...) {
  cat(sprintf("<%s>\n", class(x)[1L]))
  if (!is.null(x$method)) cat("Method:", x$method, "\n")
  if (!is.null(x$summary)) {
    cat("Summary:\n")
    print(x$summary)
  } else if (!is.null(x$tidy)) {
    cat("Tidy:\n")
    print(x$tidy)
  } else {
    str(x, max.level = 2L, give.attr = FALSE, no.list = TRUE)
  }
  invisible(x)
}

#' Generic: summary of crossStatR result
#' @param object result object
#' @param ... additional arguments
#' @export
summary.crossstatr_result <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

# ---------------------------------------------------------------------------
# Plot helper - simple base-R fallback if ggplot2 not available at runtime
# ---------------------------------------------------------------------------

.ensure_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the 'ggplot2' package.", call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Time-grid helpers
# ---------------------------------------------------------------------------

.normalize_grid <- function(x, y) {
  .assert_numeric_vector(x, "x", min_length = 2L)
  .assert_numeric_vector(y, "y", min_length = 2L)
  if (length(x) != length(y)) {
    stop("x and y must have equal length.", call. = FALSE)
  }
  list(x = x, y = y)
}

# Linear interpolation onto a common grid, with extrapolation disallowed.
.interp_to_grid <- function(x, y, target_x) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) {
    stop("Need at least 2 finite points to interpolate.", call. = FALSE)
  }
  rng <- range(x)
  tgt_in <- target_x[target_x >= rng[1] & target_x <= rng[2]]
  if (length(tgt_in) < 2L) {
    stop("Common grid does not overlap with signal range.", call. = FALSE)
  }
  stats::approx(x, y, xout = tgt_in, rule = 2)$y
}

# ---------------------------------------------------------------------------
# Bootstrap CI helper (used by many functions)
# ---------------------------------------------------------------------------

.boot_ci <- function(stats_boot, conf_level = 0.95, type = "perc") {
  alpha <- 1 - conf_level
  qs <- stats::quantile(stats_boot, probs = c(alpha / 2, 1 - alpha / 2),
                        na.rm = TRUE)
  names(qs) <- NULL
  qs
}

# ---------------------------------------------------------------------------
# Numerical Jacobian via finite differences (used by parameter_uncertainty
# and the delta method for arbitrary functions of fitted parameters).
# ---------------------------------------------------------------------------

.num_jacobian <- function(par, fn, eps = 1e-6) {
  p <- length(par)
  f0 <- fn(par)
  out <- matrix(NA_real_, nrow = length(f0), ncol = p)
  for (k in seq_len(p)) {
    p_plus <- p_minus <- par
    p_plus[k] <- p_plus[k] + eps
    p_minus[k] <- p_minus[k] - eps
    out[, k] <- (fn(p_plus) - fn(p_minus)) / (2 * eps)
  }
  out
}

# ---------------------------------------------------------------------------
# Multiple comparison helper
# ---------------------------------------------------------------------------

.adjust_p <- function(p, method = "holm") {
  if (length(p) == 0L) return(structure(numeric(0), names = character(0)))
  stats::p.adjust(p, method = method)
}

# ---------------------------------------------------------------------------
# Detection-limit math (used by robust_lod, detection_reliability)
# ---------------------------------------------------------------------------

# Inverse-predict a concentration from a calibration curve y = a + b*x
# (Currie 1968 / IUPAC 1995 simple form).  Returns list(xhat, se_xhat, ci).
.inverse_predict_linear <- function(cal, y0, conf_level = 0.95,
                                    se_resid = NULL, df_resid = NULL) {
  a <- cal$coef[1L]; b <- cal$coef[2L]
  if (abs(b) < .Machine$double.eps) {
    stop("Calibration slope is zero; cannot inverse-predict.", call. = FALSE)
  }
  x0 <- (y0 - a) / b
  if (is.null(se_resid) || is.null(df_resid)) {
    return(list(xhat = x0, se_xhat = NA_real_, ci = c(NA_real_, NA_real_)))
  }
  if (is.null(cal$x)) {
    return(list(xhat = x0, se_xhat = NA_real_, ci = c(NA_real_, NA_real_)))
  }
  n <- length(cal$x)
  xm <- mean(cal$x)
  ssx <- sum((cal$x - xm)^2)
  tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = df_resid)
  se_xhat <- (se_resid / abs(b)) *
    sqrt(1 + 1 / n + (x0 - xm)^2 / ssx)
  ci <- c(x0 - tcrit * se_xhat, x0 + tcrit * se_xhat)
  list(xhat = x0, se_xhat = se_xhat, ci = ci)
}

# ---------------------------------------------------------------------------
# Random seed utilities
# ---------------------------------------------------------------------------

.set_seed <- function(seed = NULL) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv())) {
      set.seed(seed)
    } else {
      old_seed <- .Random.seed
      set.seed(seed)
      return(invisible(old_seed))
    }
  }
  invisible(NULL)
}

.restore_seed <- function(old_seed) {
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Monte Carlo standard error helpers
# ---------------------------------------------------------------------------

# MCSE for a probability estimate p_hat based on N simulations
.mcse_p <- function(p_hat, N) {
  sqrt(p_hat * (1 - p_hat) / N)
}

# MCSE for a mean estimate (e.g. bias, MSE) using sample SD of replicates
.mcse_mean <- function(replicates) {
  n <- length(replicates)
  if (n < 2L) return(NA_real_)
  stats::sd(replicates) / sqrt(n)
}

utils::globalVariables(c("theoretical", "sample", ".crossstatr_weights"))

# Symbols used in ggplot2 data masks.
utils::globalVariables(c("theoretical", "sample", ".crossstatr_weights"))
