#' Adaptive baseline estimation
#'
#' Estimates a slowly varying baseline under a signal that contains sharp
#' peaks. Two methods are provided. Asymmetric least squares (Eilers &
#' Boelens, 2005) iteratively reweights a penalised least squares fit so the
#' baseline tracks the lower envelope of the signal; the lowess baseline
#' smooths the signal and then clips the result to lie below the data.
#'
#' The documentation recommends asymmetric least squares for chromatography, Raman, and XPS
#' data where the peak structure is sharp relative to the baseline. The
#' lowess option is intended as a sanity check; it does not try to be
#' asymmetric and will be biased upward when the signal sits well above
#' the true baseline.
#'
#' @param x numeric vector of abscissae (wavelength, time, m/z, etc.).
#' @param y numeric vector of signal intensities, same length as `x`.
#' @param method character; one of `"asls"` (asymmetric least squares,
#'   default) or `"lowess"`.
#' @param lambda numeric; smoothness penalty for `"asls"`. Larger values
#'   give smoother baselines. Typical values range from 1e5 to 1e9.
#' @param p numeric in (0, 1); asymmetry parameter for `"asls"`. Smaller
#'   values push the baseline below the signal. The default 0.01 produces
#'   baselines that follow the lower envelope.
#' @param span numeric in (0, 1]; span for `"lowess"`.
#' @param max_iter integer; maximum number of reweighting iterations.
#' @param tol numeric; convergence tolerance on the change in the baseline
#'   between iterations.
#' @param conf_level currently unused; reserved for future uncertainty
#'   quantification on the baseline. Validated for being a probability.
#'
#' @return A `crossstatr_baseline` object with elements:
#'   \item{x}{input abscissae.}
#'   \item{y}{input signal.}
#'   \item{baseline}{fitted baseline.}
#'   \item{corrected}{`y - baseline`.}
#'   \item{weights}{final asymmetric weights (uniform for `"lowess"`).}
#'   \item{converged}{logical; whether reweighting converged within
#'     `max_iter`. For `"lowess"` this is always `TRUE`.}
#'   \item{n_iter}{number of iterations performed (1 for `"lowess"`).}
#'   \item{method}{baseline method used.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Asymmetric least squares (Eilers & Boelens, 2005) minimises the
#'   penalised least squares objective
#'   \deqn{Q(\mathbf{z}) = (\mathbf{y} - \mathbf{z})^\top W
#'         (\mathbf{y} - \mathbf{z}) + \lambda\, \mathbf{z}^\top
#'         D_2^\top D_2\, \mathbf{z},}
#'   where \eqn{D_2} is the \eqn{(n-2) \times n} second-order difference
#'   operator (so \eqn{D_2 \mathbf{z}} measures the discrete curvature of
#'   \eqn{\mathbf{z}}) and \eqn{W = \mathrm{diag}(w_i)} holds the asymmetric
#'   weights \eqn{w_i = p} when \eqn{y_i > z_i} and \eqn{w_i = 1 - p}
#'   otherwise. The closed-form update at each iteration is
#'   \deqn{\mathbf{z} = (W + \lambda D_2^\top D_2)^{-1} W\, \mathbf{y},}
#'   solved via `Matrix::solve`. The package iterates until
#'   \eqn{\max_i |z_i^{\mathrm{new}} - z_i^{\mathrm{old}}| < \mathrm{tol}}
#'   or `max_iter` is reached.
#'
#'   The lowess baseline is `stats::lowess(x, y, f = span)`, aligned back
#'   to the input order, then clipped to the lower envelope by taking
#'   \eqn{\tilde z_i = \min(\hat z_i, y_i)}. This is not truly asymmetric
#'   but enforces the baseline-lies-below-signal property.
#'
#' @section Limitations:
#'
#'   - AsLS weights are updated after each solve; the last reported weights
#'     can differ from the weights used for the last fitted baseline. No
#'     uncertainty band is computed and `conf_level` is reserved for future use.
#'
#'   - ALS with large `lambda` produces highly correlated baselines; small
#'     changes in the signal can shift the whole baseline. Inspect
#'     `corrected` for systematic residual structure before downstream use.
#'   - All methods assume peaks are positive (the signal sits above the
#'     baseline). For negative peaks, multiply the signal by -1 before
#'     fitting or set `p > 0.5`.
#'   - Convergence of the asymmetric reweighting is not guaranteed for
#'     pathological signals. Always check `converged` and `n_iter`.
#'   - The lowess baseline is a poor default for signals with wide peaks;
#'     it will be biased upward and clip the corrected signal at zero.
#'
#' @references
#'   - Eilers, P. H. C. and Boelens, H. F. M. (2005). Baseline Correction
#'     with Asymmetric Least Squares Smoothing. Unpublished manuscript.
#'     \url{https://www.researchgate.net/publication/228961729_Baseline_Correction_with_Asymmetric_Least_Squares_Smoothing}.
#'
#'   - Eilers, P. H. C. (2003). A Perfect Smoother. Analytical Chemistry,
#'     75, 3631-3636.
#'     \doi{10.1021/ac034173t}.
#'
#'   - Cleveland, W. S. (1979). Robust Locally Weighted Regression and
#'     Smoothing Scatterplots. Journal of the American Statistical
#'     Association, 74, 829-836.
#'     \doi{10.1080/01621459.1979.10481038}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 200)
#' peak1 <- 5 * exp(-(x - 3)^2 / 0.1)
#' peak2 <- 3 * exp(-(x - 7)^2 / 0.05)
#' baseline_true <- 0.5 + 0.2 * x + 0.03 * x^2
#' y <- baseline_true + peak1 + peak2 + rnorm(length(x), sd = 0.1)
#' fit <- adaptive_baseline(x, y, method = "asls", lambda = 1e6, p = 0.01)
#' print(fit)
#'
#' fit_l <- adaptive_baseline(x, y, method = "lowess", span = 0.3)
#' print(fit_l)
#'
#' @export
adaptive_baseline <- function(x, y, method = c("asls", "lowess"),
                              lambda = 1e7, p = 0.01, span = 0.5,
                              max_iter = 20L, tol = 1e-6,
                              conf_level = 0.95) {
  call <- match.call()
  method <- match.arg(method)
  .assert_numeric_vector(x, "x", min_length = 5L)
  .assert_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  .assert_probability(p, "p")
  .assert_probability(span, "span")
  .assert_probability(conf_level, "conf_level")
  .assert_positive_integer(max_iter, "max_iter")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a positive finite numeric.", call. = FALSE)
  }
  if (method == "asls") {
    if (!is.numeric(lambda) || length(lambda) != 1L ||
        !is.finite(lambda) || lambda <= 0) {
      stop("'lambda' must be a positive finite numeric.", call. = FALSE)
    }
    fit <- .als_baseline(x, y, lambda = lambda, p = p,
                         max_iter = max_iter, tol = tol)
  } else {
    fit <- .lowess_baseline(x, y, span = span)
  }

  baseline <- fit$baseline
  weights  <- fit$weights
  corrected <- y - baseline

  .crossstatr_result(
    list(
      x = x,
      y = y,
      baseline = baseline,
      corrected = corrected,
      weights = weights,
      converged = fit$converged,
      n_iter = fit$n_iter,
      method = method,
      call = call
    ),
    "crossstatr_baseline"
  )
}

# Asymmetric least squares baseline (Eilers & Boelens 2005).
# Uses a sparse second-order difference matrix D2 and Matrix::solve.
# Convergence on max|z_new - z_old| < tol.
# @noRd
.als_baseline <- function(x, y, lambda, p, max_iter, tol) {
  n <- length(y)
  if (n < 5L) {
    stop("ALS requires at least 5 points.", call. = FALSE)
  }
  D2 <- .build_d2(n)
  DtD <- Matrix::t(D2) %*% D2
  w <- rep(1, n)
  z <- y
  converged <- FALSE
  n_iter <- 0L
  for (iter in seq_len(max_iter)) {
    n_iter <- iter
    W <- Matrix::Diagonal(n, w)
    A <- W + lambda * DtD
    rhs <- as.numeric(w * y)
    z_new <- as.numeric(Matrix::solve(A, rhs))
    if (max(abs(z_new - z)) < tol) {
      z <- z_new
      converged <- TRUE
      break
    }
    z <- z_new
    w <- ifelse(y > z, p, 1 - p)
  }
  list(baseline = z, weights = w, converged = converged, n_iter = n_iter)
}

# Build sparse second-order difference matrix D2 (size (n-2) x n).
# Each row picks (z[i] - 2 z[i+1] + z[i+2]).
# @noRd
.build_d2 <- function(n) {
  if (n < 3L) stop("Need at least 3 points for D2.", call. = FALSE)
  m <- n - 2L
  i <- integer(0L); j <- integer(0L); xv <- numeric(0L)
  for (k in seq_len(m)) {
    i <- c(i, rep(k, 3L))
    j <- c(j, k, k + 1L, k + 2L)
    xv <- c(xv, 1, -2, 1)
  }
  Matrix::sparseMatrix(i = i, j = j, x = xv, dims = c(m, n))
}

# Lowess baseline: smooth, then clip to the lower envelope of y.
# Convergence is trivially TRUE; weights are uniform.
# @noRd
.lowess_baseline <- function(x, y, span) {
  lo <- stats::lowess(x, y, f = span)
  ord <- order(x)
  z <- numeric(length(x))
  z[ord] <- lo$y
  # Lower envelope: baseline must not exceed the signal.
  z <- pmin(z, y)
  list(baseline = z,
       weights = rep(1, length(y)),
       converged = TRUE,
       n_iter = 1L)
}

#' @export
print.crossstatr_baseline <- function(x, ...) {
  cat(sprintf("<crossStatR::adaptive_baseline> method = %s\n", x$method))
  cat(sprintf("n = %d, converged = %s, iterations = %d\n",
              length(x$y), x$converged, x$n_iter))
  cat(sprintf("Baseline range: [%.4g, %.4g]\n",
              min(x$baseline, na.rm = TRUE),
              max(x$baseline, na.rm = TRUE)))
  cat(sprintf("Corrected range: [%.4g, %.4g]\n",
              min(x$corrected, na.rm = TRUE),
              max(x$corrected, na.rm = TRUE)))
  invisible(x)
}
