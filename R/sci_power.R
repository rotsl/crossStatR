#' Power and sample size for common designs
#'
#' Analytic power (or sample size) for the workhorse inferential designs
#' used in experimental work: two-sample t, paired t, correlation,
#' one-way ANOVA, and multiple regression (overall F test). When
#' `n_sim > 0`, the analytic power is cross-checked (or replaced) by a
#' simulation from the assumed generating model. The documentation recommends reporting the
#' simulation whenever the effect size is small or the design has
#' unequal group sizes, where the analytic formulae are known to be slightly
#' anti-conservative.
#'
#' @param design character; one of `"two_sample_t"`, `"paired_t"`,
#'   `"correlation"`, `"anova_oneway"`, `"regression"`.
#' @param effect_size numeric; Cohen's d for t-tests, the Pearson r for
#'   correlation, Cohen's f for one-way ANOVA, or Cohen's f² for regression.
#' @param n integer; sample size. For `two_sample_t` this is the per-group
#'   size when `n1`/`n2` are not given. For `anova_oneway` it is the
#'   per-group size. For `paired_t`, `correlation`, and `regression` it is
#'   the total sample size. If `NULL` and `power` is given, the function
#'   solves for `n`.
#' @param n1,n2 optional per-group sizes for `two_sample_t` (unequal-n
#'   case). If supplied, they override `n`.
#' @param alpha type-I error rate (default 0.05).
#' @param power target power (for sample-size determination). If `NULL`,
#'   power is computed from `n`.
#' @param alternative `"two.sided"` (default) or `"one.sided"`. F-tests in
#'   ANOVA/regression are inherently two-sided; this argument is silently
#'   ignored for those designs.
#' @param n_sim integer; if greater than 0, run `n_sim` simulations and
#'   return the empirical rejection rate in `sim_results`.
#' @param seed optional integer seed for reproducibility of the simulation.
#' @param k integer; number of groups for `anova_oneway` (default 3).
#' @param predictors integer; number of predictors for `regression`
#'   (default 1).
#'
#' @return A `crossstatr_power` object with elements:
#'   \item{design}{design label.}
#'   \item{effect_size}{effect size supplied.}
#'   \item{n}{the sample size (per-group for `two_sample_t`/`anova_oneway`,
#'     total otherwise). For two-sample t with unequal n, this is `n1+n2`.}
#'   \item{power}{computed or target power.}
#'   \item{alpha}{type-I error used.}
#'   \item{method}{short method label.}
#'   \item{sim_results}{NULL or a length-`n_sim` logical vector of rejection
#'     indicators.}
#'   \item{n_sim}{number of simulations run (0 if none).}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   **Two-sample t** (equal variance, equal n):
#'   \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{s_p\sqrt{2/n}}, \qquad
#'         \delta = d\sqrt{n/2}.}
#'   Power is \eqn{1 - F_{t_{2n-2,\delta}}(t_{1-\alpha/2,\,2n-2}) +
#'   F_{t_{2n-2,\delta}}(-t_{1-\alpha/2,\,2n-2})} for two-sided alternatives.
#'   The package delegates to [stats::power.t.test] when n is equal across groups and
#'   compute directly otherwise.
#'
#'   **Paired t**: same machinery with one sample of differences of size `n`
#'   and noncentrality \eqn{\delta = d_z\sqrt{n}}.
#'
#'   **Correlation** (Fisher z): under the alternative \eqn{\rho \neq 0},
#'   \deqn{z_r = \mathrm{atanh}(r), \qquad \lambda = z_r\sqrt{n - 3},}
#'   and the two-sided power is
#'   \deqn{1 - \Phi(z_{1-\alpha/2} - \lambda) + \Phi(-z_{1-\alpha/2} - \lambda).}
#'   For one-sided tests the package uses \eqn{z_{1-\alpha}} and the appropriate tail.
#'
#'   **One-way ANOVA** (k groups, n per group, total `N = n*k`):
#'   Cohen's \eqn{f^2 = \sigma_m^2/\sigma^2} where \eqn{\sigma_m^2} is the
#'   variance of the group means. The noncentrality parameter is
#'   \deqn{\lambda = N f^2 = n k f^2,}
#'   with \eqn{df_1 = k-1}, \eqn{df_2 = N - k}, and
#'   \deqn{\mathrm{power} = 1 - F_{df_1, df_2, \lambda}(F_{1-\alpha, df_1, df_2}).}
#'
#'   **Regression** (overall F test): with `p` predictors and total `n`,
#'   \eqn{f^2 = R^2/(1-R^2)}, noncentrality \eqn{\lambda = n f^2},
#'   \eqn{df_1 = p}, \eqn{df_2 = n - p - 1}, and the power is the upper
#'   tail of the noncentral F at the critical value.
#'
#'   **Simulation**: for `n_sim > 0`, the assumed model is simulated from a
#'   Gaussian with the supplied effect size (Cohen's d for t-tests, Pearson
#'   r for correlation, Cohen's f for ANOVA, Cohen's f² for regression),
#'   the test is run, and the proportion of rejections is reported as the
#'   empirical power. A finite seed makes the simulation reproducible.
#'
#' @section Limitations:
#'
#'   - Analytic formulae assume Gaussian errors, equal variances (for
#'     two-sample t and ANOVA), and the correctness of the supplied effect
#'     size. None of these are guaranteed in practice; the documentation recommends a
#'     sensitivity analysis with `n_sim > 0`.
#'   - The correlation test formula uses the Fisher-z approximation, which
#'     is mildly biased for small `n` (Bonett & Wright, 2000). For
#'     `n < 25`, consider exact permutation power.
#'   - Sample-size determination uses `uniroot` with a bracket of
#'     `[4, 1e6]`; if the target power cannot be reached within that range
#'     (e.g. effect size exactly zero), the function errors.
#'   - Simulation for `anova_oneway` distributes group means evenly; if the
#'     real pattern of means is uneven, the empirical power will differ.
#'
#' @references
#'   - Cohen, J. (1988). Statistical Power Analysis for the Behavioral
#'     Sciences. 2nd ed. Lawrence Erlbaum Associates.
#'     \doi{10.4324/9780203771587}.
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
#' # Power for two-sample t with d = 0.5, n = 30 per group
#' sci_power("two_sample_t", effect_size = 0.5, n = 30)
#'
#' # Sample size for correlation r = 0.3, target power = 0.8
#' sci_power("correlation", effect_size = 0.3, power = 0.8)
#'
#' # ANOVA power with Cohen's f = 0.25, 3 groups, 20 per group
#' sci_power("anova_oneway", effect_size = 0.25, n = 20, k = 3)
#'
#' # Cross-check with simulation
#' \donttest{
#' sci_power("correlation", effect_size = 0.3, n = 100, n_sim = 50, seed = 1)
#' }
#'
#' @export
sci_power <- function(design = c("two_sample_t", "paired_t", "correlation",
                                  "anova_oneway", "regression"),
                       effect_size, n = NULL, n1 = NULL, n2 = NULL,
                       alpha = 0.05, power = NULL,
                       alternative = c("two.sided", "one.sided"),
                       n_sim = 0L, seed = NULL,
                       k = 3L, predictors = 1L) {
  call <- match.call()
  design <- match.arg(design)
  alternative <- match.arg(alternative)
  .assert_probability(alpha, "alpha")
  if (!is.null(power)) .assert_probability(power, "power")
  if (!is.numeric(effect_size) || length(effect_size) != 1L ||
      !is.finite(effect_size)) {
    stop("'effect_size' must be a finite numeric scalar.", call. = FALSE)
  }
  .assert_positive_integer(k, "k")
  .assert_positive_integer(predictors, "predictors")
  n_sim <- as.integer(n_sim)
  if (n_sim < 0L) stop("'n_sim' must be >= 0.", call. = FALSE)

  has_n <- !is.null(n) ||
    (design == "two_sample_t" && !is.null(n1) && !is.null(n2))
  if (!has_n && is.null(power)) {
    stop("Either 'n' (or 'n1'/'n2' for two_sample_t) or 'power' must be supplied.",
         call. = FALSE)
  }
  solving_n <- !has_n && !is.null(power)
  if (n_sim > 0L && !has_n) {
    stop("Simulation requires 'n' to be specified.", call. = FALSE)
  }

  # Dispatch
  if (design == "two_sample_t") {
    res <- .power_two_sample_t(effect_size, n, n1, n2, alpha, power,
                               alternative, solving_n)
  } else if (design == "paired_t") {
    res <- .power_paired_t(effect_size, n, alpha, power, alternative,
                           solving_n)
  } else if (design == "correlation") {
    res <- .power_correlation(effect_size, n, alpha, power, alternative,
                              solving_n)
  } else if (design == "anova_oneway") {
    res <- .power_anova(effect_size, n, k, alpha, power, solving_n)
  } else { # regression
    res <- .power_regression(effect_size, n, predictors, alpha, power,
                             solving_n)
  }

  # Simulation (optional)
  sim_results <- NULL
  if (n_sim > 0L) {
    old_seed <- .set_seed(seed)
    on.exit(.restore_seed(old_seed), add = TRUE)
    sim_results <- .power_simulate(design, effect_size, res$n, n1, n2, alpha,
                                  alternative, k, predictors, n_sim)
    res$power_sim <- mean(sim_results)
  }

  out <- .crossstatr_result(
    list(design = design,
         effect_size = effect_size,
         n = res$n,
         power = res$power,
         power_sim = if (is.null(sim_results)) NULL else mean(sim_results),
         alpha = alpha,
         method = res$method,
         sim_results = sim_results,
         n_sim = n_sim,
         extra = res$extra,
         call = call),
    "crossstatr_power"
  )
  out
}

# ---- Analytic power for each design -----------------------------------------

# @noRd
.power_two_sample_t <- function(d, n, n1, n2, alpha, power, alternative,
                                solving_n) {
  if (!is.null(n1) && !is.null(n2)) {
    # Unequal n: direct formula (equal-variance Student t)
    df <- n1 + n2 - 2
    ncp <- d / sqrt(1 / n1 + 1 / n2)
    p <- .t_power(ncp, df, alpha, alternative)
    return(list(n = n1 + n2, power = p,
                method = "Two-sample t (unequal n)",
                extra = list(n1 = n1, n2 = n2, d = d)))
  }
  if (solving_n) {
    alt <- if (alternative == "two.sided") "two.sided" else "one.sided"
    pt <- stats::power.t.test(delta = d, sd = 1, sig.level = alpha,
                             power = power, type = "two.sample",
                             alternative = alt)
    n_per <- ceiling(pt$n)
    return(list(n = n_per, power = power,
                method = "Two-sample t (solve n)", extra = list(d = d)))
  }
  if (!is.null(n1)) n <- n1
  .assert_positive_integer(n, "n")
  alt <- if (alternative == "two.sided") "two.sided" else "one.sided"
  pt <- stats::power.t.test(n = n, delta = d, sd = 1, sig.level = alpha,
                           type = "two.sample", alternative = alt)
  list(n = n, power = pt$power, method = "Two-sample t",
       extra = list(d = d))
}

# @noRd
.power_paired_t <- function(d, n, alpha, power, alternative, solving_n) {
  if (solving_n) {
    alt <- if (alternative == "two.sided") "two.sided" else "one.sided"
    pt <- stats::power.t.test(delta = d, sd = 1, sig.level = alpha,
                             power = power, type = "paired",
                             alternative = alt)
    n_p <- ceiling(pt$n)
    return(list(n = n_p, power = power, method = "Paired t (solve n)",
                extra = list(d = d)))
  }
  .assert_positive_integer(n, "n")
  alt <- if (alternative == "two.sided") "two.sided" else "one.sided"
  pt <- stats::power.t.test(n = n, delta = d, sd = 1, sig.level = alpha,
                           type = "paired", alternative = alt)
  list(n = n, power = pt$power, method = "Paired t", extra = list(d = d))
}

# @noRd
.power_correlation <- function(r, n, alpha, power, alternative, solving_n) {
  if (abs(r) >= 1) stop("Effect size r must satisfy |r| < 1.", call. = FALSE)
  if (solving_n) {
    f <- function(nn) {
      .corr_power(r, nn, alpha, alternative) - power
    }
    sol <- tryCatch(stats::uniroot(f, interval = c(4, 1e6), extendInt = "upX"),
                    error = function(e) {
      stop("Could not solve for n; check effect size and target power.",
           call. = FALSE)
    })
    n_sol <- ceiling(sol$root)
    return(list(n = n_sol, power = power, method = "Correlation (solve n)",
                extra = list(r = r)))
  }
  .assert_positive_integer(n, "n")
  if (n < 4L) stop("Correlation power requires n >= 4.", call. = FALSE)
  p <- .corr_power(r, n, alpha, alternative)
  list(n = n, power = p, method = "Correlation (Fisher z)",
       extra = list(r = r))
}

# @noRd
.power_anova <- function(f_eff, n, k, alpha, power, solving_n) {
  if (f_eff < 0) stop("Cohen's f must be >= 0.", call. = FALSE)
  df1 <- k - 1
  if (solving_n) {
    ff <- function(nn) {
      N <- nn * k
      df2 <- N - k
      ncp <- N * f_eff^2
      Fcrit <- stats::qf(1 - alpha, df1, df2)
      1 - stats::pf(Fcrit, df1, df2, ncp) - power
    }
    sol <- tryCatch(stats::uniroot(ff, interval = c(2, 1e6),
                                  extendInt = "upX"),
                    error = function(e) {
      stop("Could not solve for n; check effect size and target power.",
           call. = FALSE)
    })
    n_sol <- ceiling(sol$root)
    return(list(n = n_sol, power = power,
                method = "ANOVA one-way (solve n)",
                extra = list(f = f_eff, k = k)))
  }
  .assert_positive_integer(n, "n")
  N <- n * k
  df2 <- N - k
  ncp <- N * f_eff^2
  Fcrit <- stats::qf(1 - alpha, df1, df2)
  p <- 1 - stats::pf(Fcrit, df1, df2, ncp)
  list(n = n, power = p, method = "ANOVA one-way",
       extra = list(f = f_eff, k = k, N = N, df1 = df1, df2 = df2,
                   ncp = ncp))
}

# @noRd
.power_regression <- function(f2, n, p_pred, alpha, power, solving_n) {
  if (f2 < 0) stop("Cohen's f^2 must be >= 0.", call. = FALSE)
  df1 <- p_pred
  if (solving_n) {
    ff <- function(nn) {
      df2 <- nn - p_pred - 1
      if (df2 < 1) return(-1)
      ncp <- nn * f2
      Fcrit <- stats::qf(1 - alpha, df1, df2)
      1 - stats::pf(Fcrit, df1, df2, ncp) - power
    }
    sol <- tryCatch(stats::uniroot(ff,
                                  interval = c(p_pred + 2, 1e6),
                                  extendInt = "upX"),
                    error = function(e) {
      stop("Could not solve for n; check effect size and target power.",
           call. = FALSE)
    })
    n_sol <- ceiling(sol$root)
    return(list(n = n_sol, power = power,
                method = "Regression overall F (solve n)",
                extra = list(f2 = f2, p = p_pred)))
  }
  .assert_positive_integer(n, "n")
  df2 <- n - p_pred - 1
  if (df2 < 1) stop("n must be > predictors + 1.", call. = FALSE)
  ncp <- n * f2
  Fcrit <- stats::qf(1 - alpha, df1, df2)
  p <- 1 - stats::pf(Fcrit, df1, df2, ncp)
  list(n = n, power = p, method = "Regression overall F",
       extra = list(f2 = f2, p = p_pred, df1 = df1, df2 = df2, ncp = ncp))
}

# Power of a one-sample / paired t given ncp, df, alpha, alternative.
# @noRd
.t_power <- function(ncp, df, alpha, alternative) {
  if (alternative == "two.sided") {
    tcrit <- stats::qt(1 - alpha / 2, df)
    1 - stats::pt(tcrit, df, ncp) + stats::pt(-tcrit, df, ncp)
  } else if (alternative == "greater") {
    tcrit <- stats::qt(1 - alpha, df)
    1 - stats::pt(tcrit, df, ncp)
  } else {
    tcrit <- stats::qt(1 - alpha, df)
    stats::pt(-tcrit, df, ncp)
  }
}

# Power for the correlation test based on Fisher z.
# @noRd
.corr_power <- function(r, n, alpha, alternative) {
  zr <- atanh(r)
  ncp <- zr * sqrt(n - 3)
  if (alternative == "two.sided") {
    zc <- stats::qnorm(1 - alpha / 2)
    1 - stats::pnorm(zc - ncp) + stats::pnorm(-zc - ncp)
  } else if (alternative == "greater") {
    zc <- stats::qnorm(1 - alpha)
    1 - stats::pnorm(zc - ncp)
  } else {
    zc <- stats::qnorm(1 - alpha)
    stats::pnorm(-zc - ncp)
  }
}

# ---- Simulation-based power ------------------------------------------------

# @noRd
.power_simulate <- function(design, effect_size, n, n1, n2, alpha,
                             alternative, k, predictors, n_sim) {
  rejs <- logical(n_sim)
  for (b in seq_len(n_sim)) {
    rejs[b] <- .one_sim(design, effect_size, n, n1, n2, alpha, alternative,
                        k, predictors)
  }
  rejs
}

# @noRd
.one_sim <- function(design, eff, n, n1, n2, alpha, alternative, k, predictors) {
  if (design == "two_sample_t") {
    if (!is.null(n1) && !is.null(n2)) {
      a <- n1; b <- n2
    } else {
      a <- b <- n
    }
    x <- stats::rnorm(a, mean = eff)
    y <- stats::rnorm(b, mean = 0)
    tt <- stats::t.test(x, y, alternative = alternative,
                       var.equal = TRUE)
    tt$p.value < alpha
  } else if (design == "paired_t") {
    d <- stats::rnorm(n, mean = eff)
    x <- d + stats::rnorm(n, sd = 0.5)
    y <- x - d
    tt <- stats::t.test(x, y, paired = TRUE, alternative = alternative)
    tt$p.value < alpha
  } else if (design == "correlation") {
    # Bivariate normal via Cholesky
    z1 <- stats::rnorm(n)
    z2 <- stats::rnorm(n)
    x <- z1
    y <- eff * z1 + sqrt(1 - eff^2) * z2
    ct <- stats::cor.test(x, y, alternative = alternative)
    ct$p.value < alpha
  } else if (design == "anova_oneway") {
    means_unscaled <- seq(-1, 1, length.out = k)
    vpop <- mean(means_unscaled^2)
    a <- eff / sqrt(vpop)
    means <- a * means_unscaled
    y <- unlist(lapply(seq_len(k), function(j) {
      stats::rnorm(n, mean = means[j])
    }), use.names = FALSE)
    g <- factor(rep(seq_len(k), each = n))
    aov_fit <- stats::oneway.test(y ~ g, var.equal = TRUE)
    aov_fit$p.value < alpha
  } else { # regression
    p <- predictors
    X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    beta <- c(rep(sqrt(eff / p), p))  # gives f^2 = sum(beta^2) = eff approximately
    y <- as.numeric(X %*% beta) + stats::rnorm(n)
    fit <- stats::lm(y ~ X)
    fstat <- summary(fit)$fstatistic
    1 - stats::pf(fstat[1L], fstat[2L], fstat[3L]) < alpha
  }
}

#' @export
print.crossstatr_power <- function(x, ...) {
  cat(sprintf("<crossStatR::power> design = %s\n", x$design))
  cat(sprintf("effect_size = %.4f,  n = %s,  alpha = %.3f\n",
              x$effect_size,
              if (is.null(x$n)) "solved" else as.character(x$n),
              x$alpha))
  cat(sprintf("power (analytic) = %.3f\n", x$power))
  if (!is.null(x$power_sim)) {
    cat(sprintf("power (simulation, n_sim = %d) = %.3f\n",
                x$n_sim, x$power_sim))
  }
  cat(sprintf("method: %s\n", x$method))
  invisible(x)
}
