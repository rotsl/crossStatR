#' Stage-wise uncertainty attribution framework
#'
#' Provides a structured way to attribute the final measurement
#' uncertainty to the various experimental stages: sampling, instrument
#' noise, calibration, preprocessing, replicate variation, model
#' selection, parameter estimation, and drift. The framework is explicit
#' about which attributions are mathematically additive (variance
#' decomposition under independence) and which require a more careful
#' treatment (correlated sources, multiplicative effects).
#'
#' @param stages named list of variance estimates. Names correspond to
#'   stages (e.g. `sampling`, `instrument_noise`, `calibration`,
#'   `replicate`, `drift`). Each entry is a numeric scalar giving the
#'   variance contributed by that stage. Required.
#' @param correlations optional named list of pairwise correlations
#'   between stages. Each entry is a numeric scalar with name like
#'   `"sampling__calibration"` giving the correlation between those two
#'   stages. Defaults to no correlation (full additivity).
#' @param covariances optional covariance matrix, superseding
#'   `correlations` if supplied. Must have `length(stages)` rows/cols.
#' @param model_selection_uncertainty optional numeric scalar giving the
#'   additional variance attributable to model-selection uncertainty
#'   (this is *not* variance decomposition in the strict sense; see
#'   Details).
#' @param n_mc number of Monte Carlo replicates for the joint
#'   propagation under non-zero correlations (default 5000).
#' @param conf_level level for error quantiles and Monte Carlo precision intervals.
#' @param seed optional integer seed.
#'
#' @return A `crossstatr_uncertainty_origin` object.
#'
#' @section Methodology:
#'
#'   **Additive case** (all sources independent): the total variance is
#'   \deqn{\sigma^2_{\mathrm{total}} = \sum_k \sigma^2_k,}
#'   and the contribution of stage \eqn{k} is
#'   \deqn{c_k = \sigma^2_k / \sigma^2_{\mathrm{total}}.}
#'   This is the simplest case and is the only one in which the
#'   attribution is rigorously additive.
#'
#'   **Correlated sources**: if stages \eqn{k} and \eqn{j} are
#'   correlated with \eqn{\rho_{kj}}, the total variance is
#'   \deqn{\sigma^2_{\mathrm{total}} = \sum_k \sigma^2_k + 2 \sum_{k<j}
#'     \rho_{kj}\, \sigma_k\, \sigma_j.}
#'   The implemented allocation uses the cooperative game
#'   \eqn{v(S) = \mathrm{Var}(\sum_{k \in S} X_k)}. Its Shapley
#'   allocation splits each covariance term equally between its stages:
#'   \deqn{\phi_k = \sigma_k^2 + \sum_{j \ne k} \mathrm{Cov}(X_k, X_j).}
#'   This sums to the variance of the additive total. Contributions can
#'   be negative when covariances are negative. It is a different game
#'   from the conditional-expectation Shapley effects studied by Owen and
#'   Prieur (2017) and Song et al. (2016); those papers are related reading,
#'   not a derivation of the implemented allocation.
#'
#'   Monte Carlo draws use a zero-mean multivariate normal distribution
#'   with the supplied covariance matrix and are summed additively.
#'   `total_ci` contains quantiles of that simulated total, not a
#'   confidence interval for its variance or its stage allocations.
#'   `total_prediction_interval` is the same interval with an explicit name.
#'   `total_variance_ci` is a chi-square interval for Monte Carlo variance
#'   precision using `n_mc - 1` df. `shapley_ci` uses Gaussian
#'   sample-covariance standard errors for the allocation. These intervals
#'   describe simulation precision conditional on the supplied covariance,
#'   not uncertainty in stage variances estimated from experimental data.
#'   Nonlinear propagation is not implemented. The optional
#'   `model_selection_uncertainty` value is returned separately without
#'   being added to the simulated total. Draper (1995) discusses model
#'   uncertainty generally; this function does not implement Bayesian
#'   model averaging.
#'
#' @section Limitations:
#'
#'   - The Shapley attribution is well-defined only when the variance
#'     contributions and correlations are themselves estimable. With
#'     noisy estimates of \eqn{\rho_{kj}}, the attribution is unstable.
#'   - Model-selection uncertainty is treated as a separate addendum,
#'     not folded into the variance decomposition. This is intentional;
#'     the alternative (BMA-style folding) requires priors over models.
#'   - For multiplicative uncertainty (e.g. relative standard
#'     uncertainties), the additive decomposition does not apply; use
#'     a separately specified transformation and propagation model.
#'   - The function assumes the stage variances have already been
#'     estimated (e.g. via [batch_effect_quantify()], [drift_uncertainty()],
#'     [noise_characterize()]); it does not estimate them.
#'
#' @references
#'   - Owen, A. B. and Prieur, C. (2017). On Shapley Value for Measuring
#'     Importance of Dependent Inputs. SIAM/ASA Journal on Uncertainty
#'     Quantification, 5, 986-1002.
#'     \doi{10.1137/16M1097717}.
#'
#'   - Song, E., Nelson, B. L. and Staum, J. (2016). Shapley effects for
#'     global sensitivity analysis: Theory and computation. SIAM/ASA
#'     Journal on Uncertainty Quantification, 4, 1060-1083.
#'     \doi{10.1137/15M1048070}.
#'
#'   - Draper, D. (1995). Assessment and Propagation of Model Uncertainty.
#'     Journal of the Royal Statistical Society Series B: Statistical
#'     Methodology, 57, 45-70.
#'     \doi{10.1111/j.2517-6161.1995.tb02015.x}.
#'
#'   - Joint Committee for Guides in Metrology (2008). Evaluation of
#'     measurement data - Guide to the expression of uncertainty in
#'     measurement. JCGM 100:2008. BIPM.
#'     \doi{10.59161/JCGM100-2008E}.
#'
#' @examples
#' stages <- c(sampling = 0.05, instrument = 0.10,
#'             calibration = 0.02, drift = 0.03)
#' res <- uncertainty_origin(stages)
#' print(res)
#'
#' # With correlations
#' res2 <- uncertainty_origin(stages,
#'   correlations = list(sampling__instrument = 0.3,
#'                       calibration__drift = 0.4))
#' print(res2)
#'
#' @export
uncertainty_origin <- function(stages, correlations = NULL,
                                covariances = NULL,
                                model_selection_uncertainty = NULL,
                                n_mc = 5000L, conf_level = 0.95,
                                seed = NULL) {
  call <- match.call()
  .assert_probability(conf_level, "conf_level")

  if (!is.numeric(stages) || length(stages) < 2L) {
    stop("'stages' must be a named numeric vector with at least 2 entries.",
         call. = FALSE)
  }
  .assert_positive_integer(n_mc, "n_mc")
  if (n_mc < 2L) stop("'n_mc' must be at least 2.", call. = FALSE)
  if (any(!is.finite(stages)) || any(stages < 0) || sum(stages) == 0) {
    stop("Stage variances must be finite, non-negative and have a positive sum.", call. = FALSE)
  }
  stage_names <- names(stages)
  if (is.null(stage_names)) stage_names <- paste0("V", seq_along(stages))
  names(stages) <- stage_names

  K <- length(stages)
  sig <- sqrt(stages)

  # Build covariance matrix
  if (!is.null(covariances)) {
    if (!is.matrix(covariances) || nrow(covariances) != K || ncol(covariances) != K) {
      stop("'covariances' must be a K x K matrix.", call. = FALSE)
    }
    Sigma <- covariances
    diag(Sigma) <- stages
  } else {
    Sigma <- diag(stages, nrow = K)
    if (!is.null(correlations)) {
      for (nm in names(correlations)) {
        kk <- strsplit(nm, "__")[[1L]]
        if (length(kk) != 2L || !all(kk %in% stage_names)) {
          warning(sprintf("Ignoring malformed correlation entry: %s", nm),
                   call. = FALSE)
          next
        }
        i <- match(kk[1L], stage_names); j <- match(kk[2L], stage_names)
        rho <- correlations[[nm]]
        if (!is.numeric(rho) || length(rho) != 1L || abs(rho) > 1) {
          warning(sprintf("Invalid rho in '%s'; ignoring.", nm), call. = FALSE)
          next
        }
        Sigma[i, j] <- Sigma[j, i] <- rho * sig[i] * sig[j]
      }
    }
  }

  if (!isSymmetric(Sigma) || any(!is.finite(Sigma)) ||
      min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) < -1e-10) {
    stop("The covariance matrix must be finite, symmetric and positive semidefinite.", call. = FALSE)
  }

  # Total variance: trace(Sigma) + 2 * sum off-diagonal / 2 = sum of all entries
  total_var_additive <- sum(diag(Sigma))
  total_var_correlated <- sum(Sigma)

  # Shapley allocation for the subset-sum variance game:
  # phi_k = sigma_k^2 + sum_{j != k} Cov(k, j)
  # Because each pair (k,j) contributes Cov(k,j) twice (once to phi_k, once to phi_j),
  # sum(phi_k) = sum(diag(Sigma)) + 2 * sum_{k<j} Cov(k,j) = total_var_correlated.
  shapley <- rowSums(Sigma)
  names(shapley) <- stage_names

  # Percentages
  pct_additive <- stages / total_var_additive
  pct_shapley  <- shapley / total_var_correlated

  # Monte Carlo propagation (for nonlinear transmission)
  .set_seed(seed)
  samples <- MASS::mvrnorm(n_mc, mu = rep(0, K), Sigma = Sigma)
  total_mc <- rowSums(samples)
  total_var_mc <- stats::var(total_mc)
  # MCSE
  mcse_total_var <- stats::sd(total_mc^2) / sqrt(n_mc)

  alpha <- 1 - conf_level
  chi_lo <- stats::qchisq(1 - alpha / 2, df = n_mc - 1L)
  chi_hi <- stats::qchisq(alpha / 2, df = n_mc - 1L)
  total_variance_ci <- c((n_mc - 1L) * total_var_mc / chi_lo,
                         (n_mc - 1L) * total_var_mc / chi_hi)
  total_prediction_interval <- stats::quantile(total_mc,
                                  c(alpha / 2, 1 - alpha / 2))

  # Gaussian sample-covariance precision for row sums of Sigma.
  shapley_mcse <- sqrt(pmax(diag(Sigma) * total_var_correlated + shapley^2, 0) /
                         (n_mc - 1L))
  zcrit <- stats::qnorm(1 - alpha / 2)
  shapley_ci_low <- shapley - zcrit * shapley_mcse
  shapley_ci_high <- shapley + zcrit * shapley_mcse
  shapley_mcse <- sqrt(pmax(stages * total_var_correlated + shapley^2, 0) /
                         (n_mc - 1L))
  z <- stats::qnorm(1 - alpha / 2)
  shapley_ci_low <- shapley - z * shapley_mcse
  shapley_ci_high <- shapley + z * shapley_mcse
  out <- .crossstatr_result(
    list(stages = stages,
         covariance = Sigma,
         total_variance_additive = total_var_additive,
         total_variance_correlated = total_var_correlated,
         total_variance_mc = total_var_mc,
         total_variance_mcse = mcse_total_var,
         total_variance_ci = total_variance_ci,
         total_prediction_interval = total_prediction_interval,
         total_ci = total_prediction_interval,
         total_ci = total_prediction_interval,
         interval_scope = "Monte Carlo precision conditional on supplied covariance",
         shapley_mcse = shapley_mcse,
         shapley_attribution = shapley,
         shapley_ci = cbind(shapley_ci_low, shapley_ci_high),
         percent_additive = pct_additive,
         percent_shapley = pct_shapley,
         model_selection_uncertainty = model_selection_uncertainty,
         n_mc = n_mc,
         conf_level = conf_level,
         method = "Uncertainty-origin attribution framework",
         call = call),
    "crossstatr_uncertainty_origin")
  out
}

#' @export
print.crossstatr_uncertainty_origin <- function(x, ...) {
  cat(sprintf("<crossStatR::uncertainty_origin> %s\n", x$method))
  cat(sprintf("Total variance: %.4g (additive), %.4g (correlated)\n",
              x$total_variance_additive, x$total_variance_correlated))
  cat(sprintf("Total variance (MC, n=%d): %.4g (MCSE ~ %.4g)\n",
              x$n_mc, x$total_variance_mc, x$total_variance_mcse))
  cat("\nAttribution:\n")
  df <- data.frame(
    stage = names(x$stages),
    variance = as.numeric(x$stages),
    pct_additive = as.numeric(x$percent_additive),
    pct_shapley = as.numeric(x$percent_shapley)
  )
  print(df, row.names = FALSE)
  if (!is.null(x$model_selection_uncertainty)) {
    cat(sprintf("\nModel-selection uncertainty (reported separately): %.4g\n",
                x$model_selection_uncertainty))
  }
  invisible(x)
}
