#' Variance component decomposition for nested and crossed replicates
#'
#' Decomposes the total variance of a measurement into contributions from
#' each level of replication: biological, technical, and residual. The
#' decomposition is done by a linear mixed model fit by `lme4::lmer` with
#' REML. The formula should contain random-effects terms such as
#' `(1|bio) + (1|bio:tech)` to encode the nesting or crossing structure.
#'
#' Variance components are read from `lme4::VarCorr`, and Wald-type
#' confidence intervals are constructed by treating the standard deviations
#' (square root of the diagonal of each `VarCorr` block) as approximately
#' normal. This is the standard large-sample approximation; for designs
#' with few levels per random effect, profile-likelihood CIs from
#' `confint(model, method = "profile")` are more accurate and can be
#' obtained directly from the stored `model` element.
#'
#' @param formula a model formula containing random-effects terms, e.g.
#'   `y ~ treatment + (1|bio) + (1|bio:tech)`.
#' @param data data frame containing the variables in `formula`.
#' @param type character; currently only `"lmm"` is supported.
#' @param conf_level confidence level for the variance component CIs.
#'
#' @return A `crossstatr_replicate` object with elements:
#'   \item{variance_components}{data frame with columns `component`,
#'     `variance`, `sd`, `proportion`, `ci_low`, `ci_high`. The `Residual`
#'     row is always last.}
#'   \item{model}{the underlying `lmerMod` fit.}
#'   \item{method}{backend used (`"LMM (REML)"`).}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   For the linear mixed model
#'   \deqn{y_{ij} = X_{ij}\beta + \sum_{l} b^{(l)}_{g_l(i,j)} +
#'         \epsilon_{ij},}
#'   where \eqn{b^{(l)} \sim N(0, \sigma_l^2)} are random effects with
#'   grouping factor \eqn{g_l} and \eqn{\epsilon \sim N(0, \sigma_r^2)} is
#'   the residual, the variance components are \eqn{\sigma_l^2} for each
#'   random effect and \eqn{\sigma_r^2} for the residual. The total
#'   variance is \eqn{\sigma_{\mathrm{total}}^2 = \sum_l \sigma_l^2 +
#'   \sigma_r^2}, and the proportion attributable to level \eqn{l} is
#'   \eqn{p_l = \sigma_l^2 / \sigma_{\mathrm{total}}^2}. The intra-class
#'   correlation for a single level is \eqn{\mathrm{ICC}_l =
#'   \sigma_l^2 / \sigma_{\mathrm{total}}^2}.
#'
#'   The Wald CIs treat the standard deviation estimates
#'   \eqn{\hat\sigma_l = \sqrt{\mathrm{diag}(\mathrm{VarCorr}_l)}} as
#'   approximately normal. The asymptotic variance of \eqn{\hat\sigma_l} is
#'   approximated via the delta method on the variance estimate:
#'   \deqn{\mathrm{SE}(\hat\sigma_l) \approx \frac{\hat\sigma_l}{\sqrt{2\,
#'         \nu_l}},}
#'   where \eqn{\nu_l} is the degrees of freedom associated with component
#'   \eqn{l} (number of groups minus one for random effects, residual df
#'   for the residual). The CI on the SD is
#'   \eqn{\hat\sigma_l \pm z_{1-\alpha/2}\, \mathrm{SE}(\hat\sigma_l)},
#'   clipped to be non-negative, and the CI on the variance is the square
#'   of the SD CI bounds.
#'
#' @section Limitations:
#'
#'   - The LMM assumes normally distributed random effects and residuals.
#'     For count or proportion outcomes, use `lme4::glmer` with an
#'     appropriate family directly.
#'   - Wald CIs on variance components are anti-conservative when the true
#'     variance is near zero (the boundary of the parameter space). For
#'     designs where this matters, use
#'     `confint(model, method = "profile")` directly on the stored `model`.
#'   - The degrees-of-freedom approximation \eqn{\nu_l = n_l - 1} is rough
#'     for unbalanced or crossed designs. Treat the CIs as approximate.
#'   - The proportion of variance attributable to each level is meaningful
#'     only when the model is correctly specified. Including or excluding a
#'     fixed effect changes the variance decomposition.
#'
#' @references
#'   - Bates, D., Mächler, M., Bolker, B. and Walker, S. (2015). Fitting
#'     Linear Mixed-Effects Models Using lme4. Journal of Statistical
#'     Software, 67, 1-48.
#'     \doi{10.18637/jss.v067.i01}.
#'
#'   - Searle, S. R., Casella, G. and McCulloch, C. E. (1992). Variance
#'     Components. Wiley.
#'     \doi{10.1002/9780470316856}.
#'
#'   - Pinheiro, J. C. and Bates, D. M. (2000). Mixed-Effects Models in S
#'     and S-PLUS. Springer New York.
#'     \doi{10.1007/978-1-4419-0318-1}.
#'
#' @examples
#' set.seed(1)
#' # 4 bio samples, 3 tech replicates within each bio, 5 measurements per tech
#' df <- data.frame(
#'   y = rnorm(60, mean = rep(rnorm(12, sd = 1.5), each = 5), sd = 0.5),
#'   bio = factor(rep(1:4, each = 15)),
#'   tech = factor(rep(rep(1:3, each = 5), 4)))
#' fit <- replicate_decompose(y ~ 1 + (1|bio) + (1|bio:tech), data = df)
#' print(fit)
#'
#' @export
replicate_decompose <- function(formula, data, type = "lmm",
                              conf_level = 0.95) {
  call <- match.call()
  if (!identical(type, "lmm")) {
    stop("Currently only type = 'lmm' is supported.", call. = FALSE)
  }
  .assert_data_frame(data, "data")
  .assert_probability(conf_level, "conf_level")

  formula_chr <- deparse(formula)
  if (!grepl("\\|", formula_chr)) {
    stop("The formula must contain random-effects terms such as ",
         "(1|bio).", call. = FALSE)
  }

  fit <- lme4::lmer(formula, data = data, REML = TRUE)
  vc <- .extract_varcor_lmm(fit, conf_level = conf_level)

  .crossstatr_result(
    list(
      variance_components = vc,
      model = fit,
      method = "LMM (REML)",
      call = call
    ),
    "crossstatr_replicate"
  )
}

# Extract variance components from a lmerMod, with Wald CIs on the SDs.
# df for each random effect = n_groups - 1; df for residual = N - p.
# @noRd
.extract_varcor_lmm <- function(fit, conf_level) {
  vc <- lme4::VarCorr(fit)
  comp_names <- names(vc)
  sds <- vapply(vc, function(x) attr(x, "stddev")[[1L]], numeric(1L))
  vars <- sds^2
  resid_sd <- attr(vc, "sc")
  resid_var <- resid_sd^2
  comp_names <- c(comp_names, "Residual")
  sds <- c(sds, resid_sd)
  vars <- c(vars, resid_var)
  total_var <- sum(vars)
  prop <- vars / total_var

  # Degrees of freedom per component.
  ngrps <- vapply(lme4::getME(fit, "flist"), nlevels, integer(1L))
  n_obs <- as.integer(lme4::getME(fit, "n"))
  p_fixef <- length(lme4::fixef(fit))
  df_vec <- c(as.numeric(ngrps) - 1, n_obs - p_fixef - sum(as.numeric(ngrps)))
  df_vec <- pmax(df_vec, 1)

  # Wald CI on SD: sigma_hat +/- z * sigma_hat / sqrt(2 * df).
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  sd_se <- sds / sqrt(2 * df_vec)
  sd_lo <- pmax(0, sds - z * sd_se)
  sd_hi <- sds + z * sd_se
  var_ci_lo <- sd_lo^2
  var_ci_hi <- sd_hi^2

  data.frame(
    component = comp_names,
    variance = vars,
    sd = sds,
    proportion = prop,
    ci_low = var_ci_lo,
    ci_high = var_ci_hi,
    stringsAsFactors = FALSE
  )
}

#' @export
print.crossstatr_replicate <- function(x, ...) {
  cat(sprintf("<crossStatR::replicate_decompose> method = %s\n", x$method))
  cat("\nVariance components:\n")
  print(x$variance_components, row.names = FALSE)
  invisible(x)
}
