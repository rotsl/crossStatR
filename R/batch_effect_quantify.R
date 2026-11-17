#' Quantify batch effects
#'
#' Decomposes the variance in an experimental outcome into treatment,
#' batch, and residual components, and reports the proportion of variance
#' attributable to batch. Two model parameterisations are offered: a
#' linear mixed model with a random intercept on batch (`type = "lmm"`,
#' recommended when batches are a random sample from a larger population
#' of batches), and a fixed-effect model with batch as a categorical
#' covariate (`type = "fixed"`, recommended when the observed batches
#' are the only ones of interest).
#'
#' Confounding between treatment and batch is checked by a chi-square
#' test of association. If treatment is perfectly confounded with batch
#' (each treatment level occurring in only one batch, and vice versa),
#' the function stops with an error because the treatment effect is not
#' estimable.
#'
#' @param formula a formula of the form `outcome ~ treatment`. Only one
#'   predictor is allowed.
#' @param data data frame containing the variables in `formula` and the
#'   batch column.
#' @param batch either a character string naming a column in `data`, or a
#'   vector (factor or coercible) of batch identifiers of length
#'   `nrow(data)`.
#' @param type `"lmm"` (default) for a random-intercept LMM via
#'   [lme4::lmer()], or `"fixed"` for a fixed-effect linear model via
#'   [stats::lm()].
#' @param conf_level confidence level; stored on the result for downstream
#'   reporting.
#'
#' @return A `crossstatr_batch` object with elements:
#'   \item{batch_variance}{estimated batch (random-intercept) variance.
#'     For `type = "fixed"` this is the one-way random-effects variance
#'     estimate derived from the batch sum of squares.}
#'   \item{batch_pct_variance}{batch variance as a percentage of total
#'     (batch + residual) variance.}
#'   \item{icc}{intraclass correlation coefficient
#'     \eqn{\sigma^2_{\mathrm{batch}} / (\sigma^2_{\mathrm{batch}} + \sigma^2_{\mathrm{resid}})}.}
#'   \item{confound_p}{chi-square p-value for the treatment-by-batch
#'     association test. A small p indicates partial confounding.}
#'   \item{confound_test}{list with `statistic`, `p_value`, `df` from a
#'     chi-square test of association between treatment and batch.}
#'   \item{model}{model label.}
#'   \item{treatment_effect_change}{change in the treatment coefficient
#'     when batch is added to the model.}
#'   \item{fit, base_fit}{the fitted full and treatment-only models.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   For `type = "lmm"` The package fits the random-intercept model
#'   \deqn{y_{ij} = \mu + \alpha\,x_{ij} + b_j + \epsilon_{ij}, \qquad
#'         b_j \sim N(0, \sigma^2_b), \quad \epsilon_{ij} \sim N(0, \sigma^2),}
#'   where `j` indexes batch. The intraclass correlation coefficient is
#'   \deqn{\mathrm{ICC} = \frac{\sigma^2_b}{\sigma^2_b + \sigma^2}.}
#'
#'   For `type = "fixed"` The package fits the linear model
#'   \deqn{y_{ij} = \mu + \alpha\,x_{ij} + \beta_j + \epsilon_{ij},}
#'   compare it to the treatment-only model via [stats::anova()], and
#'   derive an approximate batch variance from the batch sum of squares
#'   using the one-way random-effects decomposition (Searle, Casella and
#'   McCulloch, 1992, sec. 3.7):
#'   \deqn{\hat{\sigma}^2_b \approx \frac{\mathrm{SS}_{\mathrm{batch}} /
#'         \mathrm{df}_{\mathrm{batch}} - \hat{\sigma}^2}{n_0},}
#'   where \eqn{n_0} is the effective number of observations per batch
#'   \eqn{(N - \sum_j n_j^2 / N) / (J - 1)}.
#'
#'   Treatment-effect change: \eqn{\Delta\alpha = \alpha_{\mathrm{full}}
#'   - \alpha_{\mathrm{base}}}. A large absolute value signals that batch
#'   was confounded with treatment.
#'
#'   Confounding check: a chi-square test of independence between
#'   treatment and batch. Perfect confounding raises an error.
#'
#' @section Limitations:
#'
#'   - The LMM assumes batches are a random sample from a larger
#'     population of batches. If only a small number of batches is
#'     observed, the batch-variance estimate is unstable; consider the
#'     fixed-effect model in that case.
#'   - The fixed-effect model does not generalise to unobserved batches.
#'     The reported batch variance is a method-of-moments estimate and
#'     should be interpreted as a description of the observed batches
#'     only.
#'   - ICC for the fixed-effect path is approximate and can be negative
#'     when between-batch variation is small; the package truncates at zero.
#'   - The function does not check whether the treatment effect itself
#'     is well-estimated. A small `treatment_effect_change` is
#'     consistent with no confounding, but it does not rule out subtle
#'     confounding that the chi-square test misses.
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
#' # 3 treatments x 4 batches x 5 obs per (treatment, batch) cell,
#' # with batch-level shifts of 0, 1, -1, 0.5 respectively
#' df <- data.frame(
#'   y = rnorm(60, mean = rep(c(0, 1, 2), each = 20) +
#'                       rep(rep(c(0, 1, -1, 0.5), each = 5), times = 3)),
#'   treat = factor(rep(c("A","B","C"), each = 20)),
#'   batch = factor(rep(1:4, each = 5, times = 3))
#' )
#' r_lmm <- batch_effect_quantify(y ~ treat, data = df, batch = "batch",
#'                                type = "lmm")
#' print(r_lmm)
#'
#' r_fixed <- batch_effect_quantify(y ~ treat, data = df, batch = "batch",
#'                                  type = "fixed")
#' print(r_fixed)
#'
#' @export
batch_effect_quantify <- function(formula, data, batch,
                                   type = c("lmm", "fixed"),
                                   conf_level = 0.95) {
  call <- match.call()
  type <- match.arg(type)
  .assert_data_frame(data, "data")
  .assert_probability(conf_level, "conf_level")

  formula <- stats::formula(formula)
  terms_obj <- stats::terms(formula)
  response_var <- as.character(attr(terms_obj, "variables"))[2L]
  pred_vars <- attr(terms_obj, "term.labels")
  if (length(pred_vars) != 1L) {
    stop("Formula must have exactly one predictor (treatment).",
         call. = FALSE)
  }
  treat_var <- pred_vars[1L]
  if (!(response_var %in% names(data))) {
    stop(sprintf("Response '%s' not found in data.", response_var),
         call. = FALSE)
  }
  if (!(treat_var %in% names(data))) {
    stop(sprintf("Treatment '%s' not found in data.", treat_var),
         call. = FALSE)
  }

  # Batch handling: column name or vector.
  if (is.character(batch) && length(batch) == 1L) {
    batch_col <- batch
    if (!(batch_col %in% names(data))) {
      stop(sprintf("Batch column '%s' not in data.", batch_col),
           call. = FALSE)
    }
    data_full <- data
    batch_fac <- factor(data[[batch_col]])
  } else {
    if (length(batch) != nrow(data)) {
      stop("'batch' must be a column name or a vector of length nrow(data).",
           call. = FALSE)
    }
    batch_col <- ".batch"
    data_full <- data
    data_full$.batch <- factor(batch)
    batch_fac <- data_full$.batch
  }

  treat_fac <- factor(data_full[[treat_var]])

  # Confounding check: chi-square test of treatment vs batch.
  tab <- table(treat_fac, batch_fac)
  if (all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)) {
    stop("Treatment is perfectly confounded with batch; the treatment\n",
         "effect is not estimable. Redesign the experiment or restrict\n",
         "the analysis to within-batch comparisons.",
         call. = FALSE)
  }
  correct_2x2 <- all(dim(tab) == c(2L, 2L))
  chi_test <- suppressWarnings(
    stats::chisq.test(tab, correct = correct_2x2)
  )
  confound_p <- as.numeric(chi_test$p.value)
  chi_out <- list(statistic = unname(chi_test$statistic),
                  p_value = confound_p,
                  df = unname(chi_test$parameter),
                  table = tab)

  # Base model: treatment only.
  base_form <- stats::as.formula(paste(response_var, "~", treat_var))
  base_fit <- stats::lm(base_form, data = data_full)

  # Helper to extract treatment coefficient (skipping intercept).
  get_treat_coef <- function(fit, treat_var) {
    if (inherits(fit, "merMod")) {
      cn <- names(lme4::fixef(fit))
      idx <- which(!cn %in% "(Intercept)" & grepl(treat_var, cn))
      if (length(idx) == 0L) idx <- 2L
      return(as.numeric(lme4::fixef(fit)[idx[1L]]))
    }
    cn <- names(stats::coef(fit))
    idx <- which(!cn %in% "(Intercept)" & grepl(treat_var, cn))
    if (length(idx) == 0L) idx <- 2L
    as.numeric(stats::coef(fit)[idx[1L]])
  }

  if (type == "lmm") {
    full_form <- stats::as.formula(
      paste(response_var, "~", treat_var, "+ (1|", batch_col, ")")
    )
    lmer_fit <- tryCatch(
      lme4::lmer(full_form, data = data_full, REML = TRUE),
      error = function(e) {
        stop("lmer fit failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    vc <- lme4::VarCorr(lmer_fit)
    vc_df <- as.data.frame(vc)
    batch_row <- vc_df[vc_df$grp == batch_col, , drop = FALSE]
    batch_var_est <- if (nrow(batch_row) > 0L) batch_row$vcov[1L] else NA_real_
    resid_var_est <- attr(vc, "sc")^2
    icc <- if (is.finite(batch_var_est) && is.finite(resid_var_est) &&
              (batch_var_est + resid_var_est) > 0) {
      batch_var_est / (batch_var_est + resid_var_est)
    } else NA_real_
    batch_pct <- if (is.finite(batch_var_est) &&
                     (batch_var_est + resid_var_est) > 0) {
      100 * batch_var_est / (batch_var_est + resid_var_est)
    } else NA_real_

    treat_coef_base <- get_treat_coef(base_fit, treat_var)
    treat_coef_full <- get_treat_coef(lmer_fit, treat_var)
    treat_change <- treat_coef_full - treat_coef_base

    tidy_df <- data.frame(
      component = c("batch", "residual"),
      variance = c(batch_var_est, resid_var_est),
      pct = c(batch_pct, 100 - batch_pct),
      stringsAsFactors = FALSE
    )

    out <- list(batch_variance = batch_var_est,
                batch_pct_variance = batch_pct,
                icc = icc,
                confound_p = confound_p,
                confound_test = chi_out,
                model = "lmm (random intercept)",
                treatment_effect_change = unname(treat_change),
                fit = lmer_fit,
                base_fit = base_fit,
                conf_level = conf_level,
                tidy = tidy_df,
                call = call)
  } else {
    # Fixed-effects path: lm(outcome ~ treatment + factor(batch)).
    full_form <- stats::as.formula(
      paste(response_var, "~", treat_var, "+", batch_col)
    )
    full_fit <- stats::lm(full_form, data = data_full)
    an <- stats::anova(base_fit, full_fit)

    df_batch <- length(levels(batch_fac)) - 1L
    df_resid_full <- full_fit$df.residual
    rss_full <- sum(stats::resid(full_fit)^2)
    rss_base <- sum(stats::resid(base_fit)^2)
    batch_ss <- rss_base - rss_full
    total_ss <- rss_base

    # Approximate one-way random effects variance from batch SS.
    sigma2_resid <- rss_full / df_resid_full
    n_per_batch <- as.numeric(table(batch_fac))
    J <- length(n_per_batch)
    N <- sum(n_per_batch)
    n0 <- (N - sum(n_per_batch^2) / N) / max(1L, J - 1L)
    sigma2_batch <- max(0, (batch_ss / df_batch - sigma2_resid) / n0)
    total_var <- sigma2_batch + sigma2_resid
    icc <- if (total_var > 0) sigma2_batch / total_var else NA_real_

    batch_pct <- if (total_var > 0) 100 * sigma2_batch / total_var else NA_real_
    treat_coef_base <- get_treat_coef(base_fit, treat_var)
    treat_coef_full <- get_treat_coef(full_fit, treat_var)
    treat_change <- treat_coef_full - treat_coef_base

    tidy_df <- data.frame(
      component = c("batch", "residual"),
      variance = c(sigma2_batch, sigma2_resid),
      pct = c(if (total_var > 0) 100 * sigma2_batch / total_var else NA_real_,
              if (total_var > 0) 100 * sigma2_resid / total_var else NA_real_),
      stringsAsFactors = FALSE
    )

    out <- list(batch_variance = sigma2_batch,
                batch_pct_variance = batch_pct,
                icc = icc,
                confound_p = confound_p,
                confound_test = chi_out,
                model = "fixed effects (lm)",
                treatment_effect_change = unname(treat_change),
                fit = full_fit,
                base_fit = base_fit,
                anova = an,
                batch_ss = batch_ss,
                conf_level = conf_level,
                tidy = tidy_df,
                call = call)
  }

  .crossstatr_result(out, "crossstatr_batch")
}

#' @export
print.crossstatr_batch <- function(x, ...) {
  cat(sprintf("<crossStatR::batch_effect_quantify> %s\n", x$model))
  cat(sprintf("batch variance = %.4g  (%.1f%% of total)\n",
              x$batch_variance, x$batch_pct_variance))
  cat(sprintf("ICC = %.4g\n", x$icc))
  ct <- x$confound_test
  cat(sprintf("confounding chi-square: statistic = %.4g, df = %g, p = %g\n",
              ct$statistic, ct$df, ct$p_value))
  cat(sprintf("treatment coefficient change after adding batch: %.4g\n",
              x$treatment_effect_change))
  invisible(x)
}
