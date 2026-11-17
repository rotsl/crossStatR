#' ANOVA and related analyses
#'
#' Wraps and extends `stats::aov()` / `stats::lm()` to provide
#' effect-size reporting, multiple-comparison adjustments, and
#' Welch-type alternatives for heteroscedastic designs.
#'
#' @param formula a model formula.
#' @param data data frame.
#' @param welch logical; use Welch's ANOVA (`oneway.test`) for one-way
#'   designs with unequal variances? Default FALSE.
#' @param posthoc character; `"tukey"` (default), `"bonferroni"`, `"holm"`,
#'   or `"none"`.
#' @param conf_level confidence level for CIs.
#' @param ... additional arguments passed to `aov` / `lm`.
#'
#' @return A `crossstatr_anova` object with `aov`, `anova_table`,
#'   `effect_sizes`, `posthoc` (if requested), and `diagnostics`.
#'
#' @section Methodology:
#'
#'   One-way ANOVA partitions the total sum of squares:
#'   \deqn{SS_{\mathrm{total}} = SS_{\mathrm{between}} + SS_{\mathrm{within}},}
#'   with F-statistic
#'   \deqn{F = \frac{SS_{\mathrm{between}} / (k - 1)}
#'                {SS_{\mathrm{within}} / (N - k)}.}
#'
#'   Welch's ANOVA (Welch, 1951) replaces the within-group pooled
#'   variance with group-specific variances, and uses an
#'   Satterthwaite-adjusted denominator df.
#'
#'   Effect sizes:
#'   - Eta-squared \eqn{\eta^2 = SS_{\mathrm{between}} / SS_{\mathrm{total}}}
#'   - Partial eta-squared \eqn{\eta_p^2 = SS_{\mathrm{effect}} /
#'     (SS_{\mathrm{effect}} + SS_{\mathrm{error}})}
#'   - Omega-squared (see Carroll and Nordholm, 1975):
#'     \eqn{\omega^2 = \frac{SS_{\mathrm{between}} - (k-1) MS_{\mathrm{within}}}
#'     {SS_{\mathrm{total}} + MS_{\mathrm{within}}}}
#'
#'   Tukey HSD uses the studentized range distribution.
#'
#' @section Limitations:
#'
#'   - Welch effect sizes are descriptive sums-of-squares summaries, not
#'     heteroscedasticity-corrected inferential estimates. Factorial effects
#'     use sequential (Type I) sums of squares and depend on term order.
#'     Welch inference is restricted to one-way models. Tukey post-hoc tests
#'     are unavailable on that path; a warning accompanies such requests.
#'
#' @references
#'   - Welch, B. L. (1951). On the comparison of several mean values: An
#'     alternative approach. Biometrika, 38, 330-336.
#'     \doi{10.1093/biomet/38.3-4.330}.
#'
#'   - Kelley, T. L. (1935). An Unbiased Correlation Ratio Measure.
#'     Proceedings of the National Academy of Sciences, 21, 554-559.
#'     \doi{10.1073/pnas.21.9.554}.
#'
#'   - Carroll, R. M. and Nordholm, L. A. (1975). Sampling characteristics
#'     of Kelley's epsilon-squared and Hays' omega-squared. Educational and
#'     Psychological Measurement, 35, 541-554.
#'     \doi{10.1177/001316447503500304}.
#'
#'   - Hsu, J. C. (1996). Multiple Comparisons: Theory and Methods. Chapman
#'     & Hall.
#'     \doi{10.1201/b15074}.
#'
#' @examples
#' df <- data.frame(y = rnorm(60, mean = rep(c(5, 6, 8), each = 20)),
#'                  g = factor(rep(c("A","B","C"), each = 20)))
#' sci_anova(y ~ g, data = df)
#'
#' @export
sci_anova <- function(formula, data, welch = FALSE,
                      posthoc = c("tukey", "bonferroni", "holm", "none"),
                      conf_level = 0.95, ...) {
  call <- match.call()
  .assert_probability(conf_level, "conf_level")
  .assert_flag(welch, "welch")
  posthoc <- match.arg(posthoc)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (nrow(mf) == 0L) stop("Empty model frame.", call. = FALSE)

  # Detect one-way vs factorial
  terms_obj <- stats::terms(formula, data = data)
  pred_vars <- attr(terms_obj, "term.labels")

  if (welch && length(pred_vars) != 1L) {
    stop("Welch ANOVA requires a one-way formula.", call. = FALSE)
  }
  if (welch && length(pred_vars) == 1L) {
    wt <- stats::oneway.test(formula, data = data)
    aov_obj <- NULL
    anova_table <- data.frame(
      `Df` = wt$parameter[1L],
      denominator_df = wt$parameter[2L],
      `Sum Sq` = NA_real_,
      `Mean Sq` = NA_real_,
      `F value` = as.numeric(wt$statistic),
      `Pr(>F)` = wt$p.value,
      check.names = FALSE)
  } else {
    aov_obj <- stats::aov(formula, data = data, ...)
    anova_table <- summary(aov_obj)[[1L]]
    anova_table <- as.data.frame(anova_table)
    colnames(anova_table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  }

  # Effect sizes (one-way only; for factorial use partial eta^2 per term)
  ss_total <- sum(anova_table$`Sum Sq`, na.rm = TRUE)
  ms_within <- anova_table$`Mean Sq`[nrow(anova_table)]
  ss_within <- anova_table$`Sum Sq`[nrow(anova_table)]
  if (welch && length(pred_vars) == 1L) {
    grp <- mf[[pred_vars[1L]]]; vals <- mf[[1L]]
    gmeans <- tapply(vals, grp, mean)
    gvars  <- tapply(vals, grp, var)
    ns     <- tapply(vals, grp, length)
    grand  <- mean(vals)
    ss_between <- sum(ns * (gmeans - grand)^2)
    ss_total_w  <- sum((vals - grand)^2)
    ms_within_w <- sum((ns - 1) * gvars) / (sum(ns) - length(ns))
    k_w <- length(gmeans)
    eta2  <- ss_between / ss_total_w
    eta2_p <- ss_between / (ss_between + sum((ns - 1) * gvars))
    omega2 <- (ss_between - (k_w - 1) * ms_within_w) / (ss_total_w + ms_within_w)
    effect_sizes <- data.frame(
      effect = c("eta_squared", "partial_eta_squared", "omega_squared"),
      estimate = c(eta2, eta2_p, omega2),
      interpretation = c(
        .interpret_eta(eta2),
        .interpret_eta(eta2_p),
        .interpret_eta(omega2)
      ),
      stringsAsFactors = FALSE)
  } else if (length(pred_vars) == 1L) {
    ss_effect <- anova_table$`Sum Sq`[1L]
    k <- anova_table$Df[1L] + 1L
    eta2 <- ss_effect / ss_total
    eta2_p <- ss_effect / (ss_effect + ss_within)
    omega2 <- (ss_effect - (k - 1) * ms_within) / (ss_total + ms_within)
    effect_sizes <- data.frame(
      effect = c("eta_squared", "partial_eta_squared", "omega_squared"),
      estimate = c(eta2, eta2_p, omega2),
      interpretation = c(
        .interpret_eta(eta2),
        .interpret_eta(eta2_p),
        .interpret_eta(omega2)
      ),
      stringsAsFactors = FALSE)
  } else {
    terms_df <- anova_table[seq_len(nrow(anova_table) - 1L), , drop = FALSE]
    effect_sizes <- do.call(rbind, lapply(seq_len(nrow(terms_df)), function(i) {
      ss_e <- terms_df$`Sum Sq`[i]
      df_e <- terms_df$Df[i]
      eta2_i  <- ss_e / ss_total
      eta2_p_i <- ss_e / (ss_e + ss_within)
      omega2_i <- (ss_e - df_e * ms_within) / (ss_total + ms_within)
      data.frame(
        term = rownames(terms_df)[i],
        effect = c("eta_squared", "partial_eta_squared", "omega_squared"),
        estimate = c(eta2_i, eta2_p_i, omega2_i),
        interpretation = c(
          .interpret_eta(eta2_i),
          .interpret_eta(eta2_p_i),
          .interpret_eta(omega2_i)
        ),
        stringsAsFactors = FALSE)
    }))
  }

  # Post-hoc (only for one-way, only if aov was fit, i.e. not Welch path)
  posthoc_res <- NULL
  if (welch && posthoc != "none" && length(pred_vars) == 1L) {
    warning("Post-hoc 'Tukey' is not implemented for Welch's ANOVA. ",
            "Consider Games-Howell or pairwise Welch t-tests.",
            call. = FALSE)
  } else if (posthoc != "none" && length(pred_vars) == 1L && !is.null(aov_obj)) {
    if (posthoc == "tukey") {
      tk <- stats::TukeyHSD(aov_obj, conf.level = conf_level)
      posthoc_res <- as.data.frame(tk[[1L]])
    } else {
      # Pairwise t-tests with adjustment
      grp <- mf[[pred_vars[1L]]]
      vals <- mf[[1L]]
      pw <- stats::pairwise.t.test(vals, grp, p.adjust.method = posthoc)
      posthoc_res <- list(p_values = pw$p.value, method = posthoc)
    }
  }

  # Diagnostics
  diag <- list(
    shapiro_p = if (!is.null(aov_obj)) {
      tryCatch(stats::shapiro.test(stats::resid(aov_obj))$p.value,
               error = function(e) NA_real_)
    } else NA_real_,
    levene_p = if (length(pred_vars) == 1L) {
      .levene_test(mf[[1L]], mf[[pred_vars[1L]]])
    } else NA_real_
  )

  .crossstatr_result(
    list(aov = aov_obj, anova_table = anova_table,
         effect_sizes = effect_sizes,
         posthoc = posthoc_res,
         diagnostics = diag,
         welch = welch,
         method = if (welch) "Welch's ANOVA" else "ANOVA",
         call = call),
    "crossstatr_anova")
}

.interpret_eta <- function(e) {
  if (is.na(e)) return(NA_character_)
  if (e < 0.01) "negligible"
  else if (e < 0.06) "small"
  else if (e < 0.14) "medium"
  else "large"
}

.levene_test <- function(y, group) {
  med <- tapply(y, group, stats::median, na.rm = TRUE)
  z <- abs(y - med[as.integer(factor(group))])
  tryCatch({
    an <- stats::oneway.test(z ~ factor(group))
    an$p.value
  }, error = function(e) NA_real_)
}

# null-coalesce helper
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.crossstatr_anova <- function(x, ...) {
  cat(sprintf("<crossStatR::anova> %s\n", x$method))
  print(x$anova_table)
  cat("\nEffect sizes:\n"); print(x$effect_sizes, row.names = FALSE)
  if (!is.null(x$diagnostics$levene_p) && !is.na(x$diagnostics$levene_p))
    cat(sprintf("\nLevene's test p = %.4f\n", x$diagnostics$levene_p))
  invisible(x)
}
