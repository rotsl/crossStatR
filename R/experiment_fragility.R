#' Quantitative experimental fragility assessment
#'
#' Quantifies how sensitive a reported result is to defensible analysis
#' decisions: which replicates to include, how to handle outliers, which
#' baseline correction to use, which model to fit, and how to propagate
#' measurement uncertainty. This is a sensitivity analysis, not a
#' confirmatory test; it does not change the original result.
#'
#' @param analysis_fn function that takes a list of options and returns
#'   a numeric scalar (or a small named vector) summary statistic of
#'   interest. Required.
#' @param data data frame passed to `analysis_fn` as the `data` argument.
#' @param replicates optional character; name of the column identifying
#'   replicates. If supplied, leave-one-replicate-out is computed.
#' @param outlier_methods character vector; subset of
#'   `c("none", "iqr", "mad")`. Defaults to all three.
#' @param baseline_methods character vector; subset of
#'   `c("none", "asls", "lowess")`. Only used if `analysis_fn` accepts
#'   a `baseline` argument.
#' @param model_alternatives optional list of named alternative model
#'   specifications. Only used if `analysis_fn` accepts a `model`
#'   argument.
#' @param measurement_uncertainties optional numeric vector of SDs, used
#'   to compute a fragility index: the smallest noise level at which the
#'   result loses significance (i.e. crosses a threshold).
#' @param decision_threshold numeric scalar; the value that defines a
#'   "decision" (e.g. p = 0.05, or zero for effect sizes). Required
#'   for the fragility index.
#' @param n_boot number of bootstrap perturbations for the fragility
#'   index search (default 200).
#' @param seed optional integer seed.
#'
#' @return A `crossstatr_fragility` object.
#'
#' @section Methodology:
#'
#'   The function systematically perturbs analysis decisions and
#'   measures the change in the result:
#'
#'   - **Leave-one-replicate-out**: recompute the statistic excluding
#'     each replicate in turn. The fragility index is the number of
#'     replicates whose exclusion flips the decision.
#'
#'   - **Outlier handling**: recompute under each outlier rule (none,
#'     IQR, MAD). The fragility index is 1 if any rule flips the
#'     decision.
#'
#'   - **Baseline correction**: recompute under each baseline method.
#'     The fragility index is 1 if any baseline flips the decision.
#'
#'   - **Model selection**: recompute under each alternative model.
#'
#'   - **Measurement uncertainty**: progressively inject Gaussian noise
#'     at increasing SD, recompute, and find the smallest SD at which
#'     the decision flips. This is the **noise fragility threshold**.
#'
#'   The overall fragility score is
#'   \deqn{F = \sum_k w_k F_k,}
#'   where \eqn{F_k \in [0, 1]} is the fraction of perturbations of
#'   type \eqn{k} that flip the decision, and \eqn{w_k} are user-
#'   supplied weights (uniform by default).
#'
#'   The function reports the original result, the range across all
#'   perturbations, and the most influential decision (the one that
#'   produces the largest swing).
#'
#' @section Comparison with established alternatives:
#'
#'   - **Fragility index** (Walsh et al., 2014) in clinical trials:
#'     counts the number of events that would need to change to flip
#'     statistical significance. This function generalizes that idea
#'     to continuous outcomes and to a wider range of analysis
#'     decisions.
#'   - **Multiverse analysis** (Steegen et al., 2016): a general
#'     framework for exploring analysis-decision sensitivity. This
#'     function is a focused, quantitative implementation for
#'     experimental science.
#'
#' @section Limitations:
#'
#'   - This is exploratory sensitivity analysis, not confirmatory
#'     inference. The reported fragility score does not have a known
#'     reference distribution; do not treat it as a p-value.
#'   - The fragility index counts "decision flips" with respect to a
#'     user-supplied threshold. The threshold itself is a modeling
#'     decision and should be reported alongside the score.
#'   - The measurement-uncertainty injection assumes additive Gaussian
#'     noise; multiplicative or correlated noise is not modelled.
#'   - For continuous statistics (e.g. effect sizes), use the range
#'     and the most-influential-decision summary; the binary fragility
#'     index is less informative.
#'
#' @references
#'   - Walsh, M., Srinathan, S. K., McAuley, D. F. et al. (2014). The
#'     statistical significance of randomized controlled trial results is
#'     frequently fragile: a case for a Fragility Index. Journal of
#'     Clinical Epidemiology, 67, 622-628.
#'     \doi{10.1016/j.jclinepi.2013.10.019}.
#'
#'   - Steegen, S., Tuerlinckx, F., Gelman, A. and Vanpaemel, W. (2016).
#'     Increasing Transparency Through a Multiverse Analysis. Perspectives
#'     on Psychological Science, 11, 702-712.
#'     \doi{10.1177/1745691616658637}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = c(rnorm(15, mean = 1), rnorm(15, mean = 1.5)),
#'   grp = factor(rep(c("A", "B"), each = 15)),
#'   rep = factor(rep(1:15, times = 2)))
#' # Analysis function: returns mean difference
#' af <- function(opts) {
#'   d <- opts$data
#'   if (!is.null(opts$rep_lo)) d <- d[d$rep != opts$rep_lo, ]
#'   diff(tapply(d$y, d$grp, mean))
#' }
#' res <- experiment_fragility(
#'   analysis_fn = af,
#'   data = df,
#'   replicates = "rep",
#'   decision_threshold = 0)
#' print(res)
#'
#' @export
experiment_fragility <- function(analysis_fn, data, replicates = NULL,
                                  outlier_methods = c("none", "iqr", "mad"),
                                  baseline_methods = c("none"),
                                  model_alternatives = NULL,
                                  measurement_uncertainties = NULL,
                                  decision_threshold,
                                  n_boot = 200L, seed = NULL) {
  call <- match.call()
  if (!is.function(analysis_fn)) {
    stop("'analysis_fn' must be a function.", call. = FALSE)
  }
  if (!is.numeric(decision_threshold) || length(decision_threshold) != 1L) {
    stop("'decision_threshold' must be a single numeric value.", call. = FALSE)
  }
  .set_seed(seed)

  # Original result
  orig <- analysis_fn(list(data = data))
  orig_decision <- orig > decision_threshold

  # Leave-one-replicate-out
  loro_results <- NULL
  loro_fragility <- NA_real_
  loro_swings <- NA_real_
  if (!is.null(replicates)) {
    if (!(replicates %in% colnames(data))) {
      stop(sprintf("Replicate column '%s' not in data.", replicates),
           call. = FALSE)
    }
    reps <- unique(data[[replicates]])
    loro_results <- vapply(reps, function(r) {
      d <- data[data[[replicates]] != r, , drop = FALSE]
      analysis_fn(list(data = d, rep_lo = r))
    }, numeric(1L))
    flips <- (loro_results > decision_threshold) != orig_decision
    loro_fragility <- mean(flips)
    loro_swings <- max(abs(loro_results - orig))
  }

  # Outlier method fragility
  outlier_results <- lapply(outlier_methods, function(om) {
    d <- data
    if (om != "none" && "y" %in% colnames(d)) {
      y <- d$y
      if (om == "iqr") {
        q <- stats::quantile(y, c(0.25, 0.75))
        iqr <- q[2L] - q[1L]
        keep <- y >= q[1L] - 1.5 * iqr & y <= q[2L] + 1.5 * iqr
        d <- d[keep, , drop = FALSE]
      } else if (om == "mad") {
        m <- stats::median(y); md <- stats::mad(y)
        keep <- abs(y - m) < 3.5 * md
        d <- d[keep, , drop = FALSE]
      }
    }
    analysis_fn(list(data = d, outlier_method = om))
  })
  names(outlier_results) <- outlier_methods
  outlier_flips <- vapply(outlier_results, function(v)
    (v > decision_threshold) != orig_decision, logical(1L))
  outlier_fragility <- mean(outlier_flips)
  outlier_swings <- max(abs(unlist(outlier_results) - orig))

  # Baseline method fragility (only if analysis_fn appears to support it)
  baseline_results <- lapply(baseline_methods, function(bm) {
    res <- tryCatch(analysis_fn(list(data = data, baseline = bm)),
                    error = function(e) NA_real_)
    res
  })
  names(baseline_results) <- baseline_methods
  baseline_res_ok <- !vapply(baseline_results, is.na, logical(1L))
  baseline_flips <- vapply(baseline_results, function(v) {
    if (is.na(v)) return(NA)
    (v > decision_threshold) != orig_decision
  }, logical(1L))
  baseline_fragility <- if (any(baseline_res_ok))
    mean(baseline_flips[!is.na(baseline_flips)]) else NA_real_

  # Model alternative fragility
  model_results <- NULL
  model_fragility <- NA_real_
  if (!is.null(model_alternatives)) {
    model_results <- lapply(model_alternatives, function(m) {
      tryCatch(analysis_fn(list(data = data, model = m)),
               error = function(e) NA_real_)
    })
    model_flips <- vapply(model_results, function(v) {
      if (is.na(v)) return(FALSE)
      (v > decision_threshold) != orig_decision
    }, logical(1L))
    model_fragility <- mean(model_flips)
  }

  # Noise fragility threshold
  noise_threshold <- NA_real_
  if (!is.null(measurement_uncertainties)) {
    flips_at_sd <- vapply(measurement_uncertainties, function(sd_inject) {
      injected <- rep(FALSE, n_boot)
      for (b in seq_len(n_boot)) {
        d_noisy <- data
        if ("y" %in% colnames(d_noisy)) {
          d_noisy$y <- d_noisy$y + rnorm(nrow(d_noisy), sd = sd_inject)
        }
        res_b <- tryCatch(analysis_fn(list(data = d_noisy)),
                          error = function(e) NA_real_)
        if (is.na(res_b)) next
        injected[b] <- (res_b > decision_threshold) != orig_decision
      }
      mean(injected)
    }, numeric(1L))
    # Smallest SD where >= 50% of perturbations flip the decision
    flip_sd <- measurement_uncertainties[flips_at_sd >= 0.5]
    if (length(flip_sd) > 0L) noise_threshold <- min(flip_sd)
  }

  # Overall fragility score (uniform weights among available components)
  components <- c(
    if (!is.na(loro_fragility)) loro_fragility,
    outlier_fragility,
    if (!is.na(baseline_fragility)) baseline_fragility,
    if (!is.na(model_fragility)) model_fragility
  )
  overall <- if (length(components) > 0L) mean(components) else NA_real_

  # Most influential decision
  swing_df <- data.frame(
    decision = c("leave_one_replicate", "outlier_method",
                 "baseline_method", "model_alternative"),
    swing = c(
      if (!is.null(loro_results)) loro_swings else NA_real_,
      outlier_swings,
      if (!is.na(baseline_fragility)) max(abs(unlist(baseline_results) - orig),
                                          na.rm = TRUE) else NA_real_,
      if (!is.null(model_results)) max(abs(unlist(model_results) - orig),
                                       na.rm = TRUE) else NA_real_
    ),
    stringsAsFactors = FALSE)
  swing_df <- swing_df[order(-swing_df$swing), , drop = FALSE]

  .crossstatr_result(
    list(original_result = orig,
         original_decision = orig_decision,
         decision_threshold = decision_threshold,
         leave_one_replicate = loro_results,
         outlier_methods = unlist(outlier_results),
         baseline_methods = unlist(baseline_results[baseline_res_ok]),
         model_alternatives = unlist(model_results),
         noise_fragility_threshold = noise_threshold,
         component_fragility = c(
           replicate = loro_fragility,
           outlier = outlier_fragility,
           baseline = baseline_fragility,
           model = model_fragility),
         overall_fragility = overall,
         most_influential = swing_df,
         n_boot = n_boot,
         method = "Experimental fragility analysis",
         call = call),
    "crossstatr_fragility")
}

#' @export
print.crossstatr_fragility <- function(x, ...) {
  cat("<crossStatR::experiment_fragility> Sensitivity summary\n")
  cat(sprintf("Original result: %.4g (threshold: %.4g, decision: %s)\n",
              x$original_result, x$decision_threshold,
              if (x$original_decision) "above" else "below"))
  cat("\nComponent fragility (fraction flipping decision):\n")
  for (nm in names(x$component_fragility)) {
    v <- x$component_fragility[[nm]]
    cat(sprintf("  %-22s %s\n", nm, if (is.na(v)) "n/a" else sprintf("%.2f", v)))
  }
  cat(sprintf("\nOverall fragility score: %s\n",
              if (is.na(x$overall_fragility)) "n/a"
              else sprintf("%.3f", x$overall_fragility)))
  cat("\nMost influential decisions:\n")
  print(x$most_influential, row.names = FALSE)
  invisible(x)
}
