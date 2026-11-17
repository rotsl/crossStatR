#' Optimise sampling timepoints for an experiment
#'
#' Given a kinetic model `model_func(params, t) -> y`, chooses additional
#' sampling timepoints that maximise the expected information about the
#' parameters, given a set of fixed (already scheduled) times. Several
#' optimal-design criteria are supported: D-optimal (maximise the
#' determinant of the information matrix), A-optimal (minimise the trace of
#' its inverse), C-optimal (minimise the variance of a specific linear
#' combination of parameters), and a model-discrimination criterion that
#' maximises the detectable difference between two competing models.
#'
#' The search is greedy: starting from `fixed_times`, the function picks the
#' candidate that most improves the criterion, repeats for `n_extra` points,
#' and returns the augmented design. Greedy is fast but not guaranteed to
#' find the global optimum; it is a good practical approximation when the
#' number of points is small.
#'
#' @param model_func function `f(params, t) -> y` returning predicted
#'   responses.
#' @param params numeric vector of (assumed true) parameter values.
#' @param t_range numeric length-2 vector giving the time interval.
#' @param n_candidates integer; number of points in the candidate grid.
#' @param criterion `"d_optimal"`, `"a_optimal"`, `"c_optimal"`, or
#'   `"model_discrimination"`.
#' @param fixed_times numeric vector of already-fixed sampling times.
#' @param n_extra integer; number of additional points to choose.
#' @param conf_level confidence level (stored on the result).
#' @param competitor_model required for `criterion = "model_discrimination"`:
#'   a function `f(params, t) -> y` for the competitor model.
#' @param competitor_params parameter values for the competitor model.
#' @param c_vector required for `criterion = "c_optimal"`: a numeric vector
#'   of length `length(params)` defining the linear combination
#'   \eqn{c'\theta} whose variance is to be minimised.
#'
#' @return A `crossstatr_timepoint` object with elements:
#'   \item{optimal_times}{numeric vector of chosen sampling times
#'     (fixed_times + extra).}
#'   \item{criterion_value}{value of the criterion at the chosen design
#'     (for D-optimal, det(M); for A-optimal, -trace(M^-1); for
#'     C-optimal, -c'M^-1 c; for model_discrimination, the squared
#'     prediction difference scaled by total sensitivity).}
#'   \item{sensitivity_matrix}{the matrix `G` of dy/dparam at the chosen
#'     times.}
#'   \item{information_matrix}{M = G'G (returned for D/A/C criteria).}
#'   \item{criterion}{criterion label.}
#'   \item{model_label}{label identifying the model(s) used.}
#'   \item{call}{the matched call.}
#'
#' @section Methodology:
#'
#'   Let \eqn{G_{ij} = \partial f_i / \partial \theta_j} be the
#'   sensitivity matrix evaluated at the candidate times. The Fisher
#'   information matrix under a Gaussian noise model with unit variance is
#'   \deqn{M = G'G.}
#'
#'   D-optimal design maximises \eqn{\det(M)}, which corresponds to
#'   minimising the volume of the joint parameter confidence ellipsoid
#'   (Fedorov, 1972).
#'
#'   A-optimal design minimises \eqn{\mathrm{tr}(M^{-1})}, the sum of the
#'   parameter variances (Kiefer, 1959).
#'
#'   C-optimal design minimises \eqn{c'M^{-1}c}, the variance of a specific
#'   linear combination \eqn{c'\theta} (Silvey, 1980).
#'
#'   Model discrimination: The package maximises
#'   \deqn{\frac{\sum_i (f_1(\theta_1, t_i) - f_2(\theta_2, t_i))^2}
#'              {\sum_i (G_{1i}^2 + G_{2i}^2) + \epsilon}}
#'   where the numerator captures how far apart the models are at the
#'   chosen times and the denominator approximates the squared standard
#'   error of the difference (Lopez-Fidalgo and Tommasi, 2007). Higher
#'   values are better.
#'
#'   Sensitivities are computed by central finite differences via
#'   `.num_jacobian`.
#'
#' @section Limitations:
#'
#'   - D-optimal search completes rank first when the starting information
#'     matrix is singular. Candidates are ordered by rank, then by the product
#'     of positive squared singular values. A warning identifies this phase;
#'     the reported final criterion remains the determinant. A final singular
#'     design produces a warning and has no finite D-optimal criterion.
#'
#'   - Local optimality: the criteria are evaluated at the supplied `params`
#'     only. Bayesian or minimax designs that integrate over a prior on
#'     `params` are more robust when the assumed parameter values are
#'     uncertain; this function does not implement them.
#'   - The greedy search may not find the global optimum. For `n_extra <= 4`
#'     and small candidate grids, an exhaustive search would be feasible but
#'     is not currently implemented.
#'   - The information matrix assumes unit noise variance; for non-unit
#'     variance the criteria scale by `1 / sigma^2` and the relative
#'     ordering of designs is unchanged, so the same `optimal_times` would
#'     be selected.
#'   - For `a_optimal` and `c_optimal`, when `M` is singular (which happens
#'     early in the greedy search when fewer timepoints are chosen than
#'     there are parameters), the method falls back to the Moore-Penrose
#'     pseudo-inverse `MASS::ginv(M)`. The criterion is still well-defined
#'     and guides the search toward full-rank designs, but the literal
#'     interpretation as `trace(M^{-1})` requires a non-singular `M`.
#'
#' @references
#'   - Fedorov, V. V. (1972). Theory of Optimal Experiments. Academic
#'     Press.
#'     \url{https://shop.elsevier.com/books/theory-of-optimal-experiments/fedorov/978-0-12-250750-2}.
#'
#'   - Atkinson, A. C., Donev, A. N. and Tobias, R. D. (2007). Optimum
#'     Experimental Designs, with SAS. Oxford University Press.
#'     \doi{10.1093/oso/9780199296590.001.0001}.
#'
#' @examples
#' # First-order decay model f(p, t) = A * exp(-k * t)
#' model <- function(p, t) p[1] * exp(-p[2] * t)
#' params <- c(A = 5, k = 0.3)
#'
#' # D-optimal: 3 extra timepoints in [0, 20]
#' r <- timepoint_optimize(model, params, t_range = c(0, 20),
#'                         n_candidates = 20L, criterion = "d_optimal",
#'                         fixed_times = c(0, 20), n_extra = 3L)
#' print(r)
#'
#' # C-optimal: minimise variance of the parameter 'k' (second parameter)
#' timepoint_optimize(model, params, t_range = c(0, 20),
#'                    criterion = "c_optimal",
#'                    c_vector = c(0, 1), n_extra = 3L)
#'
#' # Model discrimination between first-order and saturation
#' sat <- function(p, t) p[1] * t / (p[2] + t)
#' timepoint_optimize(model, params, t_range = c(0, 20),
#'                    criterion = "model_discrimination",
#'                    competitor_model = sat,
#'                    competitor_params = c(Vmax = 5, Km = 3),
#'                    n_extra = 3L)
#'
#' @export
timepoint_optimize <- function(model_func, params, t_range = c(0, 100),
                                n_candidates = 10L,
                                criterion = c("d_optimal", "a_optimal",
                                              "c_optimal",
                                              "model_discrimination"),
                                fixed_times = NULL, n_extra = 3L,
                                conf_level = 0.95,
                                competitor_model = NULL,
                                competitor_params = NULL,
                                c_vector = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  if (!is.function(model_func)) {
    stop("'model_func' must be a function.", call. = FALSE)
  }
  .assert_numeric_vector(params, "params", min_length = 1L)
  .assert_numeric_vector(t_range, "t_range", min_length = 2L)
  if (length(t_range) != 2L || t_range[1L] >= t_range[2L]) {
    stop("'t_range' must be c(t_min, t_max) with t_min < t_max.",
         call. = FALSE)
  }
  .assert_positive_integer(n_candidates, "n_candidates")
  .assert_positive_integer(n_extra, "n_extra")
  .assert_probability(conf_level, "conf_level")
  if (!is.null(fixed_times)) {
    if (!is.numeric(fixed_times)) {
      stop("'fixed_times' must be a numeric vector.", call. = FALSE)
    }
    fixed_times <- sort(unique(fixed_times))
  }
  if (criterion == "c_optimal" && is.null(c_vector)) {
    stop("c_optimal requires 'c_vector'.", call. = FALSE)
  }
  if (criterion == "c_optimal" &&
      (!is.numeric(c_vector) || length(c_vector) != length(params))) {
    stop("'c_vector' must be a numeric vector of length length(params).",
         call. = FALSE)
  }
  if (criterion == "model_discrimination") {
    if (is.null(competitor_model) || !is.function(competitor_model)) {
      stop("model_discrimination requires 'competitor_model' (a function).",
           call. = FALSE)
    }
    if (is.null(competitor_params) ||
        !is.numeric(competitor_params)) {
      stop("model_discrimination requires numeric 'competitor_params'.",
           call. = FALSE)
    }
  }

  candidates <- seq(t_range[1L], t_range[2L], length.out = n_candidates)
  if (!is.null(fixed_times) && length(fixed_times) > 0L) {
    # Drop candidates that are essentially identical to fixed_times
    is_fixed <- vapply(candidates, function(c) {
      any(abs(c - fixed_times) < 1e-8)
    }, logical(1L))
    candidates <- candidates[!is_fixed]
  }

  sensitivity_at <- function(times, p, f) {
    .num_jacobian(p, function(pp) f(pp, times))
  }

  criterion_at <- function(times) {
    if (length(times) == 0L) return(-Inf)
    if (criterion == "model_discrimination") {
      pred1 <- model_func(params, times)
      pred2 <- competitor_model(competitor_params, times)
      G1 <- sensitivity_at(times, params, model_func)
      G2 <- sensitivity_at(times, competitor_params, competitor_model)
      num <- sum((pred1 - pred2)^2)
      den <- sum(G1^2) + sum(G2^2) + 1e-10
      return(num / den)
    }
    G <- sensitivity_at(times, params, model_func)
    M <- t(G) %*% G
    if (criterion == "d_optimal") {
      d <- tryCatch(det(M), error = function(e) NA_real_)
      return(if (is.finite(d) && d > 0) d else -Inf)
    } else if (criterion == "a_optimal") {
      # Pseudo-inverse handles rank-deficient M (early in greedy search)
      M_inv <- tryCatch(solve(M), error = function(e) NULL)
      if (is.null(M_inv)) M_inv <- MASS::ginv(M)
      v <- -sum(diag(M_inv))
      return(if (is.finite(v)) v else -Inf)
    } else if (criterion == "c_optimal") {
      M_inv <- tryCatch(solve(M), error = function(e) NULL)
      if (is.null(M_inv)) M_inv <- MASS::ginv(M)
      v <- -as.numeric(t(c_vector) %*% M_inv %*% c_vector)
      return(if (is.finite(v)) v else -Inf)
    }
  }

  chosen <- if (!is.null(fixed_times) && length(fixed_times) > 0L) {
    as.numeric(fixed_times)
  } else {
    numeric(0)
  }
  remaining <- candidates

  n_params <- length(params)
  design_rank <- function(times) {
    if (!length(times)) return(0L)
    qr(sensitivity_at(times, params, model_func))$rank
  }
  rank_completion <- criterion == "d_optimal" && design_rank(chosen) < n_params
  if (rank_completion) {
    warning("D-optimal search starts with a rank-deficient design; ",
            "rank completion precedes determinant optimisation.", call. = FALSE)
  }
  for (k in seq_len(n_extra)) {
    if (length(remaining) == 0L) break
    use_fallback <- criterion == "d_optimal" && design_rank(chosen) < n_params
    vals <- vapply(remaining, function(t) {
      if (use_fallback) {
        G <- sensitivity_at(c(chosen, t), params, model_func)
        singular <- svd(G)$d
        positive <- singular[singular > max(singular, 1) * 1e-7]
        if (!length(positive)) return(-Inf)
        sum(log(positive^2))
      } else criterion_at(c(chosen, t))
    }, numeric(1L))
    if (all(!is.finite(vals))) break
    if (use_fallback) {
      ranks <- vapply(remaining, function(t) design_rank(c(chosen, t)), integer(1L))
      vals[ranks < max(ranks)] <- -Inf
    }
    best_idx <- which.max(vals)
    chosen <- c(chosen, remaining[best_idx])
    remaining <- remaining[-best_idx]
  }
  chosen <- sort(unique(chosen))
  criterion_value <- criterion_at(chosen)
  if (criterion == "d_optimal" && design_rank(chosen) < n_params) {
    warning("The final design remains rank deficient; more informative times are required.", call. = FALSE)
  }

  if (criterion == "model_discrimination") {
    G_final <- sensitivity_at(chosen, params, model_func)
    M_final <- NULL
    G2_final <- sensitivity_at(chosen, competitor_params, competitor_model)
    information_matrix <- list(model_1 = t(G_final) %*% G_final,
                              model_2 = t(G2_final) %*% G2_final)
    model_label <- "user_model vs competitor"
  } else {
    G_final <- sensitivity_at(chosen, params, model_func)
    M_final <- t(G_final) %*% G_final
    information_matrix <- M_final
    model_label <- "user_model"
  }

  tidy_df <- data.frame(
    time = chosen,
    stringsAsFactors = FALSE
  )

  .crossstatr_result(
    list(optimal_times = chosen,
         criterion_value = criterion_value,
         sensitivity_matrix = G_final,
         information_matrix = information_matrix,
         criterion = criterion,
         rank_completion = rank_completion,
         model_label = model_label,
         fixed_times = fixed_times,
         n_extra = n_extra,
         tidy = tidy_df,
         call = call),
    "crossstatr_timepoint"
  )
}

#' @export
print.crossstatr_timepoint <- function(x, ...) {
  cat(sprintf("<crossStatR::timepoint_optimize> criterion = %s\n",
              x$criterion))
  cat(sprintf("optimal_times = %s\n",
              paste(sprintf("%.4g", x$optimal_times), collapse = ", ")))
  cat(sprintf("criterion_value = %.4g\n", x$criterion_value))
  invisible(x)
}
