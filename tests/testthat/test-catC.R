
test_that("signal_difference_regions basic", {
  set.seed(1)
  x <- seq(400, 800, length.out = 50)
  y1 <- matrix(replicate(8, sin((x - 400) * pi / 200) * 0.2 + rnorm(50, sd = 0.1)),
               nrow = 8, byrow = TRUE)
  y2 <- matrix(replicate(8, {
    base <- sin((x - 400) * pi / 200) * 0.2
    base + dnorm(x, 600, 15) + rnorm(50, sd = 0.1)
  }), nrow = 8, byrow = TRUE)
  res <- signal_difference_regions(x, y1, y2, n_perm = 50, seed = 1)
  expect_s3_class(res, "crossstatr_sdr")
  expect_equal(length(res$x_grid), 50L)
})

test_that("signal_difference_regions Bonferroni (n_perm=0)", {
  set.seed(2)
  x <- seq(0, 10, length.out = 30)
  y1 <- matrix(replicate(5, x + rnorm(30, sd = 0.2)), nrow = 5, byrow = TRUE)
  y2 <- matrix(replicate(5, x + 0.5 + rnorm(30, sd = 0.2)), nrow = 5, byrow = TRUE)
  res <- signal_difference_regions(x, y1, y2, n_perm = 0)
  expect_s3_class(res, "crossstatr_sdr")
})

test_that("peak_evidence basic", {
  set.seed(3)
  x <- seq(400, 800, by = 0.5)
  y <- dnorm(x, 550, 8) * 5 + 0.1 + rnorm(length(x), sd = 0.05)
  res <- peak_evidence(x, y, n_perm = 50, seed = 1)
  expect_s3_class(res, "crossstatr_peak_evidence")
  expect_true(nrow(res$peaks) >= 1L)
})

test_that("peak_evidence with explicit peak_index", {
  set.seed(4)
  x <- seq(0, 20, by = 0.1)
  y <- dnorm(x, 10, 0.5) * 3 + rnorm(length(x), sd = 0.1)
  res <- peak_evidence(x, y, peak_index = 100L, n_perm = 20, seed = 1)
  expect_equal(nrow(res$peaks), 1L)
})

test_that("calibration_transfer_risk basic", {
  set.seed(5)
  x <- seq(0.1, 10, length.out = 20)
  y_s <- 2 + 3 * x + rnorm(20, sd = 0.2)
  y_t <- 2.5 + 2.7 * x + rnorm(20, sd = 0.25)
  source_cal <- lm(y_s ~ x)
  res <- calibration_transfer_risk(source_cal,
                                    data.frame(x = x, y = y_s),
                                    data.frame(x = x, y = y_t),
                                    n_boot = 50, seed = 1)
  expect_s3_class(res, "crossstatr_cal_transfer_risk")
})

test_that("kinetic_discrimination basic", {
  set.seed(6)
  m1 <- function(p, t) p[1] * exp(-p[2] * t)
  m2 <- function(p, t) p[1] * (1 + p[3] * t) * exp(-p[2] * t)
  res <- kinetic_discrimination(
    models = list(simple = m1, extended = m2),
    starts = list(c(A = 1, k = 0.1), c(A = 1, k = 0.1, alpha = 0.05)),
    true_model = 1L, true_params = c(1, 0.1),
    times = seq(0, 30, length.out = 10), sigma = 0.05, n_sim = 20, seed = 1)
  expect_s3_class(res, "crossstatr_kinetic_discrimination")
})

test_that("experiment_fragility basic", {
  set.seed(7)
  df <- data.frame(
    y = c(rnorm(10, 1), rnorm(10, 1.5)),
    grp = factor(rep(c("A","B"), each = 10)),
    rep = factor(rep(1:10, times = 2)))
  af <- function(opts) {
    d <- opts$data
    if (!is.null(opts$rep_lo)) d <- d[d$rep != opts$rep_lo, ]
    diff(tapply(d$y, d$grp, mean))
  }
  res <- experiment_fragility(
    analysis_fn = af, data = df, replicates = "rep",
    decision_threshold = 0, n_boot = 30, seed = 1)
  expect_s3_class(res, "crossstatr_fragility")
})

test_that("drift_uncertainty basic", {
  set.seed(8)
  t <- seq(0, 10, length.out = 50)
  y <- 2 + 1.5 * t + 0.3 * t + rnorm(50, sd = 0.2)
  res <- drift_uncertainty(y, time = t, drift_model = "linear",
                            n_boot = 50, seed = 1)
  expect_s3_class(res, "crossstatr_drift_uncertainty")
})

test_that("detection_reliability basic", {
  set.seed(9)
  x <- seq(0.1, 10, length.out = 20)
  y <- 0.5 + 0.9 * x + rnorm(20, sd = 0.3)
  cal <- lm(y ~ x)
  res <- detection_reliability(cal, n_boot = 50, seed = 1)
  expect_s3_class(res, "crossstatr_detection_reliability")
  expect_equal(nrow(res$curve), 50L)
})

test_that("uncertainty_origin basic additive", {
  stages <- c(sampling = 0.05, instrument = 0.10,
              calibration = 0.02, drift = 0.03)
  res <- uncertainty_origin(stages, n_mc = 200, seed = 1)
  expect_s3_class(res, "crossstatr_uncertainty_origin")
  expect_equal(sum(res$percent_additive), 1, tolerance = 1e-9)
})

test_that("uncertainty_origin with correlations", {
  stages <- c(A = 1, B = 1)
  res <- uncertainty_origin(stages,
    correlations = list(A__B = 0.5), n_mc = 500, seed = 1)
  expect_s3_class(res, "crossstatr_uncertainty_origin")
  expect_gt(res$total_variance_correlated, res$total_variance_additive)
})
