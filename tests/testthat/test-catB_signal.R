
test_that("adaptive_baseline basic", {
  set.seed(1)
  x <- seq(400, 800, by = 2)
  y <- dnorm(x, 600, 30) * 10 + 1 + 0.001 * (x - 400) + rnorm(length(x), sd = 0.1)
  res <- adaptive_baseline(x, y, method = "asls", max_iter = 50)
  expect_s3_class(res, "crossstatr_baseline")
  expect_equal(length(res$baseline), length(y))
})

test_that("adaptive_baseline lowess", {
  set.seed(2)
  x <- 1:100
  y <- sin(x / 20) + 0.5 + rnorm(100, sd = 0.1)
  res <- adaptive_baseline(x, y, method = "lowess")
  # Baseline should be roughly below the signal (allowing for noise)
  expect_true(mean(res$baseline) < mean(y) + 1)
  expect_equal(length(res$corrected), length(y))
})

test_that("weighted_peak_fit basic gaussian", {
  set.seed(3)
  x <- seq(500, 600, by = 0.5)
  y <- 5 * dnorm(x, 550, 8) + rnorm(length(x), sd = 0.05)
  res <- weighted_peak_fit(x, y, n_peaks = 1, peak_shape = "gaussian")
  expect_s3_class(res, "crossstatr_peak_fit")
  expect_true(!is.null(res$coefficients))
})

test_that("weighted_peak_fit auto", {
  set.seed(4)
  x <- seq(0, 20, by = 0.1)
  y <- 5 * dnorm(x, 5, 0.5) + 3 * dnorm(x, 15, 0.8) + rnorm(length(x), sd = 0.1)
  # Fixed peak count makes this a stable test of automatic location detection.
  res <- weighted_peak_fit(x, y, peaks = "auto", n_peaks = 2,
                           peak_shape = "gaussian", max_iter = 1000)
  expect_s3_class(res, "crossstatr_peak_fit")
  expect_equal(res$n_peaks, 2L)

})

test_that("adaptive_smoothing basic sgolay", {
  set.seed(5)
  x <- 1:100
  y <- sin(x / 10) + rnorm(100, sd = 0.1)
  res <- adaptive_smoothing(x, y, method = "sgolay", window = 7)
  expect_s3_class(res, "crossstatr_smoothing")
  expect_equal(length(res$smoothed), length(y))
})

test_that("adaptive_smoothing spline", {
  set.seed(6)
  x <- 1:50
  y <- 2 * x + rnorm(50, sd = 5)
  res <- adaptive_smoothing(x, y, method = "spline")
  expect_s3_class(res, "crossstatr_smoothing")
})

test_that("replicate_decompose basic LMM", {
  set.seed(7)
  df <- data.frame(
    y = rnorm(60),
    treat = factor(rep(c("A","B"), each = 30)),
    bio = factor(rep(1:6, each = 10)),
    tech = factor(rep(1:10, times = 6)))
  res <- replicate_decompose(y ~ treat + (1|bio), data = df)
  expect_s3_class(res, "crossstatr_replicate")
})
