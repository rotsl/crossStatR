
test_that("drift_correct_stat linear", {
  set.seed(1)
  t <- seq(0, 100, length.out = 30)
  y <- 5 + 0.02 * t + rnorm(30, sd = 0.1)
  res <- drift_correct_stat(y, time = t, drift_model = "linear")
  expect_s3_class(res, "crossstatr_drift")
  expect_equal(length(res$y_corrected), length(y))
})

test_that("drift_correct_stat polynomial", {
  set.seed(2)
  t <- seq(0, 50, length.out = 30)
  y <- 10 + 0.5 * t - 0.005 * t^2 + rnorm(30, sd = 0.2)
  res <- drift_correct_stat(y, time = t, drift_model = "polynomial", degree = 2)
  expect_s3_class(res, "crossstatr_drift")
})

test_that("robust_lod basic", {
  set.seed(3)
  x <- seq(0, 10, length.out = 20)
  y <- 2 + 0.5 * x + rnorm(20, sd = 0.2)
  cal <- lm(y ~ x)
  res <- robust_lod(cal, method = "currie")
  expect_s3_class(res, "crossstatr_lod")
  expect_true(res$lod > 0)
})

test_that("robust_lod hubaux_vos", {
  set.seed(4)
  x <- seq(0, 5, length.out = 25)
  y <- 0.5 + 0.8 * x + rnorm(25, sd = 0.15)
  cal <- lm(y ~ x)
  res <- robust_lod(cal, method = "hubaux_vos")
  expect_true(res$lod > 0)
})

test_that("calibration_compare basic", {
  set.seed(5)
  x1 <- seq(0, 10, length.out = 20)
  y1 <- 1 + 0.5 * x1 + rnorm(20, sd = 0.2)
  x2 <- seq(0, 10, length.out = 20)
  y2 <- 1.2 + 0.45 * x2 + rnorm(20, sd = 0.2)
  d1 <- data.frame(x = x1, y = y1)
  d2 <- data.frame(x = x2, y = y2)
  res <- calibration_compare(d1, d2)
  expect_s3_class(res, "crossstatr_calibration_compare")
})

test_that("curve_similarity l2", {
  set.seed(6)
  x <- seq(0, 10, length.out = 50)
  y1 <- sin(x)
  y2 <- sin(x) + 0.05
  res <- curve_similarity(x, y1, x, y2, method = "l2")
  expect_s3_class(res, "crossstatr_curve_similarity")
})

test_that("curve_similarity correlation", {
  set.seed(7)
  x <- seq(0, 10, length.out = 30)
  y1 <- x + rnorm(30, sd = 0.5)
  y2 <- x + rnorm(30, sd = 0.5)
  res <- curve_similarity(x, y1, x, y2, method = "correlation")
  expect_true(!is.na(res$distance))
})

test_that("noise_characterize basic", {
  set.seed(8)
  x <- seq(0, 10, length.out = 100)
  y <- sin(x) + rnorm(100, sd = 0.1)
  res <- noise_characterize(x, y)
  expect_s3_class(res, "crossstatr_noise")
  expect_true(res$noise_sd > 0)
})

test_that("batch_effect_quantify basic", {
  set.seed(9)
  df <- data.frame(
    y = rnorm(60),
    treat = factor(rep(c("A","B"), each = 30)),
    batch = factor(rep(1:3, each = 20)))
  res <- batch_effect_quantify(y ~ treat, data = df, batch = "batch")
  expect_s3_class(res, "crossstatr_batch")
})

test_that("batch_effect_quantify detects confounding", {
  df <- data.frame(
    y = rnorm(20),
    treat = factor(rep(c("A","B"), each = 10)),
    batch = factor(rep(c("X","Y"), each = 10)))  # perfectly confounded
  expect_error(batch_effect_quantify(y ~ treat, data = df, batch = "batch"),
               "confound")
})
