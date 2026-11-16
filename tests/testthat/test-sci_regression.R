
test_that("sci_regression basic OLS", {
  set.seed(1)
  df <- data.frame(x = 1:20, y = 2 * (1:20) + rnorm(20, sd = 0.5))
  res <- sci_regression(y ~ x, data = df)
  expect_s3_class(res, "crossstatr_regression")
  expect_equal(res$fit_type, "ols")
  # Compare to stats::lm
  expect_equal(coef(res$model), coef(lm(y ~ x, data = df)), tolerance = 1e-9)
})

test_that("sci_regression weighted", {
  set.seed(2)
  df <- data.frame(x = 1:20, y = (1:20) + rnorm(20, sd = 1))
  w <- runif(20, 0.5, 2)
  res <- sci_regression(y ~ x, data = df, weights = w)
  expect_equal(res$fit_type, "wls")
})

test_that("sci_regression robust", {
  set.seed(3)
  df <- data.frame(x = 1:30, y = 3 * (1:30) + rnorm(30, sd = 2))
  # Add an outlier
  df$y[5] <- 100
  res <- sci_regression(y ~ x, data = df, robust = TRUE)
  expect_equal(res$fit_type, "rlm")
})

test_that("sci_regression nonlinear", {
  set.seed(4)
  x <- seq(0, 10, length.out = 30)
  y <- 5 * exp(-0.3 * x) + rnorm(30, sd = 0.2)
  df <- data.frame(x = x, y = y)
  res <- sci_regression(y ~ A * exp(-k * x), data = df,
                        nonlinear = TRUE, start = list(A = 5, k = 0.3))
  expect_equal(res$fit_type, "nls")
  expect_true(!is.null(res$coefficients))
})

test_that("sci_regression errors with bad start for nls", {
  df <- data.frame(x = 1:5, y = 1:5)
  expect_error(sci_regression(y ~ a * x, data = df, nonlinear = TRUE),
              "starting")
})

test_that("sci_correlation basic pearson", {
  set.seed(5)
  x <- rnorm(50)
  y <- 2 * x + rnorm(50, sd = 1)
  res <- sci_correlation(x, y)
  expect_s3_class(res, "crossstatr_correlation")
  expect_true(abs(res$correlation - cor(x, y)) < 1e-9)
})

test_that("sci_correlation spearman", {
  set.seed(6)
  x <- rnorm(30)
  y <- rank(x) + rnorm(30, sd = 0.5)
  res <- sci_correlation(x, y, method = "spearman")
  expect_true(!is.na(res$correlation))
})

test_that("sci_correlation handles NA", {
  x <- c(1, 2, NA, 4, 5)
  y <- c(2, 4, 6, 8, 10)
  res <- sci_correlation(x, y)
  expect_s3_class(res, "crossstatr_correlation")
})
