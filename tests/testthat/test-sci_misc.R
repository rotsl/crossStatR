
test_that("sci_normality basic", {
  set.seed(1)
  res <- sci_normality(rnorm(100))
  expect_s3_class(res, "crossstatr_normality")
  expect_equal(nrow(res$summary), 1L)
  expect_true("verdict" %in% colnames(res$summary))
})

test_that("sci_normality handles tiny n gracefully", {
  res <- sci_normality(1)
  expect_s3_class(res, "crossstatr_normality")
  expect_equal(res$summary$verdict, "too few observations")
})

test_that("sci_outliers basic", {
  set.seed(2)
  x <- c(rnorm(40), 10)
  res <- sci_outliers(x)
  expect_s3_class(res, "crossstatr_outliers")
  expect_true(any(res$flags$robust_z_flag, na.rm = TRUE))
})

test_that("sci_outliers with grouping", {
  set.seed(3)
  x <- c(rnorm(20), 8, rnorm(20), 9)
  g <- rep(c("A","B"), each = 21)
  res <- sci_outliers(x, group = g)
  expect_s3_class(res, "crossstatr_outliers")
  expect_equal(nrow(res$counts), 2L)
})

test_that("sci_bootstrap basic nonparametric", {
  set.seed(4)
  x <- rnorm(50)
  res <- sci_bootstrap(x, statistic = mean, n_replicates = 200, seed = 1)
  expect_s3_class(res, "crossstatr_bootstrap")
  expect_equal(length(res$replicates), 200L)
  expect_true(!is.null(res$ci))
})

test_that("sci_bootstrap BCa", {
  set.seed(5)
  x <- rnorm(30, mean = 5)
  res <- sci_bootstrap(x, statistic = mean, n_replicates = 200,
                       ci_type = "bca", seed = 1)
  expect_s3_class(res, "crossstatr_bootstrap")
})

test_that("sci_power two-sample t", {
  res <- sci_power(design = "two_sample_t", effect_size = 0.8, n1 = 26, n2 = 26)
  expect_s3_class(res, "crossstatr_power")
  expect_true(res$power > 0.5)
})

test_that("sci_power correlation", {
  res <- sci_power(design = "correlation", effect_size = 0.5, n = 30)
  expect_true(!is.na(res$power))
})

test_that("sci_effect_size basic d", {
  set.seed(6)
  x <- rnorm(20, 0)
  y <- rnorm(20, 0.5)
  res <- sci_effect_size(x, y, type = "cohens_d")
  expect_s3_class(res, "crossstatr_effect_size")
  expect_true(!is.na(res$estimate))
})

test_that("sci_effect_size hedges_g", {
  set.seed(7)
  x <- rnorm(20, 0)
  y <- rnorm(20, 0.5)
  res <- sci_effect_size(x, y, type = "hedges_g")
  expect_true(!is.na(res$estimate))
})
