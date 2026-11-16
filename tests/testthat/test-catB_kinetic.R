
test_that("kinetic_model_compare basic", {
  set.seed(1)
  t <- seq(0, 50, length.out = 25)
  y <- 5 * exp(-0.1 * t) + rnorm(25, sd = 0.1)
  models <- list(
    first_order = function(p, t) p[1] * exp(-p[2] * t),
    biexp = function(p, t) p[1] * exp(-p[2] * t) + p[3] * exp(-p[4] * t)
  )
  starts <- list(first_order = c(A = 5, k = 0.1),
                 biexp = c(A = 4, k = 0.1, B = 1, k2 = 0.05))
  res <- kinetic_model_compare(t, y, models = models, start = starts)
  expect_s3_class(res, "crossstatr_kinetic_compare")
})

test_that("spectral_change_test T2", {
  set.seed(2)
  n1 <- 10; n2 <- 10; m <- 30
  y1 <- matrix(rnorm(n1 * m), nrow = n1)
  y2 <- matrix(rnorm(n2 * m, mean = 0.2), nrow = n2)
  grp <- factor(rep(c("A","B"), c(n1, n2)))
  spectra <- rbind(y1, y2)
  res <- spectral_change_test(spectra, group = grp, method = "t2")
  expect_s3_class(res, "crossstatr_spectral_change")
})

test_that("spectral_change_test pointwise", {
  set.seed(3)
  n1 <- 5; n2 <- 5; m <- 20
  y1 <- matrix(rnorm(n1 * m), nrow = n1)
  y2 <- matrix(rnorm(n2 * m, mean = 0.5), nrow = n2)
  grp <- factor(rep(c("A","B"), c(n1, n2)))
  spectra <- rbind(y1, y2)
  res <- spectral_change_test(spectra, group = grp, method = "pointwise_t")
  expect_s3_class(res, "crossstatr_spectral_change")
})

test_that("parameter_uncertainty delta method", {
  params <- c(a = 1, b = 2)
  cov <- matrix(c(0.01, 0, 0, 0.04), 2, 2)
  f <- function(p) p[1]^2 + p[2]
  res <- parameter_uncertainty(params, cov, f, method = "delta")
  expect_s3_class(res, "crossstatr_param_uncert")
  expect_true(!is.na(res$estimate))
})

test_that("parameter_uncertainty monte carlo", {
  params <- c(a = 1, b = 2)
  cov <- matrix(c(0.01, 0, 0, 0.04), 2, 2)
  f <- function(p) p[1]^2 + p[2]
  res <- parameter_uncertainty(params, cov, f, method = "monte_carlo",
                                n_sim = 100, seed = 1)
  expect_s3_class(res, "crossstatr_param_uncert")
  expect_equal(length(res$sim_samples), 100L)
})

test_that("timepoint_optimize basic", {
  model <- function(p, t) p[1] * exp(-p[2] * t)
  params <- c(1, 0.1)
  res <- timepoint_optimize(model, params, t_range = c(0, 50),
                            n_candidates = 5, criterion = "d_optimal")
  expect_s3_class(res, "crossstatr_timepoint")
})

test_that("timepoint_optimize a-optimal", {
  model <- function(p, t) p[1] * exp(-p[2] * t)
  params <- c(1, 0.1)
  res <- timepoint_optimize(model, params, t_range = c(0, 50),
                            n_candidates = 5, criterion = "a_optimal")
  expect_s3_class(res, "crossstatr_timepoint")
})
