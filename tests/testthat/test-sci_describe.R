
test_that("sci_describe handles a simple vector", {
  set.seed(1)
  x <- rnorm(100, mean = 5, sd = 2)
  res <- sci_describe(x)
  expect_s3_class(res, "crossstatr_describe")
  expect_equal(res$summary$n, 100)
  expect_equal(res$summary$mean, mean(x), tolerance = 1e-9)
  expect_equal(res$summary$sd, sd(x), tolerance = 1e-9)
  expect_equal(res$summary$median, median(x), tolerance = 1e-9)
  expect_equal(res$summary$iqr, IQR(x), tolerance = 1e-9)
})

test_that("sci_describe handles grouping", {
  set.seed(2)
  x <- c(rnorm(20, 0), rnorm(20, 5))
  g <- rep(c("A", "B"), each = 20)
  res <- sci_describe(x, group = g)
  expect_equal(nrow(res$summary), 2L)
  expect_true(res$summary$mean[2L] > res$summary$mean[1L])
})

test_that("sci_describe handles technical vs biological replicates", {
  set.seed(3)
  df <- data.frame(
    y = rnorm(60),
    grp = rep(c("A","B","C"), each = 20),
    bio = rep(1:6, each = 10),
    tech = rep(1:10, times = 6)
  )
  res <- sci_describe(df$y, group = df$grp, bio_rep = df$bio, tech_rep = df$tech)
  expect_s3_class(res, "crossstatr_describe")
  expect_equal(res$n_bio_rep, 6L)
  expect_true(res$n_obs == 60L)
})

test_that("sci_describe errors on non-numeric input", {
  expect_error(sci_describe("a"), "numeric")
})

test_that("sci_describe handles NA values", {
  x <- c(1, 2, NA, 4, 5)
  res <- sci_describe(x, na.rm = TRUE)
  expect_equal(res$summary$n, 4L)
})

test_that("sci_describe handles constant data", {
  res <- sci_describe(rep(5, 10))
  expect_equal(res$summary$sd, 0)
  expect_equal(res$summary$cv, 0)
})

test_that("sci_describe handles single-observation groups", {
  x <- c(1, 2)
  g <- c("A", "B")
  res <- sci_describe(x, group = g)
  expect_true(all(is.na(res$summary$sd)))
})

test_that("print method works", {
  res <- sci_describe(1:10)
  out <- capture.output(print(res))
  expect_true(any(grepl("crossStatR", out)))
})
