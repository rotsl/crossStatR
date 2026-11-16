caltech_summary <- function() {
  read.csv(system.file("extdata", "caltech-current-summary.csv", package="crossStatR"))
}

test_that("Caltech diagnostics retain material and analysis identities", {
  d <- caltech_summary()
  expect_equal(nrow(d),58L)
  expect_equal(as.integer(table(d$analysis)),c(29L,29L))
  expect_equal(length(unique(d$sample)),29L)
  expect_equal(anyDuplicated(d[c("analysis","sample")]),0L)
  expect_equal(d$observations+d$excluded_rows,d$input_rows)
  expect_true(all(is.finite(as.matrix(d))))
})

test_that("descriptive statistics agree with base R on measured current summaries", {
  d <- subset(caltech_summary(), analysis==1)
  r <- sci_describe(d$mean_photocurrent_uA,units="microampere")
  expect_equal(r$summary$mean, mean(d$mean_photocurrent_uA),tolerance=1e-10)
  expect_equal(r$summary$sd,sd(d$mean_photocurrent_uA),tolerance=1e-10)
  expect_equal(r$summary$median,median(d$mean_photocurrent_uA),tolerance=1e-10)
})

test_that("regression agrees with lm on measured current summaries", {
  d <- subset(caltech_summary(),analysis==1)
  r <- sci_regression(mean_photocurrent_uA~sample,data=d)
  reference <- lm(mean_photocurrent_uA~sample,data=d)
  expect_equal(r$coefficients$estimate,unname(coef(reference)),tolerance=1e-10)
})
