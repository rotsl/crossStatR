
test_that("sci_compare welch gives sensible results", {
  set.seed(1)
  x <- rnorm(30, 10)
  y <- rnorm(30, 11)
  res <- sci_compare(x, y, test = "welch")
  expect_s3_class(res, "crossstatr_compare")
  expect_lt(res$p_value, 0.05)
  expect_true("Cohen_d" %in% res$effect_size$name)
})

test_that("sci_compare handles paired data", {
  set.seed(2)
  x <- rnorm(20, 10)
  y <- x + rnorm(20, 0.5)
  res <- sci_compare(x, y, paired = TRUE, test = "student")
  expect_true(!is.na(res$p_value))
  expect_equal(res$n1, 20L)
})

test_that("sci_compare permutation test", {
  set.seed(3)
  x <- rnorm(15, 0)
  y <- rnorm(15, 1)
  res <- sci_compare(x, y, test = "perm", n_perm = 100, seed = 1)
  expect_lt(res$p_value, 0.1)
})

test_that("sci_compare wilcoxon test", {
  set.seed(4)
  x <- rnorm(30, 0)
  y <- rnorm(30, 0.5)
  res <- sci_compare(x, y, test = "wilcoxon")
  expect_true(!is.na(res$p_value))
})

test_that("sci_compare errors on length mismatch in paired", {
  expect_error(sci_compare(1:5, 1:6, paired = TRUE), "equal-length")
})

test_that("sci_anova basic one-way", {
  set.seed(5)
  df <- data.frame(
    y = c(rnorm(20, 5), rnorm(20, 6), rnorm(20, 8)),
    g = factor(rep(c("A","B","C"), each = 20))
  )
  res <- sci_anova(y ~ g, data = df)
  expect_s3_class(res, "crossstatr_anova")
  expect_lt(res$anova_table$`Pr(>F)`[1L], 0.05)
  expect_true(!is.null(res$effect_sizes))
})

test_that("sci_anova welch alternative", {
  set.seed(6)
  df <- data.frame(
    y = c(rnorm(20, 5, sd = 1), rnorm(20, 6, sd = 2)),
    g = factor(rep(c("A","B"), each = 20))
  )
  res <- sci_anova(y ~ g, data = df, welch = TRUE)
  expect_true(res$welch)
})

test_that("sci_anova with tukey post-hoc", {
  set.seed(7)
  df <- data.frame(
    y = c(rnorm(15, 5), rnorm(15, 6), rnorm(15, 8)),
    g = factor(rep(c("A","B","C"), each = 15))
  )
  res <- sci_anova(y ~ g, data = df, posthoc = "tukey")
  expect_true(!is.null(res$posthoc))
})
