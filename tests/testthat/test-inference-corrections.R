test_that("normality moments agree with finite-sample moment formulas", {
  x <- c(1, 2, 3, 4, 9, 15)
  n <- length(x)
  m2 <- mean((x - mean(x))^2)
  r <- sci_normality(x)$summary
  expect_equal(r$skewness, sqrt(n*(n-1))/(n-2)*mean((x-mean(x))^3)/m2^1.5)
  expect_equal(r$kurtosis, (n-1)/((n-2)*(n-3))*((n+1)*mean((x-mean(x))^4)/m2^2-3*(n-1)))
})

test_that("robust regression inference and multi-predictor variance test are finite", {
  set.seed(901)
  d <- data.frame(x = seq(-2, 2, length.out = 60), z = rnorm(60))
  d$y <- 2 + d$x - 0.5*d$z + rnorm(60, sd = 0.2)
  r <- sci_regression(y ~ x + z, d, robust = TRUE)
  expect_true(all(is.finite(r$coefficients$p_value)))
  expect_true(all(is.finite(r$coefficients$ci_low)))
  p <- crossStatR:::.bptest(y ~ x + z, d)
  f <- lm(y ~ x + z, d)
  aux <- lm(resid(f)^2 ~ x + z, d)
  expect_equal(p, pchisq(nrow(d)*summary(aux)$r.squared, 2, lower.tail = FALSE))
})

test_that("BCa intervals agree with an independent mean-acceleration calculation", {
  x <- c(1, 2, 3, 4, 6, 12, 30)
  r <- sci_bootstrap(x, mean, ci_type = "bca", n_replicates = 999, seed = 902)
  a <- sum((x-mean(x))^3)/(6*sum((x-mean(x))^2)^1.5)
  z0 <- qnorm(mean(r$replicates < mean(x)))
  z <- qnorm(c(0.025, 0.975))
  expected <- quantile(r$replicates, pnorm(z0+(z0+z)/(1-a*(z0+z))))
  expect_equal(unname(r$ci), unname(expected))
})

test_that("paired Hedges correction and its interval scale together", {
  x <- c(1, 3, 6, 9, 12, 14)
  y <- c(2, 3, 4, 8, 9, 10)
  r <- sci_compare(x, y, paired = TRUE, test = "student")$effect_size
  J <- 1-3/(4*(length(x)-1)-1)
  expect_equal(r$estimate[2], r$estimate[1]*J)
  expect_equal(r$se[2], r$se[1]*J)
  expect_equal(r$ci_high[2]-r$ci_low[2], J*(r$ci_high[1]-r$ci_low[1]))
  expect_equal(crossStatR:::.ci_median(1:20), c(6, 15))
})

test_that("Welch degrees of freedom and factorial term effects are retained", {
  r <- sci_anova(weight ~ group, datasets::PlantGrowth, welch = TRUE, posthoc = "none")
  expected <- oneway.test(weight ~ group, datasets::PlantGrowth)
  expect_equal(unname(r$anova_table$denominator_df), unname(expected$parameter[2]))
  expect_true(all(is.finite(r$effect_sizes$estimate)))
  d <- expand.grid(a = factor(1:2), b = factor(1:3), rep = 1:5)
  set.seed(903); d$y <- as.numeric(d$a)+as.numeric(d$b)+rnorm(nrow(d))
  r <- sci_anova(y ~ a*b, d, posthoc = "none")
  expect_equal(length(unique(trimws(r$effect_sizes$term))), 3L)
  expect_error(sci_anova(y ~ a*b, d, welch = TRUE), "one-way")
})

test_that("detection uses a one-sided threshold and supplies LOD uncertainty", {
  set.seed(904); x <- 1:30; y <- 1+2*x+rnorm(30)
  f <- lm(y~x)
  r <- detection_reliability(f, n_boot = 99, seed = 905)
  expect_equal(unname(r$decision_threshold), unname(coef(f)[1]+qnorm(0.95)*summary(f)$sigma))
  expect_equal(unname(r$fpr), 0.05)
  expect_true(all(is.finite(r$lod_ci)))
  expect_lte(r$lod_ci[1], r$lod_ci[2])
  blank <- rnorm(20, 1, 1)
  r <- detection_reliability(f, blank_data = blank, n_boot = 99, seed = 905)
  expect_equal(r$decision_threshold, mean(blank)+qt(.95,19)*sd(blank)*sqrt(1+1/20))
  expect_error(detection_reliability(lm(-y~x)), "positive-slope")
})

test_that("delta propagation uses normal or supplied residual degrees of freedom", {
  f <- function(p) p[1]+2*p[2]
  V <- diag(c(0.25, 1))
  r <- parameter_uncertainty(c(1,2), V, f)
  expect_equal(r$ci_high, 5+qnorm(.975)*sqrt(4.25))
  r <- parameter_uncertainty(c(1,2), V, f, df = 8)
  expect_equal(r$ci_high, 5+qt(.975,8)*sqrt(4.25))
  expect_error(parameter_uncertainty(c(1,2), V, f, df = 0), "df")
})

test_that("uncertainty-origin intervals measure conditional simulation precision", {
  r <- uncertainty_origin(c(A=1, B=2), correlations=list(A__B=.3), n_mc=1000, seed=906)
  expect_equal(sum(r$shapley_attribution), r$total_variance_correlated)
  expect_equal(r$total_variance_ci, 999*r$total_variance_mc/qchisq(c(.975,.025),999))
  expect_equal(r$total_ci, r$total_prediction_interval)
  expect_true(all(is.finite(r$shapley_ci)))
  expect_equal(unname(r$shapley_mcse^2), (diag(r$covariance)*sum(r$covariance)+rowSums(r$covariance)^2)/999)
  expect_error(uncertainty_origin(c(A=1,B=2), n_mc=1), "at least 2")
})

test_that("D-optimal rank completion selects a full-rank three-parameter design", {
  f <- function(p,t) p[1]+p[2]*t+p[3]*t^2
  expect_warning(r <- timepoint_optimize(f, c(1,1,1), t_range=c(0,1),
                                        n_candidates=21, n_extra=3), "rank-deficient")
  expect_length(r$optimal_times, 3)
  expect_equal(qr(r$sensitivity_matrix)$rank, 3L)
  expect_gt(r$criterion_value, 0)
})

test_that("zero MAD preserves flags for deviations and handles constant data", {
  r <- sci_outliers(c(rep(1,9),10))$flags
  expect_true(r$robust_z_flag[10])
  expect_equal(r$robust_z[1:9], rep(0,9))
  r <- sci_outliers(rep(1,10))$flags
  expect_false(any(r$robust_z_flag | r$z_flag))
})

test_that("random-slope variance intervals match named lme4 profile parameters", {
  r <- replicate_decompose(Reaction ~ Days + (Days|Subject), lme4::sleepstudy)
  v <- r$variance_components
  expect_equal(v$component, c("Subject","Subject:Days","Residual"))
  expect_true(all(v$ci_method == "profile"))
  expect_true(all(v$ci_low <= v$variance & v$variance <= v$ci_high))
  expected <- as.data.frame(lme4::VarCorr(r$model))
  expect_equal(v$variance, expected$vcov[is.na(expected$var2)])
})

test_that("fragility rejects invalid original analysis results", {
  for (value in list(NA_real_, Inf, c(1,2), "1", list(1))) {
    expect_error(experiment_fragility(function(opts) value, data.frame(y=1:5),
                                      decision_threshold=0), "single finite numeric")
  }
})

test_that("weighted variance diagnostics evaluate supplied weights", {
  set.seed(909)
  d <- data.frame(x=1:30, y=1:30+rnorm(30))
  expect_true(is.finite(crossStatR:::.bptest(y~x,d,weights=seq(1,2,length.out=30))))
})

test_that("paired Hedges correction and its standard error use pair df", {
  x <- c(2,4,5,9,12,18)
  y <- c(1,2,4,5,7,11)
  r <- sci_compare(x,y,paired=TRUE)
  J <- 1-3/(4*(length(x)-1)-1)
  expect_equal(r$effect_size$estimate[2], J*mean(x-y)/sd(x-y))
  expect_equal(r$effect_size$se[2], J*r$effect_size$se[1])
  expect_equal(r$effect_size$ci_low[2], J*r$effect_size$ci_low[1])
})

test_that("periodogram uses actual sample count for non-power-of-two inputs", {
  set.seed(910)
  r <- noise_characterize(1:30, rnorm(30), method="difference")
  e <- r$residuals-mean(r$residuals)
  nfft <- nextn(length(e),2)
  p <- Mod(fft(c(e,rep(0,nfft-length(e)))))^2/length(e)
  expect_equal(r$psd$power, p[seq_len(nrow(r$psd))])
})

test_that("kinetic information criteria count residual variance", {
  set.seed(911)
  t <- seq(0,10,length.out=25)
  y <- 5*exp(-.2*t)+rnorm(25,sd=.1)
  r <- kinetic_model_compare(t,y,models=list(decay=function(p,t) p[1]*exp(-p[2]*t)),
                              start=list(decay=c(A=5,k=.2)))$comparison
  expect_equal(r$AIC, 25*log(r$RSS/25)+2*3)
  expect_equal(r$BIC, 25*log(r$RSS/25)+3*log(25))
})

test_that("drift intervals combine correction uncertainty and reference noise", {
  set.seed(912)
  t <- seq(0,10,length.out=30)
  reference <- data.frame(time=t[seq(1,30,by=3)])
  reference$y <- .2*reference$time+rnorm(nrow(reference),sd=.1)
  y <- 5+.2*t+rnorm(30,sd=.1)
  r <- drift_correct_stat(y,t,reference=reference,drift_model="linear")
  half <- qt(.975,nrow(reference)-2)*sqrt(r$se_corr^2+r$sigma_y^2)
  expect_equal(unname(r$correction_ci[,1]), r$y_corrected-half)
  expect_equal(unname(r$correction_ci[,2]), r$y_corrected+half)
})

test_that("kinetic selection reports all-candidate failures separately", {
  f <- function(p,t) p[1]^2*exp(-t)
  r <- kinetic_discrimination(list(a=f,b=f),list(c(a=0),c(a=0)),
                             true_model=1L,true_params=1,times=1:10,
                             sigma=.1,n_sim=5,seed=913)
  expect_equal(r$n_no_convergence,5L)
  expect_equal(sum(r$discrimination_probabilities),0)
  expect_true(all(is.na(r$discrimination_probabilities_clean)))
  expect_false(any(r$converged))
})

test_that("fixed batch percentages use the returned variance components", {
  set.seed(914)
  d <- expand.grid(batch=factor(1:4),treat=factor(c("A","B")),replicate=1:8)
  d$y <- as.numeric(d$batch)+.5*as.numeric(d$treat)+rnorm(nrow(d))
  r <- batch_effect_quantify(y~treat,d,batch="batch",type="fixed")
  expect_equal(r$batch_pct_variance,r$tidy$pct[1])
})

test_that("median interval order-statistic centre uses the sample midpoint", {
  expect_equal(crossStatR:::.ci_median(1:30), c(10,21))
})
