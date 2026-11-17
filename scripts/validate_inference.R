# Numerical and Monte Carlo checks of corrected inference paths.
# Run after scripts/validate.R, with the current package installed.
library(crossStatR)
rows <- list()
add <- function(method, metric, v, target) {
  stopifnot(all(is.finite(v)))
  rows[[length(rows)+1L]] <<- data.frame(method, metric, estimate=mean(v),
                                       mcse=sd(v)/sqrt(length(v)), target, n=length(v))
}
set.seed(170101)
median_cover <- replicate(5000, {
  x <- rnorm(30, mean=10)
  ci <- sci_describe(x)$summary
  ci$median_ci_low <= 10 && ci$median_ci_high >= 10
})
add("sci_describe", "median 95% interval coverage, normal mean=10 n=30", median_cover, .95)
set.seed(170102)
kurt <- replicate(2000, sci_normality(rnorm(100))$summary$kurtosis)
add("sci_normality", "mean excess kurtosis, normal n=100", kurt, 0)
set.seed(170103)
bca <- percentile <- logical(1000)
for (i in seq_len(1000)) {
  x <- rlnorm(50)
  r <- sci_bootstrap(x, mean, n_replicates=999, ci_type="bca", seed=i)
  truth <- exp(.5)
  bca[i] <- r$ci[1] <= truth && truth <= r$ci[2]
  q <- quantile(r$replicates, c(.025,.975))
  percentile[i] <- q[1] <= truth && truth <= q[2]
}
add("sci_bootstrap", "BCa 95% mean interval coverage, lognormal n=50", bca, .95)
add("sci_bootstrap", "percentile 95% mean interval coverage, lognormal n=50", percentile, .95)
set.seed(170104)
drift_cover <- replicate(500, {
  time <- seq(0, 10, length.out=30)
  reference <- data.frame(time=seq(0,10,length.out=12))
  reference$y <- .5+.1*reference$time+rnorm(12,sd=.2)
  y <- 5+.5+.1*time+rnorm(30,sd=.2)
  r <- drift_correct_stat(y,time,reference,drift_model="linear")
  truth <- 5+.5+.1*median(time)
  mean(r$correction_ci[,1] <= truth & truth <= r$correction_ci[,2])
})
add("drift_correct_stat", "mean pointwise 95% interval coverage with independent standards", drift_cover, .95)
# Built-in PlantGrowth data: no copied dataset is needed.
d <- datasets::PlantGrowth
classic <- sci_anova(weight~group,d)
welch <- sci_anova(weight~group,d,welch=TRUE,posthoc="none")
stopifnot(isTRUE(all.equal(classic$anova_table$`Pr(>F)`[1],summary(aov(weight~group,d))[[1]]$`Pr(>F)`[1])),
          isTRUE(all.equal(unname(welch$anova_table$`Pr(>F)`),unname(oneway.test(weight~group,d)$p.value))))
r <- do.call(rbind,rows)
write.csv(r,"docs/inference-validation-results.csv",row.names=FALSE)
file <- "docs/VALIDATION_REPORT.md"
old <- readLines(file)
heading <- "## Inference checks"
pos <- match(heading,old)
if (!is.na(pos)) old <- head(old,pos-1L)
lines <- c(heading,"",paste0("Local run: ",Sys.Date(),". Source: `scripts/validate_inference.R`."),"",
"The median, kurtosis, BCa and drift studies use initial seeds 170101, 170102, 170103 and 170104. BCa resampling uses seed i for dataset i, with 999 resamples. All coverage experiments have nominal level 95%. MCSE is computed from independent dataset-level metrics, including pointwise coverage averaged within each drift dataset.","",
"| Method | Metric | Estimate | MCSE | Reference | Datasets |","|---|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(r))) lines <- c(lines,sprintf("| %s | %s | %.5f | %.5f | %.2f | %d |",r$method[i],r$metric[i],r$estimate[i],r$mcse[i],r$target[i],r$n[i]))
means <- tapply(d$weight,d$group,mean)
lines <- c(lines,"","### PlantGrowth numerical example","",
 paste0("The built-in `datasets::PlantGrowth` example has group means ",paste(sprintf("%s = %.3f",names(means),means),collapse=", "),". The classical and Welch ANOVA p-values are ",sprintf("%.6f and %.6f",classic$anova_table$`Pr(>F)`[1],welch$anova_table$`Pr(>F)`[1]),"; both matched base R. This numerical comparison does not establish interval coverage or a causal treatment conclusion. Dataset provenance is recorded in `help(\"PlantGrowth\", package = \"datasets\")`."),"",
"BCa remains below nominal coverage in this skewed small-sample setting. Correcting the acceleration sign does not guarantee adequate coverage. Drift intervals combine reference-fit uncertainty and a shared noise estimate, assuming independent standards and samples with equal noise SD. Their coverage does not apply automatically to the running-median fallback or unequal reference/sample noise. Shapley precision intervals condition on supplied covariance and do not estimate uncertainty in experimental variance components.")
writeLines(c(old,"",lines),file)
