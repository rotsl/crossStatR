# Monte Carlo diagnostics for eight advanced methods; run from the package root.
library(crossStatR)
N <- 100L
rows <- list()
add <- function(method, metric, values, target, probability=FALSE) {
  stopifnot(length(values)==N, all(is.finite(values)))
  estimate <- mean(values)
  se <- sd(values)/sqrt(N)
  rows[[length(rows)+1L]] <<- data.frame(method,metric,estimate,mcse=se,target,n=N)
}
x <- seq(0, 1, length.out=30)
null <- power <- peak_null <- peak_power <- numeric(N)
transfer_bias <- transfer_coverage <- drift_bias <- drift_coverage <- numeric(N)
lod_bias <- attribution_error <- fragility <- numeric(N)
for(i in seq_len(N)) {
 set.seed(20261116+i)
 a <- matrix(rnorm(8*30,sd=0.2),8,30)
 b <- matrix(rnorm(8*30,sd=0.2),8,30)
 null[i] <- nrow(signal_difference_regions(x,a,b,n_perm=99,seed=i)$regions)>0
 power[i] <- nrow(signal_difference_regions(x,a,b+0.4,n_perm=99,seed=i)$regions)>0
 px <- seq(0,20,length.out=101)
 py <- rnorm(101,sd=0.1)
 peak_null[i] <- peak_evidence(px,py,peak_index=51L,noise_sd=0.1,n_perm=99,seed=i)$peaks$signif
 peak_power[i] <- peak_evidence(px,py+exp(-(px-10)^2/2),peak_index=51L,noise_sd=0.1,n_perm=99,seed=i)$peaks$signif
 cx <- seq(0.1,10,length.out=30)
 sy <- 2+3*cx+rnorm(30,sd=0.2)
 ty <- 2.5+2.7*cx+rnorm(30,sd=0.2)
 cal <- lm(sy~cx)
 tr <- calibration_transfer_risk(cal,data.frame(x=cx,y=sy),data.frame(x=cx,y=ty),transfer_model="slope_offset",n_boot=99,seed=i)$predictions
 transfer_bias[i] <- mean(tr$predicted_x-cx)
 transfer_coverage[i] <- mean(tr$ci_low<=cx & tr$ci_high>=cx)
 tt <- seq(0,10,length.out=50)
 dy <- 0.3*tt+rnorm(50,sd=0.2)
 dr <- drift_uncertainty(dy,time=tt,drift_model="linear",n_boot=99,seed=i)
 truth <- var(0.3*tt)/(var(0.3*tt)+0.2^2)
 drift_bias[i] <- dr$drift_fraction-truth
 drift_coverage[i] <- dr$drift_fraction_ci[1]<=truth && dr$drift_fraction_ci[2]>=truth
 dl <- detection_reliability(cal,n_boot=99,seed=i)
 lod_bias[i] <- dl$lod-(2*qnorm(0.95)*0.2/3)
 u <- uncertainty_origin(c(A=1,B=2),correlations=list(A__B=0.3),n_mc=1000,seed=i)
 attribution_error[i] <- u$total_variance_mc-(3+2*0.3*sqrt(2))
 df <- data.frame(y=c(rnorm(10,0,0.3),rnorm(10,1,0.3)),grp=rep(c("A","B"),each=10),rep=rep(1:10,2))
 af <- function(opts) {
   d <- opts$data
   if(!is.null(opts$rep_lo)) d <- d[d$rep!=opts$rep_lo,]
   unname(diff(tapply(d$y,d$grp,mean)))
 }
 fragility[i] <- experiment_fragility(af,df,replicates="rep",decision_threshold=0,n_boot=30,seed=i)$overall_fragility
}
add("signal_difference_regions","null rejection rate",null,"0.05 or lower")
add("signal_difference_regions","power for uniform shift 0.4",power,"high")
add("peak_evidence","null rejection at prespecified centre",peak_null,"0.05 or lower")
add("peak_evidence","detection at amplitude 1, noise SD 0.1",peak_power,"high")
add("calibration_transfer_risk","mean concentration bias",transfer_bias,"0")
add("calibration_transfer_risk","mean pointwise 95% interval coverage",transfer_coverage,"0.95")
add("drift_uncertainty","drift fraction bias",drift_bias,"0")
add("drift_uncertainty","95% interval coverage",drift_coverage,"0.95")
add("detection_reliability","LOD bias",lod_bias,"0")
add("uncertainty_origin","MC total variance minus analytic variance",attribution_error,"0")
add("experiment_fragility","strong effect overall fragility",fragility,"low")
m1 <- function(p,t) p[1]*exp(-p[2]*t)
m2 <- function(p,t) p[1]/(1+p[2]*t)
k <- kinetic_discrimination(list(exponential=m1,hyperbolic=m2),list(c(A=1,k=0.2),c(A=1,k=0.2)),true_model=1L,true_params=c(1,0.2),times=seq(0,20,length.out=20),sigma=0.05,n_sim=100,seed=20261116)
saveRDS(k,"/tmp/crossstatr-kinetic.rds")
r <- do.call(rbind,rows)
write.csv(r,"docs/validation-results.csv",row.names=FALSE)
lines <- c("# Validation report","","Local run: 2026-11-16. Package version: 0.99.1 (development).",paste("Environment:",R.version.string,"on",R.version$platform,"."),"",
"The script `scripts/validate.R` runs 100 independent synthetic datasets per listed metric. Dataset seeds are 20261117 through 20261216; resampling seeds are the dataset indices. Intervals use 99 bootstrap draws and permutation tests use 99 permutations. MCSE is the sample standard deviation of the per-dataset metric divided by the square root of 100. Coverage averaged across concentrations is summarised at the dataset level. Zero observed events do not establish a zero underlying event probability (a rough 95% upper bound is 3/100).",
"", "| Method | Metric | Estimate | MCSE | Reference | Datasets |","|---|---|---:|---:|---|---:|")
for(i in seq_len(nrow(r))) lines <- c(lines,sprintf("| %s | %s | %.5f | %.5f | %s | %d |",r$method[i],r$metric[i],r$estimate[i],r$mcse[i],r$target[i],r$n[i]))
lines <- c(lines,"","## Kinetic discrimination","","A separate 100-simulation experiment compares exponential and hyperbolic decay with 20 observations, amplitude 1, rate 0.2, and noise SD 0.05. Seed: 20261116.",capture.output(print(k)),"",
"## Interpretation","","These are limited simulation diagnostics, not evidence of general calibration. The transfer intervals omit new-response noise and can under-cover true concentrations. Drift coverage depends on the linear drift model. The peak procedure retains a smooth signal estimate in its permutation construction, so detection and null behaviour require independent methodological review. The uncertainty check compares Monte Carlo propagation with the supplied analytic covariance; it does not assess confidence intervals for estimated stage variances. Fragility measures sensitivity to the implemented perturbations, rather than a probability that a conclusion is correct. LOD bias is compared with the known positive-slope Gaussian calibration model. Larger simulations and independent comparisons remain necessary.")
lines <- c(lines,"", "The [methods guide](METHODS.md) distinguishes implemented variants from the published background methods. The [reference list](REFERENCES.md) supplies source identifiers. Bibliographic verification does not resolve the calibration limitations reported here.")
writeLines(lines,"docs/VALIDATION_REPORT.md")
