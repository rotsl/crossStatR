# Usage: Rscript scripts/validate_caltech.R <extracted-data-directory>
# Source: Soedarmadji and Gregoire (2021), https://doi.org/10.25989/198p-1x22
# Test data is taken from https://doi.org/10.25989/198p-1x22.
library(crossStatR)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[1L])) {
  stop("An extracted Caltech dataset directory is required.")
}
files <- sort(list.files(args[1L], pattern = "^ana__[12]__.*_interlen[.]txt$", full.names = TRUE))
stopifnot(length(files) == 58L)
checks <- 0L
assert_equal <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance=tolerance)))
  checks <<- checks + 1L
}
results <- lapply(files, function(f) {
 d <- read.delim(f, check.names=FALSE)
 stopifnot(all(c("I(A)_illdiff", "t(s)_ill") %in% names(d)), is.numeric(d[["I(A)_illdiff"]]), is.numeric(d[["t(s)_ill"]]))
 input_rows <- nrow(d)
 keep <- is.finite(d[["t(s)_ill"]]) & is.finite(d[["I(A)_illdiff"]])
 d <- d[keep, , drop=FALSE]
 t <- d[["t(s)_ill"]]
 y <- d[["I(A)_illdiff"]]
 stopifnot(length(y)>=10L, all(diff(t)>0))
 # Microampere scaling avoids absolute tolerances larger than the signals.
 uA <- y*1e6
 description <- sci_describe(uA, units="microampere")$summary
 assert_equal(description$mean, mean(uA))
 assert_equal(description$sd, sd(uA))
 assert_equal(description$n, length(y))
 fit <- sci_regression(uA ~ t, data=data.frame(uA,t))
 reference <- lm(uA ~ t)
 assert_equal(fit$coefficients$estimate, coef(reference))
 noise <- noise_characterize(t, uA, method="difference")
 # The returned SD describes successive differences, without rescaling.
 assert_equal(noise$noise_sd, sd(diff(uA)))
 smooth <- adaptive_smoothing(t,uA,method="sgolay",window=7)
 stopifnot(length(smooth$smoothed)==length(uA),all(is.finite(smooth$smoothed)))
 checks <<- checks+1L
 data.frame(analysis=as.integer(sub("ana__([12])__.*","\\1",basename(f))),
  sample=as.integer(sub(".*Sample([0-9]+)_.*","\\1",basename(f))),
  input_rows=input_rows, excluded_rows=input_rows-length(y), observations=length(y), mean_photocurrent_uA=mean(uA),sd_photocurrent_uA=sd(uA),
  slope_uA_per_s=unname(coef(reference)[2L]),difference_noise_uA=noise$noise_sd,
  lag1_acf=unname(noise$acf[1L]))
})
r <- do.call(rbind,results)
write.csv(r,"inst/extdata/caltech-current-summary.csv",row.names=FALSE)
manifest <- data.frame(file=basename(files),md5=unname(tools::md5sum(files)))
write.csv(manifest,"data-raw/caltech-input-manifest.csv",row.names=FALSE)
write.csv(data.frame(checks=checks,curves=nrow(r),observations=sum(r$observations)),"docs/caltech-checks.csv",row.names=FALSE)
# Time complete analysis passes after one warm-up, excluding file input.
inputs <- lapply(files, function(f) {
 d <- read.delim(f,check.names=FALSE)
 d[is.finite(d[["t(s)_ill"]]) & is.finite(d[["I(A)_illdiff"]]), , drop=FALSE]
})
cases <- list(
 sci_describe=function() lapply(inputs,function(d) sci_describe(d[["I(A)_illdiff"]]*1e6)),
 sci_regression=function() lapply(inputs,function(d) {
  y <- d[["I(A)_illdiff"]]*1e6; t <- d[["t(s)_ill"]]; sci_regression(y~t,data=data.frame(y,t))
 }),
 noise_characterize=function() lapply(inputs,function(d) noise_characterize(d[["t(s)_ill"]],d[["I(A)_illdiff"]]*1e6,method="difference")),
 adaptive_smoothing=function() lapply(inputs,function(d) adaptive_smoothing(d[["t(s)_ill"]],d[["I(A)_illdiff"]]*1e6,window=7)))
timing <- do.call(rbind,lapply(names(cases),function(nm) {
 invisible(cases[[nm]]())
 times <- replicate(5,unname(system.time(cases[[nm]]())["elapsed"]))
 data.frame(method=nm,curves=nrow(r),median_s=median(times),min_s=min(times),max_s=max(times))
}))
write.csv(timing,"docs/caltech-benchmark-results.csv",row.names=FALSE)
cat(sprintf("Caltech numerical checks: %d passed across %d curves and %d observations.\n",checks,nrow(r),sum(r$observations)))

source("scripts/report_caltech.R")
