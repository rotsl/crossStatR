# Reproducible elapsed-time measurements; run from the package root.
library(crossStatR)
set.seed(20261116)
rows <- list()
for (n in c(100L, 500L, 1000L)) {
  x <- seq(0, 20, length.out = n)
  y <- sin(x) + rnorm(n, sd = 0.1)
  cases <- list(
    sci_describe = function() sci_describe(y),
    sci_compare = function() sci_compare(y, y + 0.1),
    adaptive_baseline = function() adaptive_baseline(x, y, method = "asls"),
    adaptive_smoothing = function() adaptive_smoothing(x, y, window = 7),
    noise_characterize = function() noise_characterize(x, y),
    curve_similarity = function() curve_similarity(x, y, x, y + 0.1, n_boot = 50),
    peak_evidence = function() peak_evidence(x, y, peak_index = as.integer(n/2), n_perm = 99, seed = 1),
    signal_difference_regions = function() {
      a <- matrix(rnorm(8*n), 8, n)
      signal_difference_regions(x, a, a + matrix(rnorm(8*n), 8, n), n_perm = 99, seed = 1)
    })
  for (nm in names(cases)) {
    invisible(cases[[nm]]())
    times <- replicate(5, unname(system.time(cases[[nm]]())["elapsed"]))
    rows[[length(rows)+1L]] <- data.frame(method=nm, n=n, median=median(times), min=min(times), max=max(times))
  }
}
r <- do.call(rbind, rows)
write.csv(r, "docs/benchmark-results.csv", row.names=FALSE)
lines <- c("# Benchmark report", "", "Local run: 2026-11-16. Package version: 0.99.1 (development).",
 paste("Environment:", R.version.string, "on", R.version$platform, "."),
 "", "The script `scripts/benchmark.R` measures elapsed seconds after one warm-up call, with five repetitions per case. Synthetic inputs use seed 20261116. Curve tests use eight replicates per group and 99 permutations; curve similarity uses 50 bootstrap draws. Grid size is the number of points per curve. Input construction inside the signal-region case is included in its timing.",
 "", "| Method | Points | Median (s) | Minimum (s) | Maximum (s) |", "|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(r))) lines <- c(lines, sprintf("| %s | %d | %.4f | %.4f | %.4f |", r$method[i],r$n[i],r$median[i],r$min[i],r$max[i]))
writeLines(c(lines,"", "These measurements cover the listed configurations only. Times near the timer resolution are not suitable for speed comparisons. Runtime depends on input size, resampling counts, convergence, and system load. Memory, parallel execution, and GPU execution were not measured. No scaling exponent or memory bound is inferred."), "docs/BENCHMARK_REPORT.md")
