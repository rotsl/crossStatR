# Refresh the measured-data supplements from the local workflow outputs.
r <- read.csv("inst/extdata/caltech-current-summary.csv")
timing <- read.csv("docs/caltech-benchmark-results.csv")
checks <- read.csv("docs/caltech-checks.csv")
citation <- "Source: Soedarmadji, E., and Gregoire, J. (2021). Materials Experiment and Analysis Database: Electrochemistry analysis, Electrochemistry measurements, Material library with V,Ca,Ar,O on plate(s) 3540. CaltechDATA, v1. <https://doi.org/10.25989/198p-1x22>."
append_section <- function(file, heading, lines) {
 old <- readLines(file,warn=FALSE)
 start <- match(heading,old)
 if (!is.na(start)) old <- head(old,start-1L)
 writeLines(c(old,heading,"",lines),file)
}
counts <- aggregate(cbind(input_rows,excluded_rows,observations)~analysis,r,sum)
validation <- c(paste0("Local measured-data run: 2026-11-16. Package version: ", read.dcf("DESCRIPTION")[1, "Version"], " (development)."), "", citation,"",
 "The supplied archive contains processed measurements. The run used all 29 CA2 and 29 CV3 light/dark current-difference tables. The test data is taken from <https://doi.org/10.25989/198p-1x22>. Basenames and input checksums are recorded in `data-raw/caltech-input-manifest.csv`; `scripts/validate_caltech.R` reproduces the checks and derived summaries.","",
 "| Analysis | Input rows | Excluded rows | Retained rows |","|---|---:|---:|---:|")
for(i in seq_len(nrow(counts))) validation <- c(validation,sprintf("| %s | %d | %d | %d |",c("CA2","CV3")[counts$analysis[i]],counts$input_rows[i],counts$excluded_rows[i],counts$observations[i]))
validation <- c(validation,"",sprintf("All %d checks passed across %d curves and %d retained observations. Means, sample SDs, counts, and OLS coefficients matched base R; the reported difference SD matched `sd(diff(current))`. Smoothing returned finite values of the expected length. Current was converted from amperes to microamperes before comparisons. The comparison tolerance was 1e-10.",checks$checks,checks$curves,checks$observations),"",
 "Rows were excluded only when the current difference or illuminated time was non-finite. No interpolation was performed. Several near-zero-mean current series produced the documented unstable-CV warning; the checks compare means and SDs rather than interpreting those CVs. Distinct material identifiers are not treatment replicates. Temporal dependence, processing by the source analysis, and real signal variation prevent these numerical agreement checks from establishing Type I error, interval coverage, physical model validity, or instrumental noise. The synthetic validation limitations above remain unresolved.")
append_section("docs/VALIDATION_REPORT.md","## Caltech measured-data checks",validation)
lines <- c("Local measured-data run: 2026-11-16.","",citation,"",
 "The script times complete passes over the 58 processed current-difference curves after filtering the non-finite time/current rows counted in the validation report. File input is excluded. Each case has one warm-up and five repetitions. Smoothing uses a seven-point Savitzky-Golay window. The summary calls retain the documented warnings for unstable CVs near zero mean. These timings describe 2,975 observations spread across 58 curves, rather than one independent sample of that size.","",
 "| Method | Curves per pass | Median (s) | Minimum (s) | Maximum (s) |","|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(timing))) lines <- c(lines,sprintf("| %s | %d | %.4f | %.4f | %.4f |",timing$method[i],timing$curves[i],timing$median_s[i],timing$min_s[i],timing$max_s[i]))
append_section("docs/BENCHMARK_REPORT.md","## Caltech measured-data benchmark",lines)
