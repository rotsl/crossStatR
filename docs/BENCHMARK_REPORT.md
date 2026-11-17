# Benchmark report

Local run: 2026-11-17. Package version: 0.99.3 (development).
Environment: R version 4.6.2 (2026-10-29) on aarch64-apple-darwin23 .

The script `scripts/benchmark.R` measures elapsed seconds after one warm-up call, with five repetitions per case. Synthetic inputs use seed 20261116. Curve tests use eight replicates per group and 99 permutations; curve similarity uses 50 bootstrap draws. Grid size is the number of points per curve. Input construction inside the signal-region case is included in its timing.

| Method | Points | Median (s) | Minimum (s) | Maximum (s) |
|---|---:|---:|---:|---:|
| sci_describe | 100 | 0.0000 | 0.0000 | 0.0010 |
| sci_compare | 100 | 0.0000 | 0.0000 | 0.0010 |
| adaptive_baseline | 100 | 0.0030 | 0.0020 | 0.0030 |
| adaptive_smoothing | 100 | 0.0030 | 0.0030 | 0.0040 |
| noise_characterize | 100 | 0.0020 | 0.0020 | 0.0020 |
| curve_similarity | 100 | 0.0010 | 0.0000 | 0.0010 |
| peak_evidence | 100 | 0.0030 | 0.0030 | 0.0030 |
| signal_difference_regions | 100 | 0.0440 | 0.0430 | 0.0460 |
| sci_describe | 500 | 0.0010 | 0.0010 | 0.0010 |
| sci_compare | 500 | 0.0000 | 0.0000 | 0.0010 |
| adaptive_baseline | 500 | 0.0050 | 0.0040 | 0.0050 |
| adaptive_smoothing | 500 | 0.0640 | 0.0610 | 0.0650 |
| noise_characterize | 500 | 0.0020 | 0.0020 | 0.0020 |
| curve_similarity | 500 | 0.0010 | 0.0000 | 0.0010 |
| peak_evidence | 500 | 0.0040 | 0.0040 | 0.0050 |
| signal_difference_regions | 500 | 0.1970 | 0.1960 | 0.2060 |
| sci_describe | 1000 | 0.0000 | 0.0000 | 0.0010 |
| sci_compare | 1000 | 0.0000 | 0.0000 | 0.0010 |
| adaptive_baseline | 1000 | 0.0080 | 0.0080 | 0.0080 |
| adaptive_smoothing | 1000 | 0.2180 | 0.2130 | 0.2510 |
| noise_characterize | 1000 | 0.0020 | 0.0010 | 0.0020 |
| curve_similarity | 1000 | 0.0010 | 0.0010 | 0.0010 |
| peak_evidence | 1000 | 0.0060 | 0.0050 | 0.0060 |
| signal_difference_regions | 1000 | 0.3890 | 0.3850 | 0.3900 |

These measurements cover the listed configurations only. Times near the timer resolution are not suitable for speed comparisons. Runtime depends on input size, resampling counts, convergence, and system load. Memory, parallel execution, and GPU execution were not measured. No scaling exponent or memory bound is inferred.

## Caltech measured-data benchmark

Local measured-data run: 2026-11-17.

Source: Soedarmadji, E., and Gregoire, J. (2021). Materials Experiment and Analysis Database: Electrochemistry analysis, Electrochemistry measurements, Material library with V,Ca,Ar,O on plate(s) 3540. CaltechDATA, v1. <https://doi.org/10.25989/198p-1x22>.

The script times complete passes over the 58 processed current-difference curves after filtering the non-finite time/current rows counted in the validation report. File input is excluded. Each case has one warm-up and five repetitions. Smoothing uses a seven-point Savitzky-Golay window. The summary calls retain the documented warnings for unstable CVs near zero mean. These timings describe 2,975 observations spread across 58 curves, rather than one independent sample of that size.

| Method | Curves per pass | Median (s) | Minimum (s) | Maximum (s) |
|---|---:|---:|---:|---:|
| sci_describe | 58 | 0.0240 | 0.0230 | 0.0250 |
| sci_regression | 58 | 0.0350 | 0.0340 | 0.0370 |
| noise_characterize | 58 | 0.0600 | 0.0570 | 0.0640 |
| adaptive_smoothing | 58 | 0.0840 | 0.0810 | 0.0970 |
