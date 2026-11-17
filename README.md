# crossStatR

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)

`crossStatR` is an R package of statistics for experimental science. It supports biological assays, electrochemical measurements, optical and spectroscopic signals, and physics experiments. The package combines established methods (descriptive statistics, t-tests, regression, ANOVA, power analysis) with experimental diagnostics (region-level signal-difference testing, peak-evidence assessment, calibration-transfer risk, kinetic discrimination, experimental fragility, drift-attributed uncertainty, analytical detection-reliability curves, and stage-wise uncertainty-origin attribution).

This is a development version (`0.99.3`). The API may change before 1.0.0.

Local simulations found poor peak-test calibration, conservative transfer prediction intervals, and under-coverage of drift-fraction intervals in the tested scenarios. These methods remain exploratory. Results and reproducible scripts are documented in [the validation report](docs/VALIDATION_REPORT.md).

## Installation

From [R-universe](https://rotsl.r-universe.dev/crossStatR):

```r
install.packages(
  "crossStatR",
  repos = c(
    rotsl = "https://rotsl.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
```

Alternatively, from GitHub using `remotes`:

```r
install.packages("remotes")
remotes::install_github("rotsl/crossStatR")
```

Then verify it loads:

```r
library(crossStatR)
packageVersion("crossStatR")
#> [1] '0.99.3'
```

## Function catalogue

The package exports 32 user-facing functions, grouped into three categories.

| Category | Function | What it does |
|---|---|---|
| **A: Core statistics** | `sci_describe` | Descriptive statistics, with explicit biological vs. technical replicate handling |
| | `sci_compare` | Two-group comparisons (Welch, Student, paired, permutation, Wilcoxon) with effect sizes |
| | `sci_anova` | ANOVA, including one-way, two-way, and mixed-effects models |
| | `sci_regression` | Linear, weighted, and robust regression with diagnostic output |
| | `sci_correlation` | Correlation analysis with CIs (Pearson, Spearman, Kendall) |
| | `sci_normality` | Normality assessment with multiple tests and graphical diagnostics |
| | `sci_outliers` | Outlier diagnostics (IQR, MAD, Grubbs, Dixon) |
| | `sci_power` | Power and sample-size calculation for common designs |
| | `sci_effect_size` | Standardized effect sizes (Cohen's d, Hedges' g) with CIs |
| | `sci_bootstrap` | Bootstrap resampling with percentile, BCa, and normal CIs |
| **B: Domain methods** | `adaptive_baseline` | Asymmetric least squares or lowess baseline estimation for spectra |
| | `adaptive_smoothing` | Edge-preserving smoothing for noisy signals |
| | `batch_effect_quantify` | Quantify batch-to-batch variance, with mixed-effects decomposition |
| | `calibration_compare` | Compare two calibration curves via ANCOVA |
| | `curve_similarity` | Similarity between two curves sampled on possibly different grids |
| | `drift_correct_stat` | Drift correction with uncertainty quantification |
| | `kinetic_model_compare` | Compare candidate kinetic models via likelihood-ratio tests |
| | `noise_characterize` | Summarise residual variation and autocorrelation |
| | `parameter_uncertainty` | Propagate parameter uncertainty to derived quantities |
| | `replicate_decompose` | Variance-component decomposition for nested and crossed replicates |
| | `robust_lod` | Limits of detection and quantification (Currie-style) |
| | `spectral_change_test` | Test for differences across complete spectra |
| | `timepoint_optimize` | Optimise sampling timepoints for a planned experiment |
| | `weighted_peak_fit` | Weighted multi-peak fitting (Gaussian/Lorentzian) |
| **C: Advanced diagnostics** | `signal_difference_regions` | Identify contiguous regions where two groups of curves differ |
| | `peak_evidence` | Frequentist evidence (LR + permutation FPP) that a candidate feature is a true peak |
| | `calibration_transfer_risk` | Risk of transferring a calibration between instruments without recalibration |
| | `kinetic_discrimination` | Probability that a planned experiment can discriminate between kinetic models |
| | `experiment_fragility` | Sensitivity of a result to defensible analysis decisions |
| | `drift_uncertainty` | Decompose measurement variance into drift and noise components |
| | `detection_reliability` | Concentration-dependent PD/FPR/FNR curves with bootstrap bands |
| | `uncertainty_origin` | Stage-wise attribution of measurement uncertainty (Shapley) |

Every function returns an S3 object of class `crossstatr_<name>` with `print`, `summary`, `tidy`, and `as.data.frame` methods.

## Quick start

The two common starting functions are `sci_describe` (for summarizing a sample) and `sci_compare` (for comparing two groups).

```r
library(crossStatR)
set.seed(1)

# Biological assay: control vs. treatment, 8 biological replicates per group
ctrl <- rnorm(8, mean = 10, sd = 1.2)
trt  <- rnorm(8, mean = 11.5, sd = 1.2)

sci_describe(ctrl)
#> <crossStatR::describe> Descriptive statistics
#> Total observations: 8  Biological replicates: n/a
#>  group n     mean       sd      var        se   median      mad        cv
#>      1 8 10.15775 1.023833 1.048233 0.3619795 10.30789 1.214052 0.1007933
#>       iqr mean_ci_low mean_ci_high median_ci_low median_ci_high n_missing
#>  1.470133      9.3018     11.01369      8.997246       11.91434         0

sci_compare(ctrl, trt, test = "welch")
#> <crossStatR::compare> Comparison: welch
#> statistic = -2.3008,  p = 0.0387912
#> effect (mean(x) - mean(y)) = -1.4047
#> Effect sizes:
#>      name  estimate        se    ci_low    ci_high
#>   Cohen_d -1.150414 0.5601518 -2.351820 0.05099230
#>  Hedges_g -1.087664 0.5295980 -2.223539 0.04821091
```

## Examples by domain

### Biological: replicate-aware comparison

```r
set.seed(1)
df <- data.frame(
  y     = c(rnorm(12, 1.0), rnorm(12, 1.4)),
  grp   = factor(rep(c("ctrl", "trt"), each = 12)),
  bio   = factor(rep(1:12, times = 2)),
  tech  = factor(rep(1:3, times = 8))
)
# Summarise at the biological-replicate level first; do not treat
# technical replicates as independent observations.
d <- sci_describe(df$y, group = df$grp, bio_rep = df$bio, tech_rep = df$tech)
print(d)
c <- sci_compare(y ~ grp, data = df, test = "welch")
print(c)
```

### Electrochemical: drift-attributed variance

```r
set.seed(2)
t  <- seq(0, 10, length.out = 60)
y  <- 2 + 1.5 * t + 0.3 * (t - mean(t)) + rnorm(length(t), sd = 0.2)
dr <- drift_uncertainty(y, time = t, drift_model = "linear")
print(dr)
```

### Optical/spectroscopic: region-level differences

```r
set.seed(3)
wl <- seq(400, 800, length.out = 80)
g1 <- matrix(replicate(15, sin((wl - 400) * pi / 200) * 0.2 +
                          rnorm(80, sd = 0.15)), nrow = 15, byrow = TRUE)
g2 <- matrix(replicate(15, {
  base <- sin((wl - 400) * pi / 200) * 0.2 + dnorm(wl, 530, 18) * 0.8
  base + rnorm(80, sd = 0.15)
}), nrow = 15, byrow = TRUE)
sdr <- signal_difference_regions(wl, g1, g2, n_perm = 200)
print(sdr)
```

### Physics: uncertainty origin attribution

```r
# Four uncertainty sources, two correlated. Stage variances have
# already been estimated upstream (e.g. via replicate_decompose()
# and drift_uncertainty()).
stages <- c(sampling = 0.05, instrument = 0.10,
            calibration = 0.02, drift = 0.03)
uo <- uncertainty_origin(
  stages,
  correlations = list(sampling__instrument = 0.3,
                      calibration__drift = 0.2)
)
print(uo)
```

## Reproducibility

Randomness enters the package in three places: permutation tests, the bootstrap, and Monte Carlo simulation. Set the seed explicitly whenever the analyst wants a reproducible number.

```r
set.seed(42)
res <- sci_compare(x, y, test = "perm", seed = 42)
```

Seed handling varies by function. Calls that use a function-level seed may change the global random-number state; reproducible analyses should set and record seeds explicitly.

Beyond the seed, recording the analysis decisions that produced the reported result is at least as important. The documentation recommends running `experiment_fragility()` on the final result to see how sensitive the conclusion is to defensible choices: which replicates were included, how outliers were handled, which model was fit, and how much measurement noise would have to be added before the decision flips. If a single decision flips the sign of the effect, the claim is fragile and should be reported as such.

```r
af <- function(opts) {
  d <- opts$data
  if (!is.null(opts$rep_lo)) d <- d[d$rep != opts$rep_lo, ]
  diff(tapply(d$y, d$grp, mean))
}
fr <- experiment_fragility(
  analysis_fn = af, data = df, replicates = "bio",
  decision_threshold = 0, n_boot = 200, seed = 1
)
print(fr)
```

## Citation

The development version can be cited as follows:

> Rohan R (2026). crossStatR: Scientific Statistics for Biological, Electrochemical, Optical and Physics Experiments. R package version 0.99.3. https://github.com/rotsl/crossStatR

In R:

```r
citation("crossStatR")
```

## License

MIT. See [LICENSE.md](LICENSE.md) and [LICENSE](LICENSE). Copyright (c) 2026 Rohan R.

## Author

Rohan R.

Issues, bug reports, and pull requests are welcome at <https://github.com/rotsl/crossStatR>.


The test data is taken from Soedarmadji and Gregoire (2021), CaltechDATA:
<https://doi.org/10.25989/198p-1x22>.
`scripts/validate_caltech.R` reproduces the local checks and summaries.

Method descriptions and source references are available in [the methods guide](docs/METHODS.md).
