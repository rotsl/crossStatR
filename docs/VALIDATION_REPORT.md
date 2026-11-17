# Validation report

Local run: 2026-11-17. Package version: 0.99.3 (development).
Environment: R version 4.6.2 (2026-10-29) on aarch64-apple-darwin23 .

The script `scripts/validate.R` runs 100 independent synthetic datasets per listed metric. Dataset seeds are 20261117 through 20261216; resampling seeds are the dataset indices. Intervals use 99 bootstrap draws and permutation tests use 99 permutations. MCSE is the sample standard deviation of the per-dataset metric divided by the square root of 100. Coverage averaged across concentrations is summarised at the dataset level. Zero observed events do not establish a zero underlying event probability (a rough 95% upper bound is 3/100).

| Method | Metric | Estimate | MCSE | Reference | Datasets |
|---|---|---:|---:|---|---:|
| signal_difference_regions | null rejection rate | 0.06000 | 0.02387 | 0.05 or lower | 100 |
| signal_difference_regions | power for uniform shift 0.4 | 1.00000 | 0.00000 | high | 100 |
| peak_evidence | null rejection at prespecified centre | 0.23000 | 0.04230 | 0.05 or lower | 100 |
| peak_evidence | detection at amplitude 1, noise SD 0.1 | 0.10000 | 0.03015 | high | 100 |
| calibration_transfer_risk | mean concentration bias | 0.00358 | 0.00401 | 0 | 100 |
| calibration_transfer_risk | mean pointwise 95% interval coverage | 0.98767 | 0.00235 | 0.95 | 100 |
| drift_uncertainty | drift fraction bias | 0.00014 | 0.00110 | 0 | 100 |
| drift_uncertainty | 95% interval coverage | 0.33000 | 0.04726 | 0.95 | 100 |
| detection_reliability | LOD bias | -0.00317 | 0.00283 | 0 | 100 |
| uncertainty_origin | MC total variance minus analytic variance | 0.01308 | 0.01883 | 0 | 100 |
| experiment_fragility | strong effect overall fragility | 0.00000 | 0.00000 | low | 100 |

## Kinetic discrimination

A separate 100-simulation experiment compares exponential and hyperbolic decay with 20 observations, amplitude 1, rate 0.2, and noise SD 0.05. Seed: 20261116.
```text
<crossStatR::kinetic_discrimination> Kinetic discrimination (true model: exponential)
Discrimination probabilities (criterion = aic, B = 100):
  exponential: 1.000 (MCSE 0.0000)
  hyperbolic: 0.000 (MCSE 0.0000)
All-candidate failures: 0 / 100

Condition numbers:
  exponential: kappa = 4.3
  hyperbolic: kappa = 4.2
```

## Interpretation

These are limited simulation diagnostics, not evidence of general calibration. Transfer prediction intervals now include independent source and transfer response-noise draws; those terms can overlap for paired noisy responses, and coverage remains scenario dependent. Drift coverage depends on the linear drift model. The peak procedure retains a smooth signal estimate in its permutation construction, so detection and null behaviour require independent methodological review. The uncertainty check compares Monte Carlo propagation with the supplied analytic covariance; it does not assess confidence intervals for estimated stage variances. Fragility measures sensitivity to the implemented perturbations, rather than a probability that a conclusion is correct. LOD bias is compared with the known positive-slope Gaussian calibration model. Larger simulations and independent comparisons remain necessary.

The [methods guide](METHODS.md) distinguishes implemented variants from the published background methods. The [reference list](REFERENCES.md) supplies source identifiers. Bibliographic verification does not resolve the calibration limitations reported here.

## Caltech measured-data checks

Local measured-data run: 2026-11-17. Package version: 0.99.3 (development).

Source: Soedarmadji, E., and Gregoire, J. (2021). Materials Experiment and Analysis Database: Electrochemistry analysis, Electrochemistry measurements, Material library with V,Ca,Ar,O on plate(s) 3540. CaltechDATA, v1. <https://doi.org/10.25989/198p-1x22>.

The supplied archive contains processed measurements. The run used all 29 CA2 and 29 CV3 light/dark current-difference tables. The test data is taken from <https://doi.org/10.25989/198p-1x22>. Basenames and input checksums are recorded in `data-raw/caltech-input-manifest.csv`; `scripts/validate_caltech.R` reproduces the checks and derived summaries.

| Analysis | Input rows | Excluded rows | Retained rows |
|---|---:|---:|---:|
| CA2 | 870 | 29 | 841 |
| CV3 | 2163 | 29 | 2134 |

All 348 checks passed across 58 curves and 2975 retained observations. Means, sample SDs, counts, and OLS coefficients matched base R; the reported difference SD matched `sd(diff(current))`. Smoothing returned finite values of the expected length. Current was converted from amperes to microamperes before comparisons. The comparison tolerance was 1e-10.

Rows were excluded only when the current difference or illuminated time was non-finite. No interpolation was performed. Several near-zero-mean current series produced the documented unstable-CV warning; the checks compare means and SDs rather than interpreting those CVs. Distinct material identifiers are not treatment replicates. Temporal dependence, processing by the source analysis, and real signal variation prevent these numerical agreement checks from establishing Type I error, interval coverage, physical model validity, or instrumental noise. The synthetic validation limitations above remain unresolved.

## Inference checks

Local run: 2026-11-17. Source: `scripts/validate_inference.R`.

The median, kurtosis, BCa and drift studies use initial seeds 170101, 170102, 170103 and 170104. BCa resampling uses seed i for dataset i, with 999 resamples. All coverage experiments have nominal level 95%. MCSE is computed from independent dataset-level metrics, including pointwise coverage averaged within each drift dataset.

| Method | Metric | Estimate | MCSE | Reference | Datasets |
|---|---|---:|---:|---:|---:|
| sci_describe | median 95% interval coverage, normal mean=10 n=30 | 0.95940 | 0.00279 | 0.95 | 5000 |
| sci_normality | mean excess kurtosis, normal n=100 | 0.00147 | 0.01041 | 0.00 | 2000 |
| sci_bootstrap | BCa 95% mean interval coverage, lognormal n=50 | 0.91600 | 0.00878 | 0.95 | 1000 |
| sci_bootstrap | percentile 95% mean interval coverage, lognormal n=50 | 0.91400 | 0.00887 | 0.95 | 1000 |
| drift_correct_stat | mean pointwise 95% interval coverage with independent standards | 0.95033 | 0.00297 | 0.95 | 500 |

### PlantGrowth numerical example

The built-in `datasets::PlantGrowth` example has group means ctrl = 5.032, trt1 = 4.661, trt2 = 5.526. The classical and Welch ANOVA p-values are 0.015910 and 0.017393; both matched base R. This numerical comparison does not establish interval coverage or a causal treatment conclusion. Dataset provenance is recorded in `help("PlantGrowth", package = "datasets")`.

BCa remains below nominal coverage in this skewed small-sample setting. Correcting the acceleration sign does not guarantee adequate coverage. Drift intervals combine reference-fit uncertainty and a shared noise estimate, assuming independent standards and samples with equal noise SD. Their coverage does not apply automatically to the running-median fallback or unequal reference/sample noise. Shapley precision intervals condition on supplied covariance and do not estimate uncertainty in experimental variance components.
