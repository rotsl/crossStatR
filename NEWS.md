# crossStatR 0.99.3

- Refreshed local validation results, benchmarks, and documentation.

# crossStatR 0.99.2

- Updated R-universe installation documentation.

# crossStatR 0.99.1

Initial public development version. This is a pre-release: the API may change before 1.0.0. The largest changes are expected in the Category C methods (signal-difference regions, peak evidence, calibration transfer risk, kinetic discrimination, experiment fragility, drift uncertainty, detection reliability, and uncertainty-origin attribution), where the interface is still settling. Version pinning supports reproducible analyses.

## Added

32 user-facing functions, grouped into three categories.

### Category A: core statistics

- `sci_describe` - descriptive statistics with explicit biological vs. technical replicate handling
- `sci_compare` - two-group comparisons (Welch, Student, paired, permutation, Wilcoxon) with effect sizes
- `sci_anova` - ANOVA and Welch alternatives with effect sizes
- `sci_regression` - linear, weighted, and robust regression with diagnostics
- `sci_correlation` - correlation analysis with confidence intervals (Pearson, Spearman, Kendall)
- `sci_normality` - normality assessment with multiple tests
- `sci_outliers` - outlier diagnostics (IQR, MAD, z-score, robust z-score)
- `sci_power` - power and sample-size calculation for common designs
- `sci_effect_size` - standardized effect sizes (Cohen's d, Hedges' g) with CIs
- `sci_bootstrap` - bootstrap resampling with percentile, BCa, and normal CIs

### Category B: domain methods

- `adaptive_baseline` - asymmetric least squares / lowess baseline estimation for spectra
- `adaptive_smoothing` - edge-preserving smoothing for noisy signals
- `batch_effect_quantify` - quantify batch-to-batch variance via mixed-effects decomposition
- `calibration_compare` - compare two calibration curves via ANCOVA
- `curve_similarity` - similarity between two curves sampled on possibly different grids
- `drift_correct_stat` - drift correction with uncertainty quantification
- `kinetic_model_compare` - compare candidate kinetic models via likelihood-ratio tests
- `noise_characterize` - summarise residual variation and autocorrelation
- `parameter_uncertainty` - propagate parameter uncertainty to derived quantities
- `replicate_decompose` - variance-component decomposition for nested and crossed replicates
- `robust_lod` - limits of detection and quantification (Currie-style)
- `spectral_change_test` - test for differences across complete spectra
- `timepoint_optimize` - optimise sampling timepoints for a planned experiment
- `weighted_peak_fit` - weighted multi-peak fitting (Gaussian/Lorentzian)

### Category C: advanced diagnostics

- `signal_difference_regions` - statistical identification of regions where two groups of curves differ (simultaneous max-t band with Westfall-Young permutation)
- `peak_evidence` - frequentist evidence assessment for signal peaks (likelihood ratio + permutation FPP)
- `calibration_transfer_risk` - statistical framework for calibration transfer risk (TPE, OOD flag, parametric bootstrap)
- `kinetic_discrimination` - discriminability of competing kinetic models via Monte Carlo simulation
- `experiment_fragility` - quantitative experimental fragility assessment (leave-one-replicate-out, outlier, baseline, model, noise)
- `drift_uncertainty` - statistical attribution of variance to instrument drift (linear/polynomial/exponential drift models)
- `detection_reliability` - reliability analysis for analytical detection (PD/FNR/FPR curves with bootstrap CIs)
- `uncertainty_origin` - stage-wise uncertainty attribution via Shapley value for correlated sources

### Infrastructure

- S3 class system for results: `crossstatr_result` parent class with `print`, `summary`, `tidy`, `as.data.frame` methods
- Internal validation helpers (`.assert_*`) and numeric helpers (`.safe_sd`, `.ci_t`, `.cohens_d`, `.mcse_p`, etc.)
- `testthat` (3rd edition) test suite covering input validation, numerical correctness, and edge cases
- MIT license, `inst/CITATION`, `inst/extdata/` for example data
