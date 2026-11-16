
test_that("Package loads", {
  expect_true(requireNamespace("crossStatR", quietly = TRUE))
})

test_that("All 32 public functions are exported", {
  expected_exports <- c(
    "sci_describe", "sci_normality", "sci_outliers", "sci_compare",
    "sci_anova", "sci_regression", "sci_correlation", "sci_bootstrap",
    "sci_power", "sci_effect_size",
    "adaptive_baseline", "weighted_peak_fit", "adaptive_smoothing",
    "replicate_decompose", "drift_correct_stat", "robust_lod",
    "calibration_compare", "curve_similarity", "kinetic_model_compare",
    "spectral_change_test", "noise_characterize", "parameter_uncertainty",
    "timepoint_optimize", "batch_effect_quantify",
    "signal_difference_regions", "peak_evidence", "calibration_transfer_risk",
    "kinetic_discrimination", "experiment_fragility", "drift_uncertainty",
    "detection_reliability", "uncertainty_origin"
  )
  ns_exports <- getNamespaceExports("crossStatR")
  for (fn in expected_exports) {
    expect_true(fn %in% ns_exports, info = paste("missing export:", fn))
  }
})

test_that("S3 generics are exported or registered", {
  # S3 methods registered via S3method() are loaded into the methods table
  # when the package is loaded, not when the namespace is queried.
  # The cleanest test is to actually dispatch.
  obj <- structure(list(method = "test", summary = data.frame(a = 1)),
                  class = c("crossstatr_dummy", "crossstatr_result"))
  # tidy() should work because tidy.crossstatr_result is exported
  expect_true(exists("tidy.crossstatr_result", envir = asNamespace("crossStatR")))
  # Print / summary / as.data.frame methods exist
  expect_true(exists("print.crossstatr_result", envir = asNamespace("crossStatR")))
  expect_true(exists("summary.crossstatr_result", envir = asNamespace("crossStatR")))
  expect_true(exists("as.data.frame.crossstatr_result", envir = asNamespace("crossStatR")))
})
