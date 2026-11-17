Caltech electrochemistry diagnostics

Citation: Soedarmadji, Edwin, and Gregoire, John (2021).
Materials Experiment and Analysis Database: Electrochemistry analysis,
Electrochemistry measurements, Material library with V,Ca,Ar,O on plate(s) 3540.
CaltechDATA, version v1. https://doi.org/10.25989/198p-1x22
Record: https://data.caltech.edu/records/bn6dh-p0s22

caltech-current-summary.csv contains derived numerical diagnostics for all
29 CA2 and 29 CV3 interlen text tables in the supplied analysis archive.
Test data is taken from https://doi.org/10.25989/198p-1x22.
The archive is identified by
20161101.114257.copied-20161107195737833PST.zip on the cited record.

analysis: upstream stage (1 = CA2, 2 = CV3).
sample: upstream material sample identifier; not a replicate identifier.
input_rows, excluded_rows, observations: original, excluded, retained counts.
Rows are excluded only when I(A)_illdiff or t(s)_ill is non-finite.
mean_photocurrent_uA, sd_photocurrent_uA: mean and sample SD of I(A)_illdiff
multiplied by 1e6. These are upstream processed light/dark current differences,
not reconstructed differences between the other columns.
slope_uA_per_s: OLS slope against t(s)_ill.
difference_noise_uA: sample SD of successive current differences, without
rescaling by sqrt(2); genuine signal variation can contribute to this metric.
lag1_acf: lag-one autocorrelation of those successive differences.

These diagnostics test numerical behaviour on measured inputs. They do not
establish independent replication, causal effects, or inferential calibration.
The package MIT licence does not assert ownership of the upstream dataset.
