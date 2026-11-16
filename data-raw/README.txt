The script scripts/validate_caltech.R reads an extracted CaltechDATA archive
from a directory passed on the command line. It writes derived diagnostics,
a checksum manifest, numerical check totals, and runtime measurements.

Rscript scripts/validate_caltech.R <extracted-data-directory>

Citation: Soedarmadji, E., and Gregoire, J. (2021). Materials Experiment and
Analysis Database: Electrochemistry analysis, Electrochemistry measurements,
Material library with V,Ca,Ar,O on plate(s) 3540. CaltechDATA, v1.
https://doi.org/10.25989/198p-1x22

caltech-input-manifest.csv lists basenames and MD5 checksums for the 58
upstream analysis tables used in the local run. The test data is taken from https://doi.org/10.25989/198p-1x22.
