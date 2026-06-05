GaiaHST A&C supplementary CSV package
=====================================

This folder contains the CSV products that are most important to upload with the manuscript
"A Reproducible Gaia--HST Internal-Kinematics Validation Pipeline for Galactic Globular Clusters".

Primary files for the Appendix / data availability statement:
1. GaiaHST_MATCHED_RADIAL_COMPARISON.csv
2. GaiaHST_CLUSTER_QUALITY_CLASSES.csv
3. GaiaHST_GLOBAL_CALIBRATION_SUMMARY.csv
4. GaiaHST_DATA_LAYER_READINESS.csv
5. GaiaHST_ROBUSTNESS_SUMMARY.csv
6. GaiaHST_BOOTSTRAP_GLOBAL_MEAN.csv
7. GaiaHST_CLASS_BOUNDARY_SENSITIVITY.csv

Additional manuscript-support CSV files are included for input-layer summary, benchmark cases,
caution clusters, APOGEE overlap support and class-count summaries.

Expected manuscript checks:
- Matched radial bins: 680
- Validation-ready clusters: 56
- Global mean Delta_sigma: -0.002846
- Global RMS Delta_sigma: 0.132242
- Readiness classes: A=38, B=8, C=8, D=2, E=0

The table GaiaHST_CSV_MANIFEST.csv lists every CSV file, row count, column count, columns and role.