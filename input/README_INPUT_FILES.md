# Input files

This directory contains the input CSV tables used by the MATLAB workflow. The principal input layers are Gaia EDR3 profile-level products, HST/HACKS proper-motion dispersion profiles, APOGEE DR17 globular-cluster support metadata, a cluster crosswalk, and analysis configuration.

The raw HACKS text file is stored under `input/raw/`. The conversion script is `src/convert_hacks_raw_to_csv.m`.

See `metadata/source_provenance.csv` and `metadata/column_dictionary.csv` for data-layer provenance, column definitions, units, and missing-value handling.
