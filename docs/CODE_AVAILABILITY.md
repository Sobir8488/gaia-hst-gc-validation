# Code availability

The MATLAB workflow is provided in `src/` and `scripts/`. The main entry point is:

```matlab
run('scripts/run_gaia_hst_validation.m')
```

The archived CSV outputs can be checked independently with:

```bash
python3 scripts/check_csv_outputs.py
```

The source code release should be archived through Zenodo after a GitHub release is created. Record the resulting DOI in the manuscript data/code availability statement and in the GitHub repository description.
