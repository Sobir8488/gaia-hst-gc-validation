#!/usr/bin/env python3
"""Check row counts and principal summary statistics for the archived CSV outputs."""
from pathlib import Path
import csv, math, statistics, sys

root = Path(__file__).resolve().parents[1]
out = root / 'output_csv'
matched_file = out / 'GaiaHST_MATCHED_RADIAL_COMPARISON.csv'
classes_file = out / 'GaiaHST_CLUSTER_QUALITY_CLASSES.csv'

def read_csv(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

matched = read_csv(matched_file)
classes = read_csv(classes_file)
delta = [float(r['delta_frac']) for r in matched if r.get('delta_frac','') not in ('', 'NaN', 'nan')]
clusters = sorted({r['cluster_id'] for r in matched})
counts = {}
for r in classes:
    q = r['quality_class'].split('_', 1)[0]
    counts[q] = counts.get(q, 0) + 1
mean_delta = sum(delta)/len(delta)
median_delta = statistics.median(delta)
rms_delta = math.sqrt(sum(x*x for x in delta)/len(delta))
print('matched_bins', len(matched))
print('matched_clusters', len(clusters))
print('mean_delta_frac', f'{mean_delta:.6f}')
print('median_delta_frac', f'{median_delta:.6f}')
print('rms_delta_frac', f'{rms_delta:.6f}')
print('class_counts', counts)
expected = {
    'matched_bins': 680,
    'matched_clusters': 56,
    'class_counts': {'A':38, 'B':8, 'C':8, 'D':2},
}
if len(matched) != expected['matched_bins'] or len(clusters) != expected['matched_clusters'] or counts != expected['class_counts']:
    sys.exit('CSV output check failed.')
print('status OK')
