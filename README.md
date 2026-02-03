# XORI Data Analysis

Depth-dependent cross-orientation suppression in macaque V1 measured via two-photon calcium imaging. Analyzes how cross-orientation interaction metrics (M_S, M_C) vary across cortical layers using 4,785 ROIs from 28 recording sites spanning 140–518 μm depth.

## Key Findings

| Metric | Correlation with Depth | 95% CI | p-value |
|--------|----------------------|--------|---------|
| M_S (additive component) | r = −0.768 | [−0.908, −0.617] | 1.83 × 10⁻⁶ |
| M_C (response correlation) | r = +0.798 | [+0.675, +0.891] | 3.67 × 10⁻⁷ |

- M_S (cross-orientation suppression) decreases with cortical depth — superficial layers show facilitation, deep layers show suppression
- M_C (grating–plaid response correlation) increases with depth
- Effects survive partial correlation controls for ROI morphology, SNR, OSI, and spatial frequency
- Mixed-effects models confirm M_S–depth relationship at the single-ROI level (p = 7.9 × 10⁻¹⁰)

## Repository Structure

```
XORI/
├── raw_data/                  # Raw experimental data
│   ├── bm_data/               # Behavioral/metric source files
│   │   ├── site_depth.txt     # Site-to-depth mapping (28 sites)
│   │   ├── roi_osi.txt        # Orientation selectivity data (4,785 ROIs)
│   │   ├── roi_stat.txt       # ROI morphology (radius, aspect ratio, etc.)
│   │   ├── roi_hw_orth.txt    # Half-width tuning data
│   │   ├── roi_lhi.txt        # Local homogeneity index
│   │   └── roi_sf.txt         # Spatial frequency preferences
│   └── tc_data/               # Tuning curve data per site
├── metric_data/               # Computed cross-orientation metrics per site
│   ├── all_roi/               # All ROIs (primary analysis)
│   ├── cull_roi/              # SNR-filtered subsets (top 70/80/90%)
│   └── r_cull_roi/            # Low-SNR subsets (quality control)
├── depth_data/                # Depth analysis outputs (plots, ROI maps)
├── data_baseline/             # Baseline fluorescence analysis
├── data_halfwidth/            # Tuning half-width analysis
├── data_spatial/              # Spatial frequency analysis
├── data_osi/                  # Orientation selectivity analysis
├── data_size/                 # ROI size control analysis
├── data_lhi/                  # Local homogeneity analysis
├── scripts/                   # Core analysis scripts
│   ├── m_calc/                # Metric calculation
│   ├── d_calc/                # Depth correlation analysis
│   └── bm_calc/               # Baseline/metric relationship scripts
├── supplementary_analysis/    # Extended statistical analyses
│   ├── scripts/               # Supplementary analysis scripts
│   │   ├── run_all_analyses.py        # Subpopulation, partial correlations, bootstrap
│   │   ├── additional_analyses.py     # Mixed-effects models, ROI size controls
│   │   └── extended_analyses.py       # Multi-metric profiles, mediation analysis
│   ├── outputs/               # Figures and results
│   └── METHODS.md             # Draft methods section for publication
├── stat/                      # Suite2p stat files (ROI masks)
└── zz_Playground/             # Development workspace
```

## Metrics

- **M_S** — Additive component of cross-orientation interaction (plaid − grating response). Positive = facilitation, negative = suppression.
- **M_C** — Correlation between grating and plaid responses across conditions.
- **M_S_norm** — M_S normalized as percentage of baseline response.
- **M_S_ratio** — Ratio of plaid to grating response.
- **M_X** — Cross-orientation metric (additional measure).
- **SNR_g / SNR_p** — Signal-to-noise ratios for grating and plaid responses.

## Setup

```bash
# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install numpy pandas matplotlib scipy statsmodels
```

## Running Analyses

### Core depth analysis
```bash
python scripts/d_calc/depth.py
```

### Baseline, halfwidth, spatial frequency, OSI analyses
```bash
python scripts/bm_calc/baseline.py
python scripts/bm_calc/halfwidth.py
python scripts/bm_calc/spatial.py
python scripts/bm_calc/osi.py
```

### Supplementary analyses
```bash
# Subpopulation analysis, partial correlations, bootstrap
python supplementary_analysis/scripts/run_all_analyses.py

# Mixed-effects models, ROI size confound controls
python supplementary_analysis/scripts/additional_analyses.py

# Multi-metric depth profiles, mediation analysis
python supplementary_analysis/scripts/extended_analyses.py
```

## Methods

See [supplementary_analysis/METHODS.md](supplementary_analysis/METHODS.md) for a draft methods section with full statistical details, key statistics tables, and suggested figure legends.

## Data Format

Input metric files (`metric_data/all_roi/metrics_siteXXX.txt`):
```
ROI    M_S         M_C        SNR_g      SNR_p      M_S_norm   M_S_ratio   M_X
0      -43.738     0.51768    3.4882     0.93075    -15.137    0.66355     ...
1      38.941      -0.087577  0.52251    0.11513    9.7861     2.4277      ...
```

Site depth mapping (`raw_data/bm_data/site_depth.txt`):
```
site        depth
site002     168
site003     160
...
```

## Dependencies

- Python 3.10+
- NumPy
- Pandas
- Matplotlib
- SciPy
- statsmodels

## Author

Mudit Agar (rajagar@uw.edu)

## License

This repository contains research data and analysis code. Please contact the author before reuse.
