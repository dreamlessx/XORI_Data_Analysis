# XORI Data Analysis

Analysis of depth-dependent cross-orientation interactions in macaque V1 layer 2/3, measured with two-photon calcium imaging (PHP.eB-CAG-GCaMP6s). Orthogonal plaids are compared against a linear prediction (sum of component grating responses with baseline correction) to quantify normalization across ~400 μm of cortical depth using ~5,000 ROIs detected with Suite2p across 28 fields of view.

**Bair Lab** — Department of Neurobiology & Biophysics, University of Washington, Seattle, WA

## Background

The long-term goal is to combine large-scale functional imaging with dense connectomics to relate circuit wiring to population computation in primate V1. Neurons are presented with drifting sine gratings and orthogonal plaids (each 50% contrast; 4 Hz; 4 cyc/deg; 2° patch). For each ROI, a linear reference is built by shifting the single-grating tuning curve by -90° and summing with baseline correction. Deviation from this linear prediction is quantified with two complementary metrics:

- **S** (`M_S_ratio` in code) — Ratio of observed plaid response to linear prediction. Values >1 indicate facilitation (observed exceeds prediction); values <1 indicate suppression (observed falls below prediction). This is the primary metric used in current analyses. (A legacy signed-difference version, `M_S`, is also computed but no longer the focus.)
- **R** (`M_C` in code) — Pearson correlation between predicted and observed plaid tuning curves, indexing shape similarity. Higher values indicate that plaid tuning is more linearly predictable from component grating responses.

## Key Findings

| Metric | Correlation with Depth | 95% CI | p-value |
|--------|----------------------|--------|---------|
| S (suppression/facilitation) | r = -0.768 | [-0.908, -0.617] | 1.83 x 10^-6 |
| R (shape similarity) | r = +0.798 | [+0.675, +0.891] | 3.67 x 10^-7 |

- **S decreases with depth**: Positive (facilitation) in superficial layer 2/3, crossing zero, then negative (suppression) in deeper layer 2/3 -- consistent with stronger cross-orientation normalization at depth
- **R increases with depth**: Plaid tuning shapes become more linearly predictable deeper in layer 2/3, even as mean responses fall below the linear sum
- Effects are robust to SNR thresholds and hold in both dF/F and raw fluorescence units
- Mixed-effects models confirm S-depth relationship at the single-ROI level (p = 7.9 x 10^-10)
- Spatial frequency preference mediates ~40% of the S-depth effect; tuning bandwidth mediates ~21%

## Repository Structure

```
XORI/
├── raw_data/                  # Raw experimental data
│   ├── bm_data/               # ROI-level measurements
│   │   ├── site_depth.txt     # Field-of-view to depth mapping (28 sites)
│   │   ├── roi_osi.txt        # Orientation/direction selectivity (4,785 ROIs)
│   │   ├── roi_stat.txt       # ROI morphology from Suite2p (radius, aspect ratio, etc.)
│   │   ├── roi_hw_orth.txt    # Orientation tuning half-width
│   │   ├── roi_lhi.txt        # Local homogeneity index
│   │   └── roi_sf.txt         # Spatial frequency preferences
│   └── tc_data/               # Tuning curve response data per site
├── metric_data/               # Computed cross-orientation metrics (S, R) per site
│   ├── all_roi/               # All ROIs (primary analysis)
│   ├── cull_roi/              # High-SNR subsets (top 70/80/90%)
│   └── r_cull_roi/            # Low-SNR subsets (quality control)
├── depth_data/                # Depth analysis outputs (plots, ROI maps)
├── data_baseline/             # Baseline fluorescence vs depth
├── data_halfwidth/            # Tuning half-width vs depth
├── data_spatial/              # Spatial frequency preference vs depth
├── data_osi/                  # Orientation selectivity vs depth
├── data_size/                 # ROI size confound analysis
├── data_lhi/                  # Local homogeneity vs depth
├── scripts/                   # Core analysis scripts
│   ├── m_calc/                # Metric calculation from tuning curves
│   ├── d_calc/                # Depth correlation analysis
│   └── bm_calc/               # Covariate-metric relationship scripts
├── supplementary_analysis/    # Extended statistical analyses
│   ├── scripts/
│   │   ├── run_all_analyses.py        # Subpopulation splits, partial correlations, bootstrap
│   │   ├── additional_analyses.py     # Mixed-effects models, ROI size controls
│   │   └── extended_analyses.py       # Multi-metric depth profiles, mediation analysis
│   ├── outputs/               # Figures and results
│   └── METHODS.md             # Draft methods section for publication
├── stat/                      # Suite2p stat files (ROI spatial masks)
└── zz_Playground/             # Development workspace
```

## Metrics

All metrics are computed per ROI from grating and plaid tuning curve responses.

| Code name | Paper name | Description | Status |
|-----------|------------|-------------|--------|
| `M_S_ratio` | **S** | Ratio of observed plaid response to linear prediction. >1 = facilitation, <1 = suppression. | **Primary metric** |
| `M_C` | **R** | Pearson correlation between predicted and observed plaid tuning curves (shape similarity). | **Primary metric** |
| `M_S` | -- | Mean signed difference (observed - predicted). Legacy version of S, retained for comparison. | Legacy |
| `M_S_norm` | -- | Signed difference normalized as percentage of baseline fluorescence. | Legacy |
| `M_X` | -- | Additional cross-orientation metric. | Secondary |
| `SNR_g` | -- | Signal-to-noise ratio for grating responses. | Quality filter |
| `SNR_p` | -- | Signal-to-noise ratio for plaid responses. | Quality filter |

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running Analyses

### Metric calculation from tuning curves
```bash
python scripts/m_calc/all_metric.py
```

### Core depth analysis (ROI maps + depth vs metric plots)
```bash
python scripts/d_calc/depth.py
```

### Covariate analyses (baseline, halfwidth, SF, OSI, LHI, ROI size)
```bash
python scripts/bm_calc/baseline.py
python scripts/bm_calc/halfwidth.py
python scripts/bm_calc/spatial.py
python scripts/bm_calc/osi.py
python scripts/bm_calc/lhi.py
python scripts/bm_calc/size.py
```

### Supplementary statistical analyses
```bash
# Subpopulation analysis, partial correlations, bootstrap CIs
python supplementary_analysis/scripts/run_all_analyses.py

# Mixed-effects models, ROI size confound controls, publication figure
python supplementary_analysis/scripts/additional_analyses.py

# Multi-metric depth profiles, mediation analysis
python supplementary_analysis/scripts/extended_analyses.py
```

## Methods

See [supplementary_analysis/METHODS.md](supplementary_analysis/METHODS.md) for a draft methods section with statistical details, key statistics tables, and suggested figure legends.

## Data Format

Cross-orientation metric files (`metric_data/all_roi/metrics_siteXXX.txt`):
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
- NumPy, Pandas, Matplotlib, SciPy, statsmodels
- Suite2p (for ROI detection, run separately)

## License

This repository contains research data and analysis code. Contact the lab before reuse.
