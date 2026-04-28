# XORI Data Analysis

Analysis of depth-dependent cross-orientation interactions in macaque V1 layer 2/3, measured with two-photon calcium imaging (PHP.eB-CAG-GCaMP6s). Orthogonal plaids are compared against a linear prediction (sum of component grating responses with baseline correction) to quantify normalization across approximately 400 um of cortical depth using nearly 5,000 ROIs detected with Suite2p across 28 fields of view.

**Bair Lab** | Department of Neurobiology and Biophysics, University of Washington, Seattle, WA

In collaboration with the Allen Institute for Brain Science and Washington National Primate Research Center.

---

## Background

Cross-orientation suppression (COS) is a hallmark of normalization in primary visual cortex: when two gratings at orthogonal orientations are superimposed into a plaid, the response of an orientation-selective neuron typically falls below the linear sum of its responses to each grating presented alone (Morrone et al., 1982; Bonds, 1989; Carandini et al., 1997; Freeman et al., 2002). This suppression is well described by divisive normalization models (Heeger, 1992; Carandini and Heeger, 2012) and has been studied extensively in cat and macaque V1 with single-unit electrophysiology (DeAngelis et al., 1992; Carandini et al., 1997; Busse et al., 2009; Smith et al., 2006).

However, nearly all prior work measured COS one neuron at a time, leaving open how suppression and facilitation are distributed across neural populations and whether they vary systematically with cortical depth. Layer 2/3 of macaque V1 spans several hundred micrometers and contains multiple sublayers with distinct connectivity patterns, cell-type distributions, and functional properties (Callaway, 1998; Briggs, 2010). Whether cross-orientation interactions change across this depth range has not been directly tested.

The long-term goal of this project is to combine large-scale functional imaging with dense connectomics to relate circuit wiring to population computation in primate V1. Using two-photon calcium imaging in anesthetized macaque, we target wide cortical blocks so that the physiology of as many neurons as possible can be compared post-mortem to serial EM reconstructions from the same tissue.

### Approach

Neurons are presented with drifting sine gratings and orthogonal plaids (each 50% contrast; 4 Hz; 4 cyc/deg; 2 deg patch), alongside RF maps from flashed light/dark spots and measurements of direction and spatial-frequency preferences. For each ROI, a linear prediction is built by shifting the single-grating tuning curve by -90 deg and summing with baseline correction. Deviation from this linear prediction is quantified with two complementary metrics:

- **S** (`M_S_ratio` in code): Ratio of the observed plaid response to the linear prediction. Values greater than 1 indicate facilitation (observed exceeds prediction); values less than 1 indicate suppression (observed falls below prediction). This is the primary metric used in current analyses. A legacy signed-difference version (`M_S`) is also computed but is no longer the focus.
- **R** (`M_C` in code): Pearson correlation between the predicted and observed plaid tuning curves, indexing shape similarity. Higher values indicate that plaid tuning is more linearly predictable from component grating responses.

---

## Key Findings

| Metric | Correlation with Depth | 95% CI | p-value |
|--------|----------------------|--------|---------|
| S (suppression/facilitation) | r = -0.768 | [-0.908, -0.617] | 1.83 x 10^-6 |
| R (shape similarity) | r = +0.798 | [+0.675, +0.891] | 3.67 x 10^-7 |

- **S decreases with depth**: positive (facilitation) in superficial layer 2/3, crossing zero, then negative (suppression) in deeper layer 2/3, consistent with stronger cross-orientation normalization at depth.
- **R increases with depth**: plaid tuning shapes become more linearly predictable deeper in layer 2/3, even as mean responses fall below the linear sum.
- Effects are robust to SNR thresholds and hold in both dF/F and raw fluorescence units.
- Mixed-effects models confirm the S-depth relationship at the single-ROI level (p = 7.9 x 10^-10).
- Spatial frequency preference mediates approximately 40% of the S-depth effect; tuning bandwidth mediates approximately 21%.

### Depth Profile Summary

Eight metrics show coordinated depth-dependent changes across layer 2/3:

| Deeper in Layer 2/3 | Direction |
|---------------------|-----------|
| ROI size (radius, npix) | Increases |
| Baseline fluorescence | Increases |
| SNR | Increases |
| Orientation tuning bandwidth | Broader |
| Preferred spatial frequency | Lower |
| Orientation selectivity (OSI) | Lower |
| Local homogeneity (LHI) | Lower |
| Cross-orientation suppression (S) | Stronger (more suppressive) |

Superficial layer 2/3 is characterized by high SF preference, narrow tuning, high OSI, high local homogeneity, and facilitation. Deeper layer 2/3 is characterized by low SF preference, broad tuning, low OSI, low local homogeneity, and suppression.

---

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
├── paper/
│   ├── manuscript.tex         # full draft
│   ├── manuscript_stats.json  # canonical stat lock
│   ├── make_figures.py        # generates fig1, fig3..fig9
│   ├── compute_stats.py       # prints LaTeX-formatted stats
│   ├── figures/               # PDF + PNG outputs
│   └── feedback/              # PI feedback artifacts (separate from manuscript)
├── docs/PIPELINE.md           # detailed operational pipeline doc
├── Makefile                   # one-command pipeline runner
├── CLAUDE.md                  # agent-facing project context
└── README.md                  # this file
```

Detailed pipeline reference: [`docs/PIPELINE.md`](docs/PIPELINE.md).
Project context for agent-driven work: [`CLAUDE.md`](CLAUDE.md).

---

## Metrics

All metrics are computed per ROI from grating and plaid tuning curve responses.

| Code name | Paper name | Description | Status |
|-----------|------------|-------------|--------|
| `M_S_ratio` | **S** | Ratio of observed plaid response to linear prediction. >1 = facilitation, <1 = suppression. | **Primary** |
| `M_C` | **R** | Pearson correlation between predicted and observed plaid tuning curves (shape similarity). | **Primary** |
| `M_S` | | Mean signed difference (observed minus predicted). Legacy version of S, retained for comparison. | Legacy |
| `M_S_norm` | | Signed difference normalized as percentage of baseline fluorescence. | Legacy |
| `M_X` | | Additional cross-orientation suppression metric. | Secondary |
| `SNR_g` | | Signal-to-noise ratio for grating responses. | Quality filter |
| `SNR_p` | | Signal-to-noise ratio for plaid responses. | Quality filter |

---

## Setup

```bash
# preferred (uv-managed Python 3.12 venv)
make env

# or manual
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running Analyses

The Makefile drives the whole pipeline. Each target also lists its underlying
script for direct invocation when needed.

```bash
make metrics       # → scripts/m_calc/all_metric.py
make culls         # → scripts/m_calc/{cull_metric,r_cull_metric}.py
make depth         # → scripts/d_calc/depth.py
make covariates    # → scripts/bm_calc/{baseline,halfwidth,lhi,osi,size,spatial}.py
make supp          # → supplementary_analysis/scripts/*.py
make stats         # → paper/compute_stats.py (prints to stdout)
make figures       # → paper/make_figures.py (fig1, fig3..fig9)
make manuscript    # → latexmk paper/manuscript.tex (requires basictex/mactex)
make all           # full pipeline end-to-end
```

For direct script-level operation, see `docs/PIPELINE.md`.

### Supplementary statistical analyses
```bash
# Subpopulation analysis, partial correlations, bootstrap CIs
python supplementary_analysis/scripts/run_all_analyses.py

# Mixed-effects models, ROI size confound controls, publication figure
python supplementary_analysis/scripts/additional_analyses.py

# Multi-metric depth profiles, mediation analysis
python supplementary_analysis/scripts/extended_analyses.py
```

---

## Methods

See [supplementary_analysis/METHODS.md](supplementary_analysis/METHODS.md) for a draft methods section with statistical details, key statistics tables, and suggested figure legends.

---

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

---

## References

Bonds, A. B. (1989). Role of inhibition in the specification of orientation selectivity of cells in the cat striate cortex. *Visual Neuroscience*, 2(1), 41-55.

Briggs, F. (2010). Organizing principles of cortical layer 6. *Frontiers in Neural Circuits*, 4, 3.

Busse, L., Wade, A. R., and Carandini, M. (2009). Representation of concurrent stimuli by population activity in visual cortex. *Neuron*, 64(6), 931-942.

Callaway, E. M. (1998). Local circuits in primary visual cortex of the macaque monkey. *Annual Review of Neuroscience*, 21, 47-74.

Carandini, M., and Heeger, D. J. (2012). Normalization as a canonical neural computation. *Nature Reviews Neuroscience*, 13(1), 51-62.

Carandini, M., Heeger, D. J., and Movshon, J. A. (1997). Linearity and normalization in simple cells of the macaque primary visual cortex. *Journal of Neuroscience*, 17(21), 8621-8644.

DeAngelis, G. C., Robson, J. G., Ohzawa, I., and Freeman, R. D. (1992). Organization of suppression in receptive fields of neurons in cat visual cortex. *Journal of Neurophysiology*, 68(1), 144-163.

Freeman, T. C. B., Durand, S., Kiper, D. C., and Carandini, M. (2002). Suppression without inhibition in visual cortex. *Neuron*, 35(4), 759-771.

Heeger, D. J. (1992). Normalization of cell responses in cat striate cortex. *Visual Neuroscience*, 9(2), 181-197.

Morrone, M. C., Burr, D. C., and Maffei, L. (1982). Functional implications of cross-orientation inhibition of cortical visual cells. I. Neurophysiological evidence. *Proceedings of the Royal Society of London B*, 216(1204), 335-354.

Pachitariu, M., Stringer, C., Dipoppa, M., Schroeder, S., Rossi, L. F., Dalgleish, H., Carandini, M., and Harris, K. D. (2017). Suite2p: beyond 10,000 neurons with standard two-photon microscopy. *bioRxiv*, 061507.

Ringach, D. L., Shapley, R. M., and Hawken, M. J. (2002). Orientation selectivity in macaque V1: diversity and laminar dependence. *Journal of Neuroscience*, 22(13), 5639-5651.

Smith, M. A., Bair, W., and Movshon, J. A. (2006). Dynamics of suppression in macaque primary visual cortex. *Journal of Neuroscience*, 26(18), 4826-4834.

---

## Dependencies

- Python 3.10+
- NumPy, Pandas, Matplotlib, SciPy, statsmodels
- Suite2p (for ROI detection, run separately)

## License

This repository contains research data and analysis code. Contact the lab before reuse.
