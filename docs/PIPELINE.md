# XORI Data Pipeline

Operational reference for the XORI cross-orientation analysis pipeline. This is the
detailed pipeline-level doc. For project orientation and headline results, see the
top-level `README.md`.

## Project Overview

This project analyzes depth-dependent patterns in cross-orientation interactions
in macaque V1 layer 2/3 across 4,785 ROIs from 28 fields of view (140-518 um),
measured with two-photon calcium imaging.

The pipeline computes:
- Cross-orientation metrics per ROI: P (`M_S_ratio`), C (`M_C`), plus M_S, M_S_norm, M_X, SNR_g, SNR_p
- Per-site covariate-vs-metric analyses: baseline (log10), bandwidth (HW), SF (log10), OSI (Gaku and circular variance), LHI (2D and 3D), ROI size
- Depth-correlation analyses with multiple SNR-cull regimes
- Supplementary statistical analyses: mixed-effects models, partial correlations, mediation, bootstrap CIs

## Directory Structure

```
xori/
├── raw_data/
│   ├── bm_data/                # ROI-level inputs and site metadata
│   │   ├── roi_hw_orth.txt     # orientation half-width
│   │   ├── roi_lhi.txt         # local homogeneity index (2D and 3D)
│   │   ├── roi_osi.txt         # orientation selectivity (Gaku + circular variance)
│   │   ├── roi_sf.txt          # spatial frequency preference
│   │   ├── roi_stat.txt        # Suite2p morphology summary
│   │   └── site_depth.txt      # FOV depth mapping (28 sites)
│   └── tc_data/                # per-site directories of per-ROI tuning curves (xplot files)
├── metric_data/                # cross-orientation metrics computed from tc_data
│   ├── all_roi/                # 28 metrics_siteXXX.txt files (PRIMARY analysis)
│   ├── cull_roi/               # SNR-thresholded subsets
│   │   ├── per_cull/{top_70,top_80,top_90}/    # percentile splits, 28 sites each
│   │   └── thr_cull/{above_0_5,above_1_0,above_1_5}/  # 27 each (site038 legitimately empty at all thresholds)
│   └── r_cull_roi/             # reverse cull (low-SNR), 28 each
├── stat/                       # Suite2p stat_siteXXX.npy (28 files, ROI spatial masks)
├── depth_data/                 # depth correlation analysis outputs (ROI maps + scatter)
│   ├── all_roi/                # all ROIs, no filtering
│   ├── cull_roi/{per_cull,thr_cull}/
│   ├── null_roi/
│   └── r_cull_roi/{per_r_cull,thr_r_cull}/
├── data_baseline/log10/all_roi/      # baseline-vs-metrics per site + summary
├── data_halfwidth/raw/all_roi/       # bandwidth-vs-metrics per site + summary
├── data_lhi/{2d,3d}/                 # LHI-vs-metrics per site + summary
├── data_osi/                         # OSI-vs-metrics
│   ├── osi/                          # Gaku's depth-of-modulation OSI
│   └── variance/                     # circular-variance OSI (1 − ocv)
├── data_size/null_roi/               # ROI-size-vs-metrics
├── data_spatial/{all_roi,all_roi_log10}/    # SF-vs-metrics (raw and log10)
├── scripts/                          # analysis pipeline
│   ├── m_calc/
│   │   ├── all_metric.py             # primary metric calc from tc_data
│   │   ├── cull_metric.py            # SNR-cull subsets
│   │   └── r_cull_metric.py          # reverse cull
│   ├── d_calc/
│   │   └── depth.py                  # depth correlation analysis
│   └── bm_calc/
│       ├── baseline.py
│       ├── halfwidth.py
│       ├── lhi.py
│       ├── osi.py
│       ├── size.py
│       └── spatial.py
├── supplementary_analysis/
│   ├── METHODS.md                    # draft methods detail
│   ├── README.md
│   ├── outputs/                      # 17 result subdirs (mixed_effects, mediation, etc.)
│   └── scripts/
│       ├── run_all_analyses.py       # subpopulation, partial corr, bootstrap CIs
│       ├── additional_analyses.py    # mixed-effects, ROI size confound
│       └── extended_analyses.py      # multi-metric depth profile, mediation
├── paper/                            # manuscript + publication figures
│   ├── manuscript.tex                # 755 lines, full draft
│   ├── manuscript.pdf
│   ├── manuscript_stats.json         # CANONICAL numbers (do not regenerate casually)
│   ├── references.bib                # 32 cites, abbrvnat
│   ├── make_figures.py               # generates fig1, fig3..fig9
│   ├── compute_stats.py              # prints LaTeX-formatted stats to stdout
│   ├── figures/                      # PDF + PNG for fig1, fig3..fig9
│   └── feedback/                     # PI feedback artifacts (not part of manuscript)
│       ├── feedback.{md,pdf,tex}
│       ├── make_feedback_figures.py
│       └── figures/                  # feedback_*.{pdf,png}
├── docs/
│   └── PIPELINE.md                   # this file
├── Makefile                          # one-command pipeline runner
├── README.md                         # publication-facing project overview
├── requirements.txt
└── .venv/                            # uv-managed Python 3.12 (gitignored)
```

### Site038 caveat

Site038, the shallowest site at 140 um, contains 134 ROIs but **all** have SNR_g < 0.5.
The thresholded-cull script (`scripts/m_calc/cull_metric.py:240`) correctly logs
`"no ROIs meet threshold (skipped)"` and writes no output for this site under any
of `thr_cull/{above_0_5, above_1_0, above_1_5}/`. This is documented behavior, not
data exclusion. Site038 is fully retained in `metric_data/all_roi/`, the percentile
culls (`per_cull/top_70/top_80/top_90`), the reverse culls, and all depth and
covariate analyses.

## Data Files

### Input Data Format

#### `site_depth.txt`
```
site        depth
--------------------
site002     168
site003     160
...
```

#### `metrics_siteXXX.txt` (in metric_data folders)
```
ROI    M_S         M_C        SNR_g      SNR_p      M_S_norm   M_S_ratio
---------------------------------------------------------------------------
0      -43.738     0.51768    3.4882     0.93075    -15.137    0.66355
1      38.941      -0.087577  0.52251    0.11513    9.7861     2.4277
...
```

Columns:
- **ROI**: ROI index
- **M_S**: Metric S (Additive Curve) - measures additive response
- **M_C**: Metric C (Curve Association) - measures correlation/association
- **SNR_g**: Signal-to-noise ratio for gratings
- **SNR_p**: Signal-to-noise ratio for plaids
- **M_S_norm**: Metric S normalized (% of baseline)
- **M_S_ratio**: Metric S ratio

#### `roi_osi.txt`
```
site        roi         bsf         bdr         dsi         osi         dcv         ocv
------------------------------------------------------------------------------------------------
2           0           3.84        300         0.205091    0.980977    0.866445    0.160741
2           1           3.84        300         0.137665    0.789709    0.945801    0.389402
...
```

Columns:
- **site**: Recording site number (integer)
- **roi**: ROI index (integer)
- **bsf**: Best spatial frequency (cyc/deg)
- **bdr**: Best direction (degrees, 0-360)
- **dsi**: Direction selectivity index (0-1)
- **osi**: Orientation selectivity index - Gaku's method (0-1)
- **dcv**: Direction circular variance (0-1)
- **ocv**: Orientation circular variance (0-1)

**OSI Calculation Methods**:
- **Gaku's OSI** (column 6): `(max_response - min_response) / (max_response + min_response)`
- **Circular Variance OSI**: `1 - ocv` (uses column 8)

## Setup

### 1. Activate Virtual Environment
```bash
source .venv/bin/activate
```

### 2. Install Dependencies
```bash
pip install numpy pandas matplotlib scipy pathlib
```

## Running the Analysis

### Main Analysis Scripts

#### 1. Baseline Analysis
```bash
python scripts/bm_calc/baseline_analysis.py
```
**Analysis Type**: Log10-transformed baseline for all ROIs

**Output**: `data_baseline/log10/all_roi/` folder with:
- `metric_site/m_s/`: M_S vs log10(baseline) per site
- `metric_site/m_s_norm/`: M_S_norm vs log10(baseline) per site
- `metric_site/m_s_r/`: M_S_ratio vs log10(baseline) per site
- `metric_site/m_c/`: M_C vs log10(baseline) per site
- `pearson_site/`: R-value vs depth plots (m_s.png, m_c.png, etc.)
- `depth_site/`: Baseline vs depth plot

**Key Features**:
- Uses log10-transformed baseline values
- Individual ROI scatter (light blue) + site means (orange)
- Linear regression on site means
- Publication-quality plots with stats boxes

#### 2. Half-Width Analysis
```bash
python scripts/bm_calc/halfwidth_analysis.py
```
**Analysis Type**: Orthogonal (normalized) half-width for null ROIs only

**Output**: `data_halfwidth/orthogonal/null_roi/` folder with:
- `metric_site/m_s/`: M_S vs HW_norm per site
- `metric_site/m_s_norm/`: M_S_norm vs HW_norm per site
- `metric_site/m_s_r/`: M_S_ratio vs HW_norm per site
- `metric_site/m_c/`: M_C vs HW_norm per site
- `pearson_site/`: R-value vs depth plots
- `depth_site/`: Half-width vs depth plot

**Key Features**:
- Uses normalized half-width values only
- Filters to null ROIs (non-null values in roi_hw_orth.txt)
- Y-axis limits: 40-75 for normalized half-width
- Smart legend placement (avoids data and stats box)

#### 3. Spatial Frequency Analysis
```bash
python scripts/bm_calc/spatial_analysis.py
```
**Analysis Type**: Spatial frequency analysis for all ROIs

**Output**: `data_spatial/all_roi/` folder with:
- `metric_site/m_s/`: M_S vs SF per site
- `metric_site/m_s_norm/`: M_S_norm vs SF per site
- `metric_site/m_s_r/`: M_S_ratio vs SF per site
- `metric_site/m_c/`: M_C vs SF per site
- `pearson_site/`: R-value vs depth plots with custom titles
  - "Metric S (SF): R-value vs. Depth (N=X)"
  - "Metric C (SF): R-value vs. Depth (N=X)"
- `depth_site/`: SF vs depth plot

**Key Features**:
- Auto-scaling y-axis with intelligent padding
- Dynamic tick intervals (0.1 or 0.2 based on range)
- Stats box in lower left corner
- Smart legend placement

#### 4. Orientation Selectivity Index (OSI) Analysis
```bash
python scripts/bm_calc/osi_analysis.py
```
**Analysis Type**: OSI analysis using both circular variance and Gaku's methods for all ROIs

**Output**: `data_osi/` folder with two method subfolders:

**Method 1 - Circular Variance**: `data_osi/variance/` folder with:
- `metric_site/m_s/`: M_S vs OSI_cv per site
- `metric_site/m_s_norm/`: M_S_norm vs OSI_cv per site
- `metric_site/m_s_r/`: M_S_ratio vs OSI_cv per site
- `metric_site/m_c/`: M_C vs OSI_cv per site
- `pearson_site/`: R-value vs depth plots for OSI_cv
- `depth_site/`: OSI_cv vs depth plot

**Method 2 - Gaku's Method**: `data_osi/osi/` folder with:
- `metric_site/m_s/`: M_S vs OSI per site
- `metric_site/m_s_norm/`: M_S_norm vs OSI per site
- `metric_site/m_s_r/`: M_S_ratio vs OSI per site
- `metric_site/m_c/`: M_C vs OSI per site
- `pearson_site/`: R-value vs depth plots for OSI
- `depth_site/`: OSI vs depth plot

**Key Features**:
- Analyzes relationship between orientation selectivity and cross-orientation suppression
- Circular variance method: OSI_cv = 1 - ocv (orientation circular variance)
- Gaku's method: OSI = (max - min) / (max + min) depth of modulation
- Separate analysis for each method to compare approaches
- Custom titles indicating which OSI method is used

**Scientific Context**:
This analysis addresses whether neurons that are more selective for orientation (high OSI) show different patterns of cross-orientation suppression (M_S, M_C) compared to broadly tuned neurons (low OSI). This helps understand if suppression mechanisms depend on tuning properties.

#### 5. Comprehensive Depth Analysis
```bash
python scripts/d_calc/depth.py
```
**Output**: `depth_data/` folder with complete hierarchical analysis including:
- ROI maps (colored by metric values) for each metric at each site
- Depth vs metric plots for all filtering conditions

### Data Preparation

Before running OSI analysis, you need to prepare the OSI data file:

```bash
# Navigate to playground scripts
cd zz_Playground/Scripts/

# Run the cleaning script (requires cell_sine_dosi_docv.csv in ../Data/)
python clean.py
```

This will:
1. Read `zz_Playground/Data/cell_sine_dosi_docv.csv`
2. Clean and format the data
3. Output `raw_data/bm_data/roi_osi.txt`
4. Show summary statistics for both OSI methods

### Publication Figure Generation
```bash
python scripts/XORI\ Data\ Analysis/final_graph.py
```
**Output**: `final graphs/` folder with poster-quality figures:
- `Metric S_vs_depth_poster.png/svg`
- `Metric C_vs_depth_poster.png/svg`
- `Metric S Ratio_vs_depth_poster.png/svg`

## Analysis Workflows

### Workflow 1: Focused Analysis (Simplified)
This workflow uses the streamlined scripts for specific analyses:

1. **Baseline Analysis** (all ROIs, log10)
   ```bash
   python scripts/bm_calc/baseline_analysis.py
   ```
   - Examines relationship between baseline fluorescence and depth
   - Uses log10 transformation to handle wide range of baseline values

2. **Half-Width Analysis** (null ROIs, normalized)
   ```bash
   python scripts/bm_calc/halfwidth_analysis.py
   ```
   - Analyzes tuning width properties for null ROIs
   - Uses normalized half-width values

3. **Spatial Frequency Analysis** (all ROIs)
   ```bash
   python scripts/bm_calc/spatial_analysis.py
   ```
   - Examines spatial frequency preferences across depth

### Workflow 2: Metric-Based Analysis
1. Ensure metric files exist in `metric_data/all_roi/`
2. Run depth analysis → Generates all ROI maps and correlations
3. Generate final graphs → Publication-ready figures

### Workflow 3: SNR-Based Filtering
Compare results across different SNR filters:
- `all_roi/`: Baseline (all ROIs)
- `cull_roi/thr_cull/above_1_0/`: High-quality ROIs (SNR ≥ 1.0)
- `r_cull_roi/thr_r_cull/below_1_0/`: Low-quality ROIs (SNR < 1.0)

## Analysis Design Decisions

### Why Log10 Baseline?
- Baseline fluorescence values span orders of magnitude
- Log transformation normalizes the distribution
- Better reveals relationships across the full dynamic range

### Why Orthogonal (Normalized) Half-Width?
- Normalized values account for differences in tuning curve properties
- More consistent across different recording conditions
- Better for comparing across depths

### Why Null ROIs for Half-Width?
- Half-width data only available for subset of ROIs
- Null ROIs represent specific response type
- Prevents bias from missing data

### Why All ROIs for Baseline and Spatial Frequency?
- These measurements available for all ROIs
- Maximum statistical power
- Represents complete population

## Key Metrics Explained

### M_S (Metric S / Additive Curve)
- Measures **additive response** properties
- Positive values: Super-additive (plaid response > sum of gratings)
- Negative values: Sub-additive (plaid response < sum of gratings)
- Zero: Perfectly additive

### M_C (Metric C / Curve Association)
- Measures **correlation between grating and plaid responses**
- Range: -1 to +1 (Pearson correlation)
- High positive: Responses are similar
- Near zero: Independent responses

### M_S_norm
- M_S normalized as percentage of baseline
- Accounts for baseline fluorescence differences

### M_S_ratio
- Ratio-based version of M_S
- Often analyzed in log2 space for symmetry
- Values > 1: Super-additive
- Values < 1: Sub-additive
- Value = 1: Additive

## Output File Structure

### Per-Site Metric Plots (`metric_site/`)
For each analysis type, organized by metric:
```
m_s/
  site002_m_s_vs_baseline.png
  site003_m_s_vs_baseline.png
  ...
m_s_norm/
  site002_m_s_norm_vs_baseline.png
  ...
m_s_r/
  site002_m_s_ratio_vs_baseline.png
  ...
m_c/
  site002_m_c_vs_baseline.png
  ...
```

### Pearson Correlation Plots (`pearson_site/`)
R-values vs depth for each metric:
```
m_s.png/svg          # S vs baseline/HW/SF correlation
m_s_norm.png/svg     # S_norm correlation
m_s_r.png/svg        # S_ratio correlation
m_c.png/svg          # C correlation
```
- Filled circles: p < 0.05 (significant)
- Hollow circles: p ≥ 0.05 (non-significant)
- Dashed line at r = 0

### Depth Plots (`depth_site/`)
```
baseline_vs_depth.png/svg    # Baseline analysis
halfwidth_vs_depth.png/svg   # Half-width analysis
sf_vs_depth.png/svg          # Spatial frequency analysis
```
- Light blue: Individual ROIs
- Orange: Site means
- Black dashed: Linear regression
- White stats box: Equation, r, p, N

### ROI Maps (`roi_map/` in depth_data)
Each site gets 4 ROI maps (one per metric):
```
01_roi_map_site002.png
02_roi_map_site003.png
...
```
- Color coding: Red = positive values, Blue = negative values
- Saturation at ±2σ (standard deviations)

## Visualization Standards

All publication-quality plots include:
- **Font sizes**: 
  - 33pt for axes labels
  - 22pt for R-value plot titles/labels
  - 16pt for annotations and stats
  - 14pt for legend text
- **Markers**: 
  - Individual ROIs: 20pt light blue with alpha=0.4
  - Site means: 100pt orange circles with black edges (1.5pt width)
  - Pearson plots: 150pt blue circles
- **Line width**: 
  - Regression: 2.5pt
  - Axes spines: 1.5pt
  - Legend frame: 1.5pt
- **DPI**: 150 for PNG, vector for SVG
- **Stats annotation**: White box with black border, rounded corners

## Data Filtering Options

### SNR-Based Filtering

#### Percentage Filters
- **Top 70%**: Keep 70% of ROIs with highest SNR
- **Top 80%**: Keep 80% of ROIs with highest SNR
- **Top 90%**: Keep 90% of ROIs with highest SNR
- **Bottom 10-30%**: Keep lowest SNR ROIs (quality control)

#### Threshold Filters
- **Above 0.5**: SNR ≥ 0.5
- **Above 1.0**: SNR ≥ 1.0 (recommended for clean data)
- **Above 1.5**: SNR ≥ 1.5 (high-quality only)
- **Below thresholds**: For identifying low-quality ROIs

### Null ROI Filtering
- ROIs identified as "null" (non-responsive) in separate analysis
- Used specifically for half-width analysis
- Ensures consistent comparison across sites

## Troubleshooting

### Common Issues

**Issue**: "File not found" errors
**Solution**: Ensure you're running from the XORI root directory:
```bash
cd ~/XORI_Analysis/XORI
source .venv/bin/activate
python scripts/bm_calc/baseline_analysis.py
```

**Issue**: "No data" warnings in output
**Solution**: Check that metric files exist in `metric_data/all_roi/metrics_siteXXX.txt`

**Issue**: ROI map IndexError
**Solution**: Ensure `stat/stat_siteXXX.npy` files exist and match metric files

**Issue**: Empty plots
**Solution**: Verify `raw_data/bm_data/site_depth.txt` exists and sites match metric files

**Issue**: Baseline analysis finds no ROIs
**Solution**: Check that `raw_data/tc_data/siteXXX/` folders contain xplot files

**Issue**: Half-width analysis skips all sites
**Solution**: Verify `raw_data/bm_data/roi_hw_orth.txt` has non-null values in hw_norm column

**Issue**: Spatial frequency analysis missing data
**Solution**: Check `raw_data/bm_data/roi_sf.txt` for matching site numbers

## File Naming Conventions

- Sites: `siteXXX` where XXX is zero-padded (e.g., `site002`, `site038`)
- Metrics files: `metrics_siteXXX.txt`
- Stat files: `stat_siteXXX.npy`
- Output plots: Named by analysis type and metric
- Site numbers: Integer format in data files (2, 3, ...), converted to site names automatically

## Script Locations

### Updated Analysis Scripts (Simplified)
- `scripts/bm_calc/baseline_analysis.py` - Log10 baseline, all ROIs only
- `scripts/bm_calc/halfwidth_analysis.py` - Orthogonal HW, null ROIs only  
- `scripts/bm_calc/spatial_analysis.py` - Spatial frequency, all ROIs only

### Comprehensive Scripts
- `scripts/d_calc/depth.py` - Full depth analysis with all filter options
- `scripts/XORI Data Analysis/final_graph.py` - Publication figure generation

## Contact

Bair Lab, Department of Neurobiology & Biophysics, University of Washington

## TODO List
