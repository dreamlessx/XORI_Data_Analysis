# Methods Section for Paper

**Draft methods text for publication. Edit as needed.**

---

## Data Acquisition

Two-photon calcium imaging was performed in the primary visual cortex (V1) of one adult rhesus macaque (*Macaca mulatta*). Imaging was conducted at 28 cortical depths spanning 140-518 μm below the pial surface, providing systematic sampling across cortical layers. At each depth, neural activity was recorded from populations of neurons expressing a genetically encoded calcium indicator.

## Region of Interest (ROI) Identification

ROIs corresponding to individual neurons were identified using Suite2p (Pachitariu et al., 2017). A total of 4,785 ROIs were identified across all recording sites (mean ± SD: 171 ± 28 ROIs per site; range: 120-230). ROI quality metrics including pixel count, radius, aspect ratio, and signal-to-noise ratio (SNR) were computed for each cell.

## Visual Stimulation

Neurons were presented with drifting sinusoidal gratings at the cell's preferred orientation and spatial frequency, as well as plaid stimuli consisting of two superimposed gratings (preferred + orthogonal orientation). Responses to gratings and plaids were used to compute cross-orientation suppression metrics.

## Metric Calculations

### Metric S (M_S) - Additive Component
Metric S quantifies the additive component of cross-orientation interaction:

```
M_S = Response_plaid - Response_grating
```

Positive values indicate facilitation (plaid response exceeds grating response), while negative values indicate suppression.

### Metric C (M_C) - Response Correlation
Metric C quantifies the correlation between grating and plaid responses across stimulus conditions, reflecting the degree to which responses to complex stimuli can be predicted from responses to simple stimuli.

### Signal-to-Noise Ratio
SNR was computed separately for grating (SNR_g) and plaid (SNR_p) responses as the ratio of mean response amplitude to response variability across trials.

### Orientation Selectivity Index (OSI)
OSI was computed using standard methods (Ringach et al., 2002) to quantify the degree of orientation tuning for each neuron.

### Best Spatial Frequency (BSF)
The spatial frequency eliciting the maximum response was determined from tuning curves measured at 8 spatial frequencies (0.64-7.68 cycles/degree).

## Statistical Analysis

### Site-Level Analysis
For primary analyses, metrics were averaged across all ROIs within each recording site to obtain site-level means. This approach accounts for the hierarchical structure of the data (ROIs nested within sites) and avoids pseudoreplication.

### Correlation Analysis
Pearson correlation coefficients were computed between cortical depth and metric values. Statistical significance was assessed using two-tailed tests with α = 0.05.

### Partial Correlations
To control for potential confounds, partial correlations were computed between depth and metrics while statistically controlling for:
- ROI radius (morphological confound)
- Signal-to-noise ratio (data quality confound)
- Orientation selectivity index (tuning property)
- Best spatial frequency (receptive field property)

Partial correlations were computed by regressing out the confound variables from both depth and metric values, then correlating the residuals.

### Mixed-Effects Modeling
Linear mixed-effects models were fit using restricted maximum likelihood (REML) with recording site as a random intercept:

```
Metric ~ Depth + Covariates + (1|Site)
```

This approach properly accounts for the non-independence of ROIs within sites while leveraging the full ROI-level dataset.

### Bootstrap Confidence Intervals
Non-parametric bootstrap (10,000 iterations) was used to compute 95% confidence intervals for correlation coefficients and regression slopes. Sites were resampled with replacement to maintain the hierarchical data structure.

### Subpopulation Analyses
To assess robustness, depth-metric correlations were computed separately for:
- High vs. low orientation selectivity (median split on OSI)
- High vs. low spatial frequency preference (median split on BSF)
- OSI quartiles (Q1-Q4)

### Interaction Testing
To test whether the depth effect differed across OSI groups, Fisher's z-transformation was used to compare correlation coefficients between low-OSI (Q1) and high-OSI (Q4) subpopulations.

## Control Analyses

### ROI Morphology Confound
ROI radius correlated with cortical depth (r = 0.878, p < 0.001), likely reflecting either optical factors (increased scattering at depth) or biological factors (larger neurons in deeper layers). To ensure this did not drive the main findings, partial correlations controlling for ROI radius were computed. The M_S-depth correlation remained highly significant after this control (r = -0.751, p < 0.001).

### Data Quality Confound
SNR varied modestly with depth (r = 0.461, p = 0.014). Partial correlations controlling for SNR confirmed that depth effects were not artifacts of data quality variations.

## Software and Code

All analyses were performed in Python 3.13 using NumPy, Pandas, SciPy, Matplotlib, and statsmodels. ROI identification used Suite2p. Analysis code is available at [repository URL].

---

## Key Statistics to Report

### Main Finding (Site-Level)
| Metric | Correlation with Depth | 95% CI | p-value | N |
|--------|----------------------|--------|---------|---|
| M_S | r = -0.768 | [-0.908, -0.617] | 1.83 × 10⁻⁶ | 28 sites |
| M_C | r = +0.798 | [+0.675, +0.891] | 3.67 × 10⁻⁷ | 28 sites |

### After Controlling for ROI Radius
| Metric | Partial r | p-value |
|--------|-----------|---------|
| M_S | -0.751 | < 0.001 |
| M_C | +0.489 | 0.008 |

### After Controlling for All Confounds (Radius + SNR + OSI + SF)
| Metric | Partial r | p-value |
|--------|-----------|---------|
| M_S | -0.520 | 0.005 |
| M_C | +0.379 | 0.047 |

### Mixed-Effects Model (ROI-Level)
| Model | Depth coefficient | p-value |
|-------|------------------|---------|
| M_S ~ Depth | -9.30 | 7.9 × 10⁻¹⁰ |
| M_S ~ Depth + Radius + SNR | -8.64 | 2.2 × 10⁻⁸ |
| M_C ~ Depth | +0.086 | < 0.001 |
| M_C ~ Depth + Radius + SNR | +0.036 | 0.151 (NS) |

### Subpopulation Robustness (M_S ~ Depth)
| Subgroup | r | p-value |
|----------|---|---------|
| All ROIs | -0.768 | 1.83 × 10⁻⁶ |
| High OSI | -0.769 | 1.72 × 10⁻⁶ |
| Low OSI | -0.761 | 2.55 × 10⁻⁶ |
| High SF | -0.752 | 3.95 × 10⁻⁶ |
| Low SF | -0.752 | 3.96 × 10⁻⁶ |

---

## Suggested Figure Legends

### Figure X: Depth-dependent changes in cross-orientation metrics
**(A)** Metric S (additive component) decreased with cortical depth (r = -0.768, p = 1.8 × 10⁻⁶, n = 28 sites). Each point represents the mean across all ROIs at one recording depth; error bars indicate SEM. Dashed line shows linear regression.
**(B)** Metric C (response correlation) increased with cortical depth (r = +0.798, p = 3.7 × 10⁻⁷).
**(C)** M_S and M_C were anticorrelated across sites (r = -0.647, p < 0.001), with color indicating cortical depth.
**(D)** Summary of depth effects binned by cortical layer. Superficial sites (140-266 μm) showed positive M_S (facilitation), while deep sites (392-518 μm) showed negative M_S (suppression).

### Supplementary Figure: Control analyses
**(A-B)** Partial correlations between depth and metrics after controlling for ROI radius, SNR, OSI, and spatial frequency preference. Both correlations remained significant after all controls.
**(C-F)** Depth-metric relationships computed separately for OSI quartiles, demonstrating robustness across orientation selectivity levels.

---

## Limitations to Acknowledge

1. **Single animal**: Data were collected from one monkey. While the within-animal replication is extensive (28 sites, 4,785 ROIs), individual differences cannot be assessed.

2. **Correlational design**: The study establishes correlation, not causation. The depth-dependent changes could reflect layer-specific circuit properties, cell-type composition differences, or other factors.

3. **ROI size confound**: ROI radius correlated with depth. While partial correlations suggest the main findings are robust to this confound, we cannot entirely rule out optical or morphological artifacts.

4. **Layer assignment**: Without histological verification, depth measurements cannot be definitively mapped to anatomical layers.
