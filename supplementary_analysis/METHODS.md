# Methods Section for Paper

**Draft methods text for publication. Edit as needed.**

---

## Data Acquisition

Two-photon calcium imaging was performed in the primary visual cortex (V1) of one adult rhesus macaque (*Macaca mulatta*) expressing GCaMP6s (PHP.eB-CAG-GCaMP6s), under anesthesia. Imaging was conducted across 28 fields of view spanning approximately 400 μm of layer 2/3 depth (140-518 μm below the pial surface). At each depth, neural activity was recorded from populations of neurons across wide cortical blocks.

## Region of Interest (ROI) Identification

ROIs corresponding to individual neurons were identified using Suite2p (Pachitariu et al., 2017). A total of 4,785 ROIs were identified across all recording sites (mean +/- SD: 171 +/- 28 ROIs per site; range: 120-230). ROI quality metrics including pixel count, radius, aspect ratio, and signal-to-noise ratio (SNR) were computed for each cell.

## Visual Stimulation

Neurons were presented with drifting sinusoidal gratings (50% contrast; 4 Hz; 4 cyc/deg; 2 degree patch) and orthogonal plaid stimuli formed by summing two gratings at the preferred and orthogonal orientations. Receptive-field maps were obtained from flashed light/dark spots, and direction and spatial-frequency preferences were measured separately.

## Metric Calculations

### Metric S - Suppression/Facilitation Strength

For each ROI, a linear prediction was constructed by shifting the single-grating tuning curve by -90 degrees and summing with baseline correction. Metric S is computed as the ratio of the observed plaid response to this linear prediction:

```
S = Response_observed_plaid / Response_linear_prediction
```

Values greater than 1 indicate facilitation (observed plaid response exceeds the linear prediction), while values less than 1 indicate suppression (observed response falls below the linear prediction). A value of 1 indicates perfect linearity.

Note: A legacy signed-difference version of S (`M_S` in code) is also computed but is no longer the primary metric. The ratio form (`M_S_ratio` in code) is used for all current analyses because it normalizes for differences in baseline response magnitude across ROIs.

### Metric R - Shape Similarity

Metric R is the Pearson correlation between the predicted (linear sum) and observed plaid tuning curves. Higher values indicate that the plaid tuning curve shape is more linearly predictable from the component grating responses. This metric captures shape similarity independent of overall amplitude differences captured by S.

### Signal-to-Noise Ratio
SNR was computed separately for grating (SNR_g) and plaid (SNR_p) responses as the ratio of mean response amplitude to response variability across trials.

### Orientation Selectivity Index (OSI)
OSI was computed using standard methods (Ringach et al., 2002) to quantify the degree of orientation tuning for each neuron.

### Best Spatial Frequency (BSF)
The spatial frequency eliciting the maximum response was determined from tuning curves measured at 8 spatial frequencies (0.64-7.68 cycles/degree).

## Statistical Analysis

### Site-Level Analysis
For primary analyses, metrics were averaged across all ROIs within each field of view to obtain site-level means. This approach accounts for the hierarchical structure of the data (ROIs nested within sites) and avoids pseudoreplication.

### Correlation Analysis
Pearson correlation coefficients were computed between cortical depth and metric values. Statistical significance was assessed using two-tailed tests with alpha = 0.05.

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
ROI radius correlated with cortical depth (r = 0.878, p < 0.001), likely reflecting either optical factors (increased scattering at depth) or biological factors (larger neurons in deeper layers). To ensure this did not drive the main findings, partial correlations controlling for ROI radius were computed. The S-depth correlation remained highly significant after this control (r = -0.751, p < 0.001).

### Data Quality Confound
SNR varied modestly with depth (r = 0.461, p = 0.014). Partial correlations controlling for SNR confirmed that depth effects were not artifacts of data quality variations.

### Ca2+-to-Fluorescence Nonlinearities
Depth-dependent changes in Ca2+-to-fluorescence transduction and optical factors could in principle bias metrics. The robustness of findings in both dF/F and raw fluorescence units mitigates but does not eliminate this concern.

## Software and Code

All analyses were performed in Python 3.13 using NumPy, Pandas, SciPy, Matplotlib, and statsmodels. ROI identification used Suite2p. Analysis code is available at https://github.com/dreamlessx/XORI_Data_Analysis.

---

## Key Statistics to Report

Note: In the codebase, S is stored as `M_S` and R is stored as `M_C`.

### Main Finding (Site-Level)
| Metric | Correlation with Depth | 95% CI | p-value | N |
|--------|----------------------|--------|---------|---|
| S | r = -0.768 | [-0.908, -0.617] | 1.83 x 10^-6 | 28 sites |
| R | r = +0.798 | [+0.675, +0.891] | 3.67 x 10^-7 | 28 sites |

### After Controlling for ROI Radius
| Metric | Partial r | p-value |
|--------|-----------|---------|
| S | -0.751 | < 0.001 |
| R | +0.489 | 0.008 |

### After Controlling for All Confounds (Radius + SNR + OSI + SF)
| Metric | Partial r | p-value |
|--------|-----------|---------|
| S | -0.520 | 0.005 |
| R | +0.379 | 0.047 |

### Mixed-Effects Model (ROI-Level)
| Model | Depth coefficient | p-value |
|-------|------------------|---------|
| S ~ Depth | -9.30 | 7.9 x 10^-10 |
| S ~ Depth + Radius + SNR | -8.64 | 2.2 x 10^-8 |
| R ~ Depth | +0.086 | < 0.001 |
| R ~ Depth + Radius + SNR | +0.036 | 0.151 (NS) |

### Subpopulation Robustness (S ~ Depth)
| Subgroup | r | p-value |
|----------|---|---------|
| All ROIs | -0.768 | 1.83 x 10^-6 |
| High OSI | -0.769 | 1.72 x 10^-6 |
| Low OSI | -0.761 | 2.55 x 10^-6 |
| High SF | -0.752 | 3.95 x 10^-6 |
| Low SF | -0.752 | 3.96 x 10^-6 |

---

## Suggested Figure Legends

### Figure X: Depth-dependent changes in cross-orientation metrics
**(A)** S (suppression/facilitation strength) decreased with cortical depth (r = -0.768, p = 1.8 x 10^-6, n = 28 fields of view). Each point represents the mean across all ROIs at one imaging depth; error bars indicate SEM. Dashed line shows linear regression. Positive S indicates facilitation; negative S indicates suppression relative to the linear prediction.
**(B)** R (shape similarity) increased with cortical depth (r = +0.798, p = 3.7 x 10^-7). Higher R indicates that plaid tuning is more linearly predictable from component grating responses.
**(C)** S and R were anticorrelated across sites (r = -0.647, p < 0.001), with color indicating cortical depth.
**(D)** Summary of depth effects binned by layer 2/3 depth. Superficial sites (140-266 um) showed positive S (facilitation), while deep sites (392-518 um) showed negative S (suppression).

### Supplementary Figure: Control analyses
**(A-B)** Partial correlations between depth and metrics after controlling for ROI radius, SNR, OSI, and spatial frequency preference. Both correlations remained significant after all controls.
**(C-F)** Depth-metric relationships computed separately for OSI quartiles, demonstrating robustness across orientation selectivity levels.

---

## Limitations to Acknowledge

1. **Single animal**: Data were collected from one macaque. While the within-animal replication is extensive (28 fields of view, 4,785 ROIs), individual differences cannot be assessed.

2. **Correlational design**: The study establishes correlation, not causation. The depth-dependent changes could reflect layer-specific circuit properties, cell-type composition differences, or other factors.

3. **ROI size confound**: ROI radius correlated with depth. While partial correlations suggest the main findings are robust to this confound, we cannot entirely rule out optical or morphological artifacts.

4. **Depth-dependent optical factors**: Ca2+-to-fluorescence nonlinearities and depth-dependent optical scattering could bias metrics, though robustness in both dF/F and raw units mitigates this concern.

5. **Layer assignment**: Without histological verification, depth measurements cannot be definitively mapped to anatomical sublayers within layer 2/3.

6. **Stimulus mismatch**: Per-ROI stimulus parameters (orientation, SF) were not individually optimized for all ROIs, which may contribute variance.
