# Methods

**Draft methods text for publication. Edit as needed.**

---

## Data Acquisition

Two-photon calcium imaging was performed in the primary visual cortex (V1) of anesthetized rhesus macaque (*Macaca mulatta*) expressing GCaMP6s (PHP.eB-CAG-GCaMP6s). Imaging was conducted across 28 fields of view spanning approximately 400 um of layer 2/3 depth (140-518 um below the pial surface). At each depth, neural activity was recorded from populations of neurons across wide cortical blocks. The long-term goal of this preparation is to compare population physiology with dense connectomic reconstructions from serial EM of the same tissue.

## Region of Interest (ROI) Identification

ROIs corresponding to individual neurons were identified using Suite2p (Pachitariu et al., 2017). A total of 4,785 ROIs were identified across all recording sites (mean +/- SD: 171 +/- 28 ROIs per site; range: 120-230). ROI quality metrics including pixel count, radius, aspect ratio, and signal-to-noise ratio (SNR) were computed for each cell.

## Visual Stimulation

Neurons were presented with drifting sinusoidal gratings (50% contrast; 4 Hz; 4 cyc/deg; 2 deg patch) and orthogonal plaid stimuli formed by summing two gratings at the preferred and orthogonal orientations. Receptive-field maps were obtained from flashed light/dark spots, and direction and spatial-frequency preferences were measured separately. Spatial frequency tuning was measured at 8 frequencies (0.64-7.68 cycles/degree).

## Metric Calculations

### Metric S: Suppression/Facilitation Strength

For each ROI, a linear prediction was constructed by shifting the single-grating tuning curve by -90 degrees and summing with baseline correction. Metric S is computed as the ratio of the observed plaid response to this linear prediction:

```
S = Sum(R_i^plaid) / Sum(Pred_i)
```

where R_i^plaid is the observed plaid response at direction i and Pred_i is the corresponding linear prediction. Values greater than 1 indicate facilitation (observed plaid response exceeds the linear prediction), while values less than 1 indicate suppression (observed response falls below the linear prediction). A value of 1 indicates perfect linearity.

This ratio form (`M_S_ratio` in code) normalizes for differences in baseline response magnitude across ROIs and is used for all current analyses. A legacy signed-difference version (`M_S` in code), computed as the mean of (observed - predicted) across directions, is also retained for comparison but is no longer the primary metric.

### Metric R: Shape Similarity

Metric R (`M_C` in code) is the Pearson correlation between the predicted (linear sum) and observed plaid tuning curves. Higher values indicate that the plaid tuning curve shape is more linearly predictable from the component grating responses. This metric captures shape similarity independent of overall amplitude differences captured by S.

### Signal-to-Noise Ratio
SNR was computed separately for grating (SNR_g) and plaid (SNR_p) responses as the ratio of mean response amplitude to response variability across trials.

### Orientation Selectivity Index (OSI)
OSI was computed using standard methods (Ringach et al., 2002) to quantify the degree of orientation tuning for each neuron.

### Best Spatial Frequency (BSF)
The spatial frequency eliciting the maximum response was determined from tuning curves measured at 8 spatial frequencies (0.64-7.68 cycles/degree).

### Local Homogeneity Index (LHI)
LHI quantifies the similarity of orientation preference between neighboring ROIs within each field of view, providing a measure of local functional organization.

### Orientation Tuning Half-Width
Half-width at half-height of the grating tuning curve at the orthogonal orientation, measuring the sharpness of orientation selectivity.

## Statistical Analysis

### Site-Level Analysis
For primary analyses, metrics were averaged across all ROIs within each field of view to obtain site-level means (n = 28). This approach accounts for the hierarchical structure of the data (ROIs nested within sites) and avoids pseudoreplication.

### Correlation Analysis
Pearson correlation coefficients were computed between cortical depth and metric values. Statistical significance was assessed using two-tailed tests with alpha = 0.05.

### Partial Correlations
To control for potential confounds, partial correlations were computed between depth and metrics while statistically controlling for:
- ROI radius (morphological confound)
- Signal-to-noise ratio (data quality confound)
- Orientation selectivity index (tuning property)
- Best spatial frequency (receptive field property)

Partial correlations were computed by regressing out the confound variables from both depth and metric values, then correlating the residuals.

### Mediation Analysis
To identify which covariates mediate the S-depth relationship, partial correlations were computed controlling for each covariate individually. The percent reduction in correlation strength relative to the zero-order correlation quantifies the mediating influence of each variable.

### Mixed-Effects Modeling
Linear mixed-effects models were fit using restricted maximum likelihood (REML) with recording site as a random intercept:

```
Metric ~ Depth + Covariates + (1|Site)
```

This approach properly accounts for the non-independence of ROIs within sites while leveraging the full ROI-level dataset (n = 4,785).

### Bootstrap Confidence Intervals
Non-parametric bootstrap (10,000 iterations) was used to compute 95% confidence intervals for correlation coefficients and regression slopes. Sites were resampled with replacement to maintain the hierarchical data structure.

### Subpopulation Analyses
To assess robustness, depth-metric correlations were computed separately for:
- High vs. low orientation selectivity (median split on OSI)
- High vs. low spatial frequency preference (median split on BSF)
- OSI quartiles (Q1-Q4)

### Interaction Testing
To test whether the depth effect differed across OSI groups, Fisher's z-transformation was used to compare correlation coefficients between low-OSI (Q1) and high-OSI (Q4) subpopulations.

### Within-Site Correlation Analyses
Within individual fields of view, correlations between metrics (S, R) and covariates (SF, baseline, halfwidth, LHI, OSI) were computed at the single-ROI level to test whether between-site depth trends are recapitulated within local populations.

## Control Analyses

### ROI Morphology Confound
ROI radius correlated with cortical depth (r = 0.878, p < 0.001), likely reflecting either optical factors (increased scattering at depth) or biological factors (larger cell bodies in deeper sublayers). To ensure this did not drive the main findings, partial correlations controlling for ROI radius were computed. The S-depth correlation remained highly significant after this control (r = -0.751, p < 0.001).

### Data Quality Confound
SNR varied modestly with depth (r = 0.461, p = 0.014). Partial correlations controlling for SNR confirmed that depth effects were not artifacts of data quality variations.

### Calcium-to-Fluorescence Nonlinearities
Depth-dependent changes in calcium-to-fluorescence transduction and optical factors could in principle bias metrics. The robustness of findings in both dF/F and raw fluorescence units mitigates but does not eliminate this concern.

### Modeling
Initial linear-nonlinear Gabor modeling augmented with divisive normalization and simple fluorescence transduction recapitulates the qualitative pattern: more negative S with stronger normalization and higher R with narrower tuning. This supports an interpretation of stronger cross-orientation suppression deeper in layer 2/3 rather than an optical or transduction artifact.

## Software and Code

All analyses were performed in Python 3.10+ using NumPy, Pandas, SciPy, Matplotlib, and statsmodels. ROI identification used Suite2p. Analysis code is available at https://github.com/dreamlessx/XORI_Data_Analysis.

---

## Key Statistics to Report

Note: In the codebase, S is stored as `M_S_ratio` (primary) or `M_S` (legacy), and R is stored as `M_C`.

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

### Mediation of the S-Depth Relationship
| Control Variable | Partial r | % Reduction |
|-----------------|-----------|-------------|
| None (baseline) | -0.768 | |
| + SF | -0.463 | 40% |
| + Halfwidth | -0.605 | 21% |
| + LHI | -0.697 | 9% |
| + SNR | -0.727 | 5% |
| + OSI | -0.740 | 4% |

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

2. **Correlational design**: The study establishes correlation, not causation. The depth-dependent changes could reflect layer-specific circuit properties, cell-type composition differences, or other factors that covary with depth.

3. **ROI size confound**: ROI radius correlated with depth. While partial correlations suggest the main findings are robust to this confound, optical or morphological artifacts cannot be entirely ruled out.

4. **Depth-dependent optical factors**: Calcium-to-fluorescence nonlinearities and depth-dependent optical scattering could bias metrics, though robustness in both dF/F and raw units mitigates this concern.

5. **Layer assignment**: Without histological verification, depth measurements cannot be definitively mapped to anatomical sublayers within layer 2/3.

6. **Stimulus mismatch**: Per-ROI stimulus parameters (orientation, SF) were not individually optimized for all ROIs, which may contribute variance.

7. **R-depth confound**: The R-depth correlation loses significance in the ROI-level mixed model when controlling for radius and SNR (p = 0.151), indicating partial confounding by morphological factors. The S-depth effect, by contrast, survives all controls.

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
