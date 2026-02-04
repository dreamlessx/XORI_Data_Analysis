# Supplementary Statistical Analyses

Extended analyses examining robustness, confounds, and mechanistic mediators of the depth-dependent cross-orientation interaction findings. All analyses use site-level averages (n = 28 fields of view) unless otherwise noted.

**Notation:** S refers to `M_S_ratio` (primary metric) or `M_S` (legacy signed difference) in code; R refers to `M_C` in code. See the main [README](../README.md) for metric definitions.

---

## 1. Subpopulation Analysis

**Output:** `outputs/subpopulation/`

| Subgroup | N ROIs | S ~ Depth (r) | p-value | R ~ Depth (r) | p-value |
|----------|--------|---------------|---------|---------------|---------|
| All ROIs | 4,785 | -0.768 | 1.83e-06 | +0.798 | 3.67e-07 |
| High OSI | 2,393 | -0.769 | 1.72e-06 | +0.868 | 2.18e-09 |
| Low OSI | 2,392 | -0.761 | 2.55e-06 | +0.684 | 6.03e-05 |
| High SF | 3,738 | -0.752 | 3.95e-06 | +0.832 | 4.08e-08 |
| Low SF | 1,047 | -0.752 | 3.96e-06 | +0.534 | 3.44e-03 |

The S-depth effect is consistent across all subgroups (r ~ -0.75 to -0.77). The R-depth effect is stronger in high-OSI cells (r = +0.87 vs +0.68) and high-SF cells (r = +0.83 vs +0.53).

---

## 2. Partial Correlation Analysis

**Output:** `outputs/partial_correlation/`

| Control Variables | S ~ Depth (r) | p-value | R ~ Depth (r) | p-value |
|-------------------|---------------|---------|---------------|---------|
| None | -0.768 | 1.83e-06 | +0.798 | 3.67e-07 |
| SNR | -0.727 | 1.20e-05 | +0.804 | 2.57e-07 |
| OSI | -0.740 | 6.90e-06 | +0.768 | 1.83e-06 |
| SF | -0.473 | 1.10e-02 | +0.632 | 3.13e-04 |
| All (SNR + OSI + SF) | -0.406 | 3.20e-02 | +0.652 | 1.73e-04 |

Both correlations remain significant after all controls. SF preference produces the largest reduction in the S-depth correlation (-0.77 to -0.47), consistent with partial mediation.

---

## 3. Within-Site Variance

**Output:** `outputs/within_site/`

| Statistic | Value |
|-----------|-------|
| Mean ROIs per site | 170.9 |
| Range | 120 - 230 |
| Mean within-site SD (S) | 23.12 |
| Mean within-site SD (R) | 0.327 |

Substantial within-site variance is expected at the single-neuron level. The between-site depth trend remains clearly visible in error bar plots.

---

## 4. S-R Relationship

**Output:** `outputs/ms_mc_relationship/`

| Level | r | p-value |
|-------|---|---------|
| Site-level | -0.647 | 1.98e-04 |
| ROI-level | -0.110 | 2.52e-14 |

S and R are anticorrelated: as S decreases with depth (more suppressive), R increases (plaid tuning becomes more linearly predictable). This is consistent with a transition from facilitative to normalized cross-orientation responses across layer 2/3.

---

## 5. ROI Morphology Confound

**Output:** `outputs/morphology_controls/`

| Comparison | r | p-value |
|------------|---|---------|
| Depth vs ROI radius | +0.878 | < 0.001 |
| Depth vs ROI npix | +0.856 | < 0.001 |
| Depth vs SNR | +0.461 | 0.014 |
| S vs ROI radius | -0.526 | 0.004 |
| R vs ROI radius | +0.725 | < 0.001 |

ROI radius increases with imaging depth (likely reflecting optical scattering and/or larger cell bodies at depth). This creates a potential confound chain: Depth -> ROI size -> Metrics. Addressed directly in Analysis 8 below.

---

## 6. Bootstrap Statistics

**Output:** `outputs/bootstrap_stats/`

10,000 bootstrap iterations, resampling sites with replacement.

| Metric | Observed | 95% CI | Excludes zero |
|--------|----------|--------|---------------|
| S ~ Depth (r) | -0.768 | [-0.908, -0.617] | Yes |
| S slope | -0.0818 | [-0.1131, -0.0566] | Yes |
| R ~ Depth (r) | +0.798 | [+0.675, +0.891] | Yes |
| R slope | +0.00076 | [+0.00055, +0.00099] | Yes |

---

## 7. ROI-Level Scatter

**Output:** `outputs/roi_scatter/`

| Level | S ~ Depth (r) | R ~ Depth (r) |
|-------|---------------|---------------|
| Site means (n = 28) | -0.768 | +0.798 |
| All ROIs (n = 4,785) | -0.344 | +0.246 |

ROI-level correlations are weaker (more noise) but highly significant (p < 10^-60). Site averaging provides stable estimates without inflating the effect.

---

## 8. ROI Size Confound Control

**Output:** `outputs/roi_size_control/`

Partial correlations controlling for ROI morphology.

| Control | S ~ Depth (r) | p-value | R ~ Depth (r) | p-value |
|---------|---------------|---------|---------------|---------|
| None | -0.768 | < 0.001 | +0.798 | < 0.001 |
| + Radius | **-0.751** | < 0.001 | +0.489 | 0.008 |
| + Radius, SNR | **-0.736** | < 0.001 | +0.539 | 0.003 |
| + All confounds | **-0.520** | 0.005 | +0.379 | 0.047 |

The S-depth correlation is minimally affected by ROI radius control (-0.77 to -0.75). The R-depth correlation is more affected (+0.80 to +0.49) but remains significant.

---

## 9. Mixed-Effects Models

**Output:** `outputs/mixed_effects/`

Linear mixed-effects models with random intercept for site, fit at the ROI level (n = 4,785).

| Model | Depth coefficient | p-value |
|-------|------------------|---------|
| S ~ Depth | -9.30 | 7.9 x 10^-10 |
| S ~ Depth + Radius + SNR | -8.64 | 2.2 x 10^-8 |
| R ~ Depth | +0.086 | < 0.001 |
| R ~ Depth + Radius + SNR | +0.036 | 0.151 (NS) |

The S-depth effect survives all controls at the ROI level. The R-depth effect loses significance when controlling for radius and SNR in the mixed model, indicating partial confounding.

---

## 10. OSI x Depth Interaction

**Output:** `outputs/osi_interaction/`

| OSI Quartile | OSI Range | S ~ Depth (r) | R ~ Depth (r) |
|--------------|-----------|---------------|---------------|
| Q1 (Low) | 0.02 - 0.44 | -0.734*** | +0.455* |
| Q2 | 0.44 - 0.66 | -0.719*** | +0.733*** |
| Q3 | 0.66 - 0.84 | -0.759*** | +0.887*** |
| Q4 (High) | 0.84 - 1.00 | -0.729*** | +0.779*** |

Fisher's z-test comparing Q1 vs Q4 for R: z = 1.95, p = 0.051 (not significant). The S-depth effect is stable across all OSI levels.

---

## 11. Publication-Ready Summary Figure

**Output:** `outputs/summary_figure/`

Four-panel figure (PNG 300 dpi, PDF, SVG):
- **A:** S vs Depth with SEM error bars
- **B:** R vs Depth with SEM error bars
- **C:** S vs R colored by depth
- **D:** Depth-binned summary (superficial / middle / deep)

---

## 12. Halfwidth (Tuning Width) vs Depth

**Output:** `outputs/halfwidth_analysis/`

| Correlation | r | p-value |
|-------------|---|---------|
| HW_raw ~ Depth | +0.818 | < 0.0001 |
| HW_norm ~ Depth | +0.869 | < 0.0001 |
| HW ~ S | -0.597 | 0.0008 |
| HW ~ R | +0.652 | 0.0002 |

Orientation tuning width increases with depth. Neurons with broader tuning show stronger suppression (lower S) and higher shape similarity (higher R).

---

## 13. Local Homogeneity Index (LHI) vs Depth

**Output:** `outputs/lhi_analysis/`

| Correlation | r | p-value |
|-------------|---|---------|
| LHI2 ~ Depth | -0.634 | 0.0003 |
| LHI3 ~ Depth | -0.927 | < 0.0001 |
| LHI2 ~ S | +0.453 | 0.0155 |
| LHI2 ~ R | -0.578 | 0.0013 |

Local homogeneity decreases with depth, correlating with metrics in the expected direction.

---

## 14. Spatial Frequency Preference vs Depth

**Output:** `outputs/sf_analysis/`

| Correlation | r | p-value |
|-------------|---|---------|
| SF ~ Depth | -0.792 | < 0.0001 |
| SF ~ S | +0.723 | < 0.0001 |
| SF ~ R | -0.631 | 0.0003 |
| OSI ~ Depth | -0.490 | 0.008 |

Preferred SF decreases with depth. Superficial high-SF neurons show facilitation; deep low-SF neurons show suppression.

---

## 15. Multi-Metric Correlation Matrix

**Output:** `outputs/correlation_matrix/`

All metrics correlating significantly with depth (site-level):

| Metric | r with Depth | p-value | Direction |
|--------|-------------|---------|-----------|
| S | -0.768 | < 0.0001 | Decreases |
| R | +0.798 | < 0.0001 | Increases |
| M_X | +0.443 | 0.018 | Increases |
| SNR | +0.461 | 0.014 | Increases |
| Halfwidth | +0.818 | < 0.0001 | Increases |
| LHI | -0.634 | 0.0003 | Decreases |
| SF | -0.792 | < 0.0001 | Decreases |
| OSI | -0.490 | 0.008 | Decreases |

A coherent laminar profile: superficial layer 2/3 is characterized by high SF preference, narrow tuning, high OSI, and facilitation; deeper layer 2/3 by low SF, broad tuning, low OSI, and suppression.

---

## 16. Mediation Analysis

**Output:** `outputs/mediation_analysis/`

Partial correlations identifying which covariates mediate the S-Depth relationship.

| Control Variable | Partial r | % Reduction |
|-----------------|-----------|-------------|
| None (baseline) | -0.768 | |
| + SF | -0.463 | 40% |
| + Halfwidth | -0.605 | 21% |
| + LHI | -0.697 | 9% |
| + SNR | -0.727 | 5% |
| + OSI | -0.740 | 4% |

SF preference is the largest single mediator (40%), followed by tuning width (21%). The S-depth effect remains significant (p < 0.05) after controlling for all covariates.

---

## 17. Comprehensive Depth Profile

**Output:** `outputs/depth_profile/`

Four-panel summary figure:
- **A:** All metrics vs depth (z-scored overlay)
- **B:** Forest plot of all depth correlations with 95% CIs
- **C:** Metrics by depth bin (superficial / middle / deep)
- **D:** Summary statistics table

---

## File Structure

```
supplementary_analysis/
├── README.md                    # This file
├── METHODS.md                   # Draft methods section for publication
├── scripts/
│   ├── run_all_analyses.py      # Analyses 1-7
│   ├── additional_analyses.py   # Analyses 8-11
│   └── extended_analyses.py     # Analyses 12-17
└── outputs/
    ├── subpopulation/           # 1. OSI/SF subgroup splits
    ├── partial_correlation/     # 2. Partial correlation forest plot
    ├── within_site/             # 3. Error bar plots
    ├── ms_mc_relationship/      # 4. S vs R scatter
    ├── morphology_controls/     # 5. ROI size confound
    ├── bootstrap_stats/         # 6. Bootstrap distributions
    ├── roi_scatter/             # 7. ROI-level scatter
    ├── roi_size_control/        # 8. Radius-controlled partials
    ├── mixed_effects/           # 9. Mixed-effects model output
    ├── osi_interaction/         # 10. OSI quartile analysis
    ├── summary_figure/          # 11. Publication figure
    ├── halfwidth_analysis/      # 12. Tuning width vs depth
    ├── lhi_analysis/            # 13. LHI vs depth
    ├── sf_analysis/             # 14. SF preference vs depth
    ├── correlation_matrix/      # 15. Multi-metric heatmap
    ├── depth_profile/           # 16. Comprehensive depth profile
    └── mediation_analysis/      # 17. Mediation forest plot
```

---

## Summary of Key Results

1. **S-depth effect is robust**: Survives all partial correlation controls, mixed-effects modeling with covariates, and subpopulation splits. Consistent across OSI quartiles.
2. **R-depth effect is partially confounded**: Significant at the site level but loses significance in the ROI-level mixed model when controlling for radius and SNR. Should be interpreted with caution.
3. **SF preference is the primary mediator**: Accounts for ~40% of the S-depth relationship. Tuning width accounts for an additional ~21%.
4. **Coherent laminar profile**: Eight metrics show coordinated depth-dependent changes consistent with a systematic transformation of cross-orientation processing across layer 2/3.
