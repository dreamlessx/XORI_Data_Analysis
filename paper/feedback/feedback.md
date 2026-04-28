# Feedback for Wyeth — Additional Statistical Analyses

## The Big Picture

We have a strong main finding: **P decreases with depth** (facilitation → suppression). But a reviewer will attack it from every angle. These analyses preemptively shut down every line of criticism.

Here's the logic: **one main result, four walls of defense.**

---

## 1. The Main Result

**P vs Depth** (r = -0.706, p = 2.72e-5, n = 28 sites)

- Superficial L2/3 (140-266 μm): P > 1 → plaid response exceeds linear prediction (facilitation)
- Deep L2/3 (392-518 μm): P < 1 → plaid response falls below prediction (suppression)
- Crossover at ~350 μm

**C vs Depth** (r = +0.798, p = 3.67e-7)

- Deeper neurons have MORE predictable plaid tuning shapes, even though they're more suppressed
- This is the normalization signature: gain goes down but shape is preserved

**See**: `paper/figures/fig3_depth.pdf` — Panel A (P) and Panel B (C)

---

## 2. Partial Correlations — "Is it just a confound?"

### What it does
Removes the influence of a confound variable from both P and depth, then checks if they still correlate.

### Why we need it
Reviewer says: "ROI size increases with depth. Maybe bigger ROIs just average over more neuropil and that drives lower P. Your finding is an artifact."

### The answer
After removing ALL confounds (ROI radius + SNR + OSI + SF) simultaneously:
- P-depth partial r = -0.520, p = 0.005 → **still significant**
- C-depth partial r = +0.379, p = 0.047 → **still significant**

No single confound or combination of confounds kills the effect.

### How it works (step by step)
1. Regress depth on the confound → get depth residuals (what's left of depth after removing the confound)
2. Regress P on the confound → get P residuals
3. Correlate the residuals
4. If still significant → the confound doesn't explain the P-depth relationship

### The table

| Controls removed | Partial r (P-depth) | p-value |
|---|---|---|
| None | -0.706 | 2.7e-5 |
| ROI radius | -0.751 | < 0.001 |
| All (radius + SNR + OSI + SF) | -0.520 | 0.005 |

**See**: `supplementary_analysis/outputs/partial_correlation/`

---

## 3. Mixed-Effects Models — "You only have 28 data points"

### What it does
Uses all 4,785 ROIs instead of 28 site means, while properly handling the nested structure (ROIs within sites).

### Why we need it
Reviewer says: "n = 28 sites is a small sample. How do you know this isn't driven by 2-3 outlier sites?"

### The answer
At the single-ROI level with site as a random effect:
- P ~ Depth: β = -9.30, **p = 7.9e-10**
- P ~ Depth + Radius + SNR: β = -8.64, **p = 2.2e-8** → survives controls

### How it works (plain English)
The model is: **P = baseline + depth_effect × depth + site_adjustment + noise**

- Each site gets its own baseline (the "random intercept") — so if one site is generally brighter or dimmer, that's absorbed
- The depth_effect (β) asks: across ALL sites, does increasing depth systematically decrease P?
- It's like a within-subjects design — each site is a "subject," and depth is the treatment

### Why it's better than just the site-level correlation
- Uses all 4,785 data points, not 28 averages
- Properly accounts for the fact that 170 ROIs from one site aren't 170 independent observations
- Can include covariates (radius, SNR) directly in the model

### Key detail for Wyeth
C loses significance in the mixed model after controlling for radius + SNR (p = 0.151). This means the C-depth relationship is partially confounded by morphological factors. **P survives everything.** That's why P is the stronger metric.

**See**: `supplementary_analysis/outputs/mixed_effects/mixed_effects_summary.txt`

---

## 4. Mediation Analysis — "WHY does P change with depth?"

### What it does
Identifies which covariates explain (mediate) the P-depth relationship. Turns a descriptive finding into a mechanistic story.

### Why we need it
Reviewer says: "OK the effect is real, but what drives it? Is it just that SF changes with depth?"

### The answer

| Control for | Partial r | % of P-depth effect explained |
|---|---|---|
| Nothing (baseline) | -0.706 | — |
| Spatial frequency | -0.424 | **40%** |
| Half-width (bandwidth) | -0.558 | **21%** |
| LHI | -0.642 | 9% |
| SNR | -0.671 | 5% |
| OSI | -0.678 | 4% |

### How it works
1. Start with the raw P-depth correlation: r = -0.706
2. Control for SF → r drops to -0.424
3. The drop: (0.706 - 0.424) / 0.706 = **40% reduction**
4. That 40% is the portion of the P-depth effect that "goes through" SF

### The story this tells
- SF is the biggest single mediator (40%) — makes biological sense:
  - Deeper neurons prefer lower SF
  - Lower-SF neurons have larger receptive fields
  - Larger receptive fields → broader normalization pools → stronger suppression
- Bandwidth mediates 21% — broader tuning at depth → less selective normalization
- **But 60% of the effect is NOT explained by SF** — there's something else going on at depth beyond just SF preference. Could be circuit architecture, inhibitory cell density, connectivity patterns...

### Why this matters for the paper
This is the difference between "we found a correlation" and "we found a correlation AND we can explain part of the mechanism." The unexplained 60% is what motivates the connectomics follow-up.

**See**: `supplementary_analysis/outputs/mediation_analysis/`

---

## 5. Subpopulation Splits — "Is it universal?"

### What it does
Splits neurons into groups and checks if P-depth holds in each.

### Why we need it
Reviewer says: "Maybe only low-OSI neurons drive this. Or maybe it's only a high-SF phenomenon."

### The answer
The P-depth effect holds in:
- High OSI neurons ✓
- Low OSI neurons ✓
- High SF preference neurons ✓
- Low SF preference neurons ✓

It's universal across neuron types. Not driven by one subclass.

### How it works
1. Median split all ROIs by OSI → two groups
2. Compute P-depth correlation in each group separately
3. Both groups show the effect → it's not driven by one extreme

**See**: `supplementary_analysis/outputs/subpopulation/`

---

## Summary: What to Tell Wyeth

> "The P-depth correlation is the main finding. I ran four additional analyses to make sure it's airtight:
>
> 1. **Partial correlations** — not driven by confounds (survives controlling for radius, SNR, OSI, SF)
> 2. **Mixed-effects model** — holds at the single-neuron level, p = 7.9e-10
> 3. **Mediation** — SF explains 40%, bandwidth explains 21%, but 60% is unexplained (that's the interesting part)
> 4. **Subpopulation splits** — holds for all neuron types
>
> The mediation result gives us a mechanistic story: SF partially drives the depth gradient, but there's something else — probably circuit architecture that changes with depth. That's what the connectomics will address."

---

## Where These Go in the Paper

- **Main text (Results)**: Mixed-effects confirmation + mediation (SF = 40%). These are interesting, not just controls.
- **Supplementary/Methods**: Partial correlations, subpopulation splits, bootstrap CIs. Important for rigor but not narratively interesting.
- **Discussion**: The unexplained 60% motivates the connectomics angle.

---

## Existing Figures for These Analyses

| Analysis | Figure location |
|---|---|
| Main P & C vs depth | `paper/figures/fig3_depth.pdf` |
| Partial correlations forest plot | `supplementary_analysis/outputs/partial_correlation/` |
| Mixed-effects summary | `supplementary_analysis/outputs/mixed_effects/` |
| Mediation bar chart | `supplementary_analysis/outputs/mediation_analysis/` |
| Subpopulation splits | `supplementary_analysis/outputs/subpopulation/` |
| ROI-level scatter | `supplementary_analysis/outputs/roi_scatter/` |
| Bootstrap distributions | `supplementary_analysis/outputs/bootstrap_stats/` |
| Morphology confounds | `supplementary_analysis/outputs/morphology_controls/` |
| Multi-metric depth profile | `paper/figures/fig9_depth_profile.pdf` |
| Publication summary composite | `supplementary_analysis/outputs/summary_figure/` |
