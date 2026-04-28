# XORI: cross-orientation suppression in macaque V1

Two-photon calcium imaging of cross-orientation suppression across cortical depth
in macaque V1 layer 2/3. ~4,800 ROIs, 28 fields of view, 140-518 um. First-author
manuscript in final-push phase.

**Bair Lab** (UW Neurobiology and Biophysics) in collaboration with Allen Institute
and WaNPRC.

## Vault landmark

`~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Dreamless_Machine/03-Research/BAIR-Lab/`
- `Plan.md`. Manuscript state, figure inventory, key numbers, next 3 steps. Source
  of truth for project status.
- `Notes.md`. Working notes.
- `BAIR-Lab.md`. Hub.

## Separate, unrelated project

V1 weight signatures lives at `~/v1-weight-signatures/`. That is mechanistic
interpretability work on pretrained vision nets vs V1 receptive fields. It is
NOT this paper, and it is NOT what "xori" refers to. Different paper, different
authorship, different scope, different agent. Do not conflate. Do not write to
`~/v1-weight-signatures/` from this workspace.

## Repo layout

```
~/xori/
├── raw_data/                  # ROI-level inputs, tuning curves
│   ├── bm_data/               # site_depth, roi_osi, roi_stat, roi_hw_orth, roi_lhi, roi_sf
│   └── tc_data/               # per-site tuning-curve responses
├── metric_data/               # cross-ori metrics S, R, etc, per ROI per site
│   ├── all_roi/               # primary
│   ├── cull_roi/              # high-SNR subsets
│   └── r_cull_roi/            # low-SNR subsets (QC)
├── depth_data/                # depth analysis outputs (plots, ROI maps)
├── data_baseline,_halfwidth,_lhi,_osi,_size,_spatial/
│                              # per-covariate analysis outputs
├── scripts/                   # core analysis
│   ├── m_calc/                # metric calculation from tuning curves
│   ├── d_calc/                # depth correlation analysis
│   └── bm_calc/               # covariate-metric relationships
├── stat/                      # Suite2p stat files (ROI spatial masks)
├── supplementary_analysis/    # extended stats (mixed-effects, partial corr, mediation)
└── paper/                     # manuscript + figures
    ├── manuscript.tex
    ├── references.bib
    ├── manuscript_stats.json  # canonical numbers (do NOT regenerate casually)
    ├── make_figures.py        # generates fig1, fig3..fig9
    ├── compute_stats.py       # prints LaTeX-formatted stats to stdout
    └── figures/               # PDF + PNG outputs
```

## Quick commands

All commands assume cwd = repo root unless noted.

```bash
# activate env (Python 3.12 venv with numpy/pandas/matplotlib/scipy/statsmodels)
source .venv/bin/activate

# regenerate all figures (fig1, fig3..fig9). Figure 2 not yet implemented.
python paper/make_figures.py

# print canonical stats (LaTeX-formatted, suitable for tables)
python paper/compute_stats.py

# build manuscript PDF (requires MacTeX or BasicTeX, currently NOT installed)
cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscript.tex
```

## Site038 caveat (data, not exclusion)

Site038 is the shallowest site at 140 µm. All 134 ROIs have SNR_g < 0.5. The
SNR-thresholded cull script (`scripts/m_calc/cull_metric.py:240`) correctly
logs `"no ROIs meet threshold (skipped)"` and writes no files for site038
under any of `metric_data/cull_roi/thr_cull/{above_0_5,above_1_0,above_1_5}/`.
This is documented behavior, not data corruption or analytic exclusion.

Site038 IS retained in `metric_data/all_roi/`, the percentile culls
(`per_cull/{top_70,top_80,top_90}/`), the reverse culls, and all depth and
covariate analyses. **Do not exclude site038** from any analysis without
explicit ask.

## Setup blockers

- **LaTeX is not installed.** No `latexmk`/`pdflatex`/`bibtex` on PATH. Install
  MacTeX or BasicTeX before claiming manuscript-ready. Until then,
  `paper/manuscript.pdf` is whatever the upstream repo committed.

## Key metrics (from README + manuscript)

| Code name    | Paper name | Description                                        |
|--------------|------------|----------------------------------------------------|
| `M_S_ratio`  | **S, P**   | observed/predicted plaid response ratio. Primary.  |
| `M_C`        | **R, C**   | Pearson r between observed and predicted curves.   |
| `M_S`        | (legacy)   | signed difference observed minus predicted.        |
| `M_S_norm`   | (legacy)   | signed diff normalized to baseline F.              |
| `M_X`        | secondary  | additional cross-ori metric.                       |
| `SNR_g/p`    | filter     | grating/plaid SNR.                                 |

S < 1 = suppression, S > 1 = facilitation. log2(P) used for site-level stats.

## Headline numbers

- S vs depth: r = -0.768, [-0.908, -0.617], p = 1.83e-6, n=28 sites (log2-transformed).
- R vs depth: r = +0.798, [+0.675, +0.891], p = 3.67e-7.
- Mixed-effects S~Depth at single-ROI level: coef = -9.30, p = 7.9e-10, n=4785.
- SF mediates ~40%, bandwidth ~21% of the S-depth effect.

## Final push (current focus)

1. **Figure 2 (model schematic).** Only missing figure. Methods §Modeling specifies
   Gabor pair (sf=1.0, SD_par=0.4, SD_orth=0.3) + Maxwell-difference temporal filter
   (s1=20ms, s2=40ms) + space-time separable filter outputs.
2. **Integrate fig9_depth_profile into manuscript.tex.** Open question: replace
   fig3D depth-binned panel, or stand alone? Resolve before committing.
3. **Internal review pass.** Verify every `\ref` resolves, every figure has
   `\includegraphics`, every numerical claim matches `manuscript_stats.json`,
   every reference in `references.bib` is cited.

## Operator trigger

Post-commit hook at `.git/hooks/post-commit` calls
`~/.local/lib/tri/bin/operator-trigger.sh`. Operator refreshes `BAIR-Lab/{Plan,
Notes,BAIR-Lab}.md` after each commit unless the subject contains `[no-operator]`
or `[operator]` (the latter is set when Operator itself commits, to prevent loops).

## Voice and tri usage

Inherited from `~/.claude/CLAUDE.md` and `~/.claude/VOICE.md`. For manuscript
writing and Methods/Results edits, use Register 2 (academic), modeled on Bair's
JNeurosci style. For figure captions: declarative, quantification-first, no
"interestingly" or "notably."
