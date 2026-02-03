#!/usr/bin/env python3
"""
Additional Analyses for XORI Project
- ROI size confound controls
- Mixed-effects modeling
- OSI interaction deep-dive
- Publication-ready summary figure
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress, ttest_ind
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try to import statsmodels for mixed effects
try:
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("[WARNING] statsmodels not installed - skipping mixed-effects analysis")

# Set up paths
BASE_DIR = Path("/Users/muditagar/XORI_Analysis/XORI")
RAW_DATA = BASE_DIR / "raw_data" / "bm_data"
METRIC_DATA = BASE_DIR / "metric_data" / "all_roi"
OUTPUT_DIR = BASE_DIR / "supplementary_analysis" / "outputs"

# ============================================================
# DATA LOADING (same as before)
# ============================================================

def load_all_data():
    """Load and merge all data"""
    # Site depths
    depths = {}
    with open(RAW_DATA / "site_depth.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                depths[parts[0]] = float(parts[1])

    # OSI data
    osi_data = []
    with open(RAW_DATA / "roi_osi.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 8:
                try:
                    osi_data.append({
                        'site': int(parts[0]), 'roi': int(parts[1]),
                        'bsf': float(parts[2]), 'osi': float(parts[5])
                    })
                except: pass
    osi_df = pd.DataFrame(osi_data)

    # Stat data
    stat_data = []
    with open(RAW_DATA / "roi_stat.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 14:
                try:
                    stat_data.append({
                        'site': int(parts[0]), 'roi': int(parts[1]),
                        'npix': int(parts[2]), 'radius': float(parts[6])
                    })
                except: pass
    stat_df = pd.DataFrame(stat_data)

    # Metric data
    metric_data = []
    for mf in METRIC_DATA.glob("metrics_site*.txt"):
        site_name = mf.stem.replace("metrics_", "")
        site_num = int(site_name.replace("site", ""))
        if site_name not in depths: continue
        depth = depths[site_name]

        with open(mf, 'r') as f:
            for line in f.readlines()[2:]:
                parts = line.strip().split()
                if len(parts) >= 8:
                    try:
                        metric_data.append({
                            'site': site_num, 'site_name': site_name,
                            'roi': int(parts[0]), 'depth': depth,
                            'm_s': float(parts[1]), 'm_c': float(parts[2]),
                            'snr_g': float(parts[3]), 'snr_p': float(parts[4])
                        })
                    except: pass
    metric_df = pd.DataFrame(metric_data)

    # Merge
    df = metric_df.merge(osi_df, on=['site', 'roi'], how='left')
    df = df.merge(stat_df, on=['site', 'roi'], how='left')

    return df, depths

# ============================================================
# ANALYSIS 1: ROI SIZE CONFOUND CONTROL
# ============================================================

def run_roi_size_control(df):
    """Control for ROI size in depth correlations"""
    print("\n" + "="*60)
    print("ANALYSIS: ROI SIZE CONFOUND CONTROL")
    print("="*60)

    out_dir = OUTPUT_DIR / "roi_size_control"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'radius': 'mean',
        'npix': 'mean',
        'snr_g': 'mean',
        'osi': 'mean',
        'bsf': 'mean'
    }).reset_index().dropna()

    # Partial correlation function
    def partial_corr(x, y, covariates):
        from numpy.linalg import lstsq
        X_cov = np.column_stack([np.ones(len(covariates)), covariates])
        beta_x, _, _, _ = lstsq(X_cov, x, rcond=None)
        x_resid = x - X_cov @ beta_x
        beta_y, _, _, _ = lstsq(X_cov, y, rcond=None)
        y_resid = y - X_cov @ beta_y
        return pearsonr(x_resid, y_resid)

    results = {}

    # Simple correlations
    r_ms, p_ms = pearsonr(site_data['depth'], site_data['m_s'])
    r_mc, p_mc = pearsonr(site_data['depth'], site_data['m_c'])
    results['simple'] = {'r_ms': r_ms, 'p_ms': p_ms, 'r_mc': r_mc, 'p_mc': p_mc}

    print("\n--- Simple Correlations ---")
    print(f"M_S ~ Depth: r = {r_ms:+.3f}, p = {p_ms:.4f}")
    print(f"M_C ~ Depth: r = {r_mc:+.3f}, p = {p_mc:.4f}")

    # Control for radius only
    cov_radius = site_data['radius'].values.reshape(-1, 1)
    r_ms_rad, p_ms_rad = partial_corr(site_data['depth'].values, site_data['m_s'].values, cov_radius)
    r_mc_rad, p_mc_rad = partial_corr(site_data['depth'].values, site_data['m_c'].values, cov_radius)
    results['control_radius'] = {'r_ms': r_ms_rad, 'p_ms': p_ms_rad, 'r_mc': r_mc_rad, 'p_mc': p_mc_rad}

    print("\n--- Controlling for ROI Radius ---")
    print(f"M_S ~ Depth | Radius: r = {r_ms_rad:+.3f}, p = {p_ms_rad:.4f}")
    print(f"M_C ~ Depth | Radius: r = {r_mc_rad:+.3f}, p = {p_mc_rad:.4f}")

    # Control for radius + SNR
    cov_rad_snr = np.column_stack([site_data['radius'], site_data['snr_g']])
    r_ms_rs, p_ms_rs = partial_corr(site_data['depth'].values, site_data['m_s'].values, cov_rad_snr)
    r_mc_rs, p_mc_rs = partial_corr(site_data['depth'].values, site_data['m_c'].values, cov_rad_snr)
    results['control_radius_snr'] = {'r_ms': r_ms_rs, 'p_ms': p_ms_rs, 'r_mc': r_mc_rs, 'p_mc': p_mc_rs}

    print("\n--- Controlling for ROI Radius + SNR ---")
    print(f"M_S ~ Depth | Radius+SNR: r = {r_ms_rs:+.3f}, p = {p_ms_rs:.4f}")
    print(f"M_C ~ Depth | Radius+SNR: r = {r_mc_rs:+.3f}, p = {p_mc_rs:.4f}")

    # Control for ALL confounds
    cov_all = np.column_stack([site_data['radius'], site_data['snr_g'], site_data['osi'], site_data['bsf']])
    r_ms_all, p_ms_all = partial_corr(site_data['depth'].values, site_data['m_s'].values, cov_all)
    r_mc_all, p_mc_all = partial_corr(site_data['depth'].values, site_data['m_c'].values, cov_all)
    results['control_all'] = {'r_ms': r_ms_all, 'p_ms': p_ms_all, 'r_mc': r_mc_all, 'p_mc': p_mc_all}

    print("\n--- Controlling for Radius+SNR+OSI+SF ---")
    print(f"M_S ~ Depth | All: r = {r_ms_all:+.3f}, p = {p_ms_all:.4f}")
    print(f"M_C ~ Depth | All: r = {r_mc_all:+.3f}, p = {p_mc_all:.4f}")

    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    labels = ['No control', '+ Radius', '+ Radius, SNR', '+ All confounds']
    keys = ['simple', 'control_radius', 'control_radius_snr', 'control_all']

    # M_S
    ax = axes[0]
    r_vals = [results[k]['r_ms'] for k in keys]
    p_vals = [results[k]['p_ms'] for k in keys]
    colors = ['forestgreen' if p < 0.05 else 'indianred' for p in p_vals]

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, r_vals, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('Partial Correlation (r)', fontsize=14)
    ax.set_title('M_S ~ Depth\n(green = p < 0.05)', fontsize=14)
    ax.set_xlim(-1, 0.2)

    for i, (r, p) in enumerate(zip(r_vals, p_vals)):
        ax.text(0.05, i, f'r={r:.2f}, p={p:.3f}', va='center', ha='left', fontsize=11, fontweight='bold')

    # M_C
    ax = axes[1]
    r_vals = [results[k]['r_mc'] for k in keys]
    p_vals = [results[k]['p_mc'] for k in keys]
    colors = ['forestgreen' if p < 0.05 else 'indianred' for p in p_vals]

    bars = ax.barh(y_pos, r_vals, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('Partial Correlation (r)', fontsize=14)
    ax.set_title('M_C ~ Depth\n(green = p < 0.05)', fontsize=14)
    ax.set_xlim(-0.2, 1)

    for i, (r, p) in enumerate(zip(r_vals, p_vals)):
        ax.text(r + 0.05, i, f'r={r:.2f}, p={p:.3f}', va='center', ha='left', fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_dir / 'roi_size_control.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'roi_size_control.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'roi_size_control.png'}")

    return results

# ============================================================
# ANALYSIS 2: MIXED-EFFECTS MODEL
# ============================================================

def run_mixed_effects(df):
    """Run mixed-effects model with site as random effect"""
    print("\n" + "="*60)
    print("ANALYSIS: MIXED-EFFECTS MODELING")
    print("="*60)

    if not HAS_STATSMODELS:
        print("[SKIPPED] statsmodels not available")
        return None

    out_dir = OUTPUT_DIR / "mixed_effects"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prepare data
    model_df = df[['site_name', 'depth', 'm_s', 'm_c', 'snr_g', 'osi', 'bsf', 'radius']].dropna()

    # Standardize predictors for better convergence
    for col in ['depth', 'snr_g', 'osi', 'bsf', 'radius']:
        model_df[f'{col}_z'] = (model_df[col] - model_df[col].mean()) / model_df[col].std()

    results = {}

    # Model 1: M_S ~ Depth (random intercept for site)
    print("\n--- Model 1: M_S ~ Depth (random intercept) ---")
    try:
        model1 = smf.mixedlm("m_s ~ depth_z", model_df, groups=model_df["site_name"])
        fit1 = model1.fit(method='powell')
        print(fit1.summary().tables[1])
        results['ms_simple'] = {
            'depth_coef': fit1.params['depth_z'],
            'depth_pval': fit1.pvalues['depth_z'],
            'aic': fit1.aic
        }
        print(f"\nDepth coefficient: {fit1.params['depth_z']:.4f}")
        print(f"Depth p-value: {fit1.pvalues['depth_z']:.4e}")
    except Exception as e:
        print(f"[ERROR] {e}")
        results['ms_simple'] = None

    # Model 2: M_S ~ Depth + Radius + SNR (random intercept)
    print("\n--- Model 2: M_S ~ Depth + Radius + SNR (random intercept) ---")
    try:
        model2 = smf.mixedlm("m_s ~ depth_z + radius_z + snr_g_z", model_df, groups=model_df["site_name"])
        fit2 = model2.fit(method='powell')
        print(fit2.summary().tables[1])
        results['ms_controlled'] = {
            'depth_coef': fit2.params['depth_z'],
            'depth_pval': fit2.pvalues['depth_z'],
            'radius_coef': fit2.params['radius_z'],
            'radius_pval': fit2.pvalues['radius_z'],
            'aic': fit2.aic
        }
        print(f"\nDepth coefficient: {fit2.params['depth_z']:.4f}, p = {fit2.pvalues['depth_z']:.4e}")
        print(f"Radius coefficient: {fit2.params['radius_z']:.4f}, p = {fit2.pvalues['radius_z']:.4e}")
    except Exception as e:
        print(f"[ERROR] {e}")
        results['ms_controlled'] = None

    # Model 3: M_C ~ Depth (random intercept)
    print("\n--- Model 3: M_C ~ Depth (random intercept) ---")
    try:
        model3 = smf.mixedlm("m_c ~ depth_z", model_df, groups=model_df["site_name"])
        fit3 = model3.fit(method='powell')
        print(fit3.summary().tables[1])
        results['mc_simple'] = {
            'depth_coef': fit3.params['depth_z'],
            'depth_pval': fit3.pvalues['depth_z'],
            'aic': fit3.aic
        }
    except Exception as e:
        print(f"[ERROR] {e}")
        results['mc_simple'] = None

    # Model 4: M_C ~ Depth + Radius + SNR
    print("\n--- Model 4: M_C ~ Depth + Radius + SNR (random intercept) ---")
    try:
        model4 = smf.mixedlm("m_c ~ depth_z + radius_z + snr_g_z", model_df, groups=model_df["site_name"])
        fit4 = model4.fit(method='powell')
        print(fit4.summary().tables[1])
        results['mc_controlled'] = {
            'depth_coef': fit4.params['depth_z'],
            'depth_pval': fit4.pvalues['depth_z'],
            'radius_coef': fit4.params['radius_z'],
            'radius_pval': fit4.pvalues['radius_z'],
            'aic': fit4.aic
        }
    except Exception as e:
        print(f"[ERROR] {e}")
        results['mc_controlled'] = None

    # Save results summary
    with open(out_dir / 'mixed_effects_summary.txt', 'w') as f:
        f.write("MIXED-EFFECTS MODEL RESULTS\n")
        f.write("="*60 + "\n\n")

        f.write("Model: metric ~ fixed_effects + (1|site)\n")
        f.write("Random intercept for each recording site\n\n")

        for name, res in results.items():
            f.write(f"\n{name}:\n")
            if res:
                for k, v in res.items():
                    f.write(f"  {k}: {v}\n")
            else:
                f.write("  [FAILED]\n")

    print(f"\n[SAVED] {out_dir / 'mixed_effects_summary.txt'}")

    return results

# ============================================================
# ANALYSIS 3: OSI INTERACTION DEEP-DIVE
# ============================================================

def run_osi_interaction(df):
    """Detailed analysis of OSI × Depth interaction"""
    print("\n" + "="*60)
    print("ANALYSIS: OSI × DEPTH INTERACTION")
    print("="*60)

    out_dir = OUTPUT_DIR / "osi_interaction"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Quartile split on OSI
    df['osi_quartile'] = pd.qcut(df['osi'], q=4, labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'])

    # Compute correlations for each quartile
    results = {}
    print("\n--- M_S ~ Depth by OSI Quartile ---")

    for q in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
        subset = df[df['osi_quartile'] == q]
        site_means = subset.groupby('site_name').agg({
            'depth': 'first', 'm_s': 'mean', 'm_c': 'mean'
        }).reset_index()

        r_ms, p_ms = pearsonr(site_means['depth'], site_means['m_s'])
        r_mc, p_mc = pearsonr(site_means['depth'], site_means['m_c'])

        results[q] = {
            'n_rois': len(subset),
            'r_ms': r_ms, 'p_ms': p_ms,
            'r_mc': r_mc, 'p_mc': p_mc,
            'osi_range': (subset['osi'].min(), subset['osi'].max())
        }

        sig_ms = "***" if p_ms < 0.001 else "**" if p_ms < 0.01 else "*" if p_ms < 0.05 else ""
        sig_mc = "***" if p_mc < 0.001 else "**" if p_mc < 0.01 else "*" if p_mc < 0.05 else ""
        print(f"{q}: M_S r={r_ms:+.3f}{sig_ms}, M_C r={r_mc:+.3f}{sig_mc} (n={len(subset)} ROIs)")

    # Create figure
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))

    colors = {'Q1 (Low)': 'blue', 'Q2': 'green', 'Q3': 'orange', 'Q4 (High)': 'red'}

    for idx, q in enumerate(['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']):
        subset = df[df['osi_quartile'] == q]
        site_means = subset.groupby('site_name').agg({
            'depth': 'first', 'm_s': 'mean', 'm_c': 'mean'
        }).reset_index()

        # M_S
        ax = axes[0, idx]
        ax.scatter(site_means['depth'], site_means['m_s'], c=colors[q], s=80,
                   edgecolors='black', alpha=0.7)
        slope, intercept, _, _, _ = linregress(site_means['depth'], site_means['m_s'])
        x_fit = np.linspace(site_means['depth'].min(), site_means['depth'].max(), 100)
        ax.plot(x_fit, slope * x_fit + intercept, '--', color=colors[q], linewidth=2)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Depth (μm)', fontsize=11)
        ax.set_ylabel('M_S', fontsize=11)
        r, p = results[q]['r_ms'], results[q]['p_ms']
        ax.set_title(f'{q}\nOSI: {results[q]["osi_range"][0]:.2f}-{results[q]["osi_range"][1]:.2f}\nr={r:.3f}, p={p:.3f}', fontsize=10)

        # M_C
        ax = axes[1, idx]
        ax.scatter(site_means['depth'], site_means['m_c'], c=colors[q], s=80,
                   edgecolors='black', alpha=0.7)
        slope, intercept, _, _, _ = linregress(site_means['depth'], site_means['m_c'])
        ax.plot(x_fit, slope * x_fit + intercept, '--', color=colors[q], linewidth=2)
        ax.set_xlabel('Depth (μm)', fontsize=11)
        ax.set_ylabel('M_C', fontsize=11)
        r, p = results[q]['r_mc'], results[q]['p_mc']
        ax.set_title(f'r={r:.3f}, p={p:.3f}', fontsize=10)

    axes[0, 0].set_ylabel('M_S (a.u.)', fontsize=12, fontweight='bold')
    axes[1, 0].set_ylabel('M_C (r-value)', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_dir / 'osi_quartile_interaction.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'osi_quartile_interaction.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'osi_quartile_interaction.png'}")

    # Test for interaction: is the M_C depth effect stronger in high OSI?
    print("\n--- Testing OSI × Depth Interaction ---")
    q1_site = df[df['osi_quartile'] == 'Q1 (Low)'].groupby('site_name').agg({'depth': 'first', 'm_c': 'mean'}).reset_index()
    q4_site = df[df['osi_quartile'] == 'Q4 (High)'].groupby('site_name').agg({'depth': 'first', 'm_c': 'mean'}).reset_index()

    r_q1, _ = pearsonr(q1_site['depth'], q1_site['m_c'])
    r_q4, _ = pearsonr(q4_site['depth'], q4_site['m_c'])

    # Fisher z-test for difference in correlations
    n1, n2 = len(q1_site), len(q4_site)
    z1 = np.arctanh(r_q1)
    z4 = np.arctanh(r_q4)
    se_diff = np.sqrt(1/(n1-3) + 1/(n2-3))
    z_diff = (z4 - z1) / se_diff
    from scipy.special import erf
    p_diff = 2 * (1 - 0.5 * (1 + erf(abs(z_diff) / np.sqrt(2))))

    print(f"M_C ~ Depth correlation in Q1 (low OSI): r = {r_q1:.3f}")
    print(f"M_C ~ Depth correlation in Q4 (high OSI): r = {r_q4:.3f}")
    print(f"Difference test: z = {z_diff:.3f}, p = {p_diff:.4f}")

    if p_diff < 0.05:
        print(">>> SIGNIFICANT: M_C depth effect is stronger in high-OSI cells!")
    else:
        print(">>> Not significant difference between OSI groups")

    results['interaction_test'] = {'z': z_diff, 'p': p_diff}

    return results

# ============================================================
# ANALYSIS 4: PUBLICATION-READY SUMMARY FIGURE
# ============================================================

def create_summary_figure(df):
    """Create publication-ready summary figure"""
    print("\n" + "="*60)
    print("CREATING PUBLICATION SUMMARY FIGURE")
    print("="*60)

    out_dir = OUTPUT_DIR / "summary_figure"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': ['mean', 'sem'],
        'm_c': ['mean', 'sem']
    }).reset_index()
    site_data.columns = ['site_name', 'depth', 'm_s_mean', 'm_s_sem', 'm_c_mean', 'm_c_sem']
    site_data = site_data.sort_values('depth')

    # Statistics
    r_ms, p_ms = pearsonr(site_data['depth'], site_data['m_s_mean'])
    r_mc, p_mc = pearsonr(site_data['depth'], site_data['m_c_mean'])

    slope_ms, int_ms, _, _, _ = linregress(site_data['depth'], site_data['m_s_mean'])
    slope_mc, int_mc, _, _, _ = linregress(site_data['depth'], site_data['m_c_mean'])

    # Create figure
    fig = plt.figure(figsize=(14, 12))

    # Panel A: M_S vs Depth
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.errorbar(site_data['depth'], site_data['m_s_mean'], yerr=site_data['m_s_sem'],
                 fmt='o', markersize=10, color='#2E86AB', markeredgecolor='black',
                 markeredgewidth=1.5, capsize=4, capthick=1.5, elinewidth=1.5,
                 label='Site mean ± SEM')

    x_fit = np.linspace(site_data['depth'].min(), site_data['depth'].max(), 100)
    ax1.plot(x_fit, slope_ms * x_fit + int_ms, '--', color='black', linewidth=2.5, alpha=0.7)
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    ax1.set_xlabel('Cortical Depth (μm)', fontsize=14, fontweight='bold')
    ax1.set_ylabel(r'Metric S$_\Delta$ (a.u.)', fontsize=14, fontweight='bold')
    ax1.set_title('A', fontsize=18, fontweight='bold', loc='left')

    # Stats box
    stats_txt = f'r = {r_ms:.3f}\np = {p_ms:.1e}\nn = {len(site_data)} sites'
    ax1.text(0.05, 0.05, stats_txt, transform=ax1.transAxes, fontsize=12,
             verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white',
             edgecolor='black', alpha=0.9))

    ax1.tick_params(axis='both', labelsize=12)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.5)

    # Panel B: M_C vs Depth
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.errorbar(site_data['depth'], site_data['m_c_mean'], yerr=site_data['m_c_sem'],
                 fmt='s', markersize=10, color='#E94F37', markeredgecolor='black',
                 markeredgewidth=1.5, capsize=4, capthick=1.5, elinewidth=1.5,
                 label='Site mean ± SEM')

    ax2.plot(x_fit, slope_mc * x_fit + int_mc, '--', color='black', linewidth=2.5, alpha=0.7)

    ax2.set_xlabel('Cortical Depth (μm)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Metric C (r-value)', fontsize=14, fontweight='bold')
    ax2.set_title('B', fontsize=18, fontweight='bold', loc='left')

    stats_txt = f'r = {r_mc:.3f}\np = {p_mc:.1e}\nn = {len(site_data)} sites'
    ax2.text(0.05, 0.95, stats_txt, transform=ax2.transAxes, fontsize=12,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white',
             edgecolor='black', alpha=0.9))

    ax2.tick_params(axis='both', labelsize=12)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.5)

    # Panel C: M_S vs M_C colored by depth
    ax3 = fig.add_subplot(2, 2, 3)
    scatter = ax3.scatter(site_data['m_s_mean'], site_data['m_c_mean'],
                          c=site_data['depth'], cmap='viridis', s=150,
                          edgecolors='black', linewidths=1.5)

    r_ms_mc, p_ms_mc = pearsonr(site_data['m_s_mean'], site_data['m_c_mean'])
    slope, intercept, _, _, _ = linregress(site_data['m_s_mean'], site_data['m_c_mean'])
    x_fit = np.linspace(site_data['m_s_mean'].min(), site_data['m_s_mean'].max(), 100)
    ax3.plot(x_fit, slope * x_fit + intercept, '--', color='black', linewidth=2, alpha=0.7)

    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax3.axvline(x=0, color='gray', linestyle='--', alpha=0.3)

    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Depth (μm)', fontsize=12)

    ax3.set_xlabel(r'Metric S$_\Delta$ (a.u.)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Metric C (r-value)', fontsize=14, fontweight='bold')
    ax3.set_title('C', fontsize=18, fontweight='bold', loc='left')

    stats_txt = f'r = {r_ms_mc:.3f}\np = {p_ms_mc:.3f}'
    ax3.text(0.95, 0.05, stats_txt, transform=ax3.transAxes, fontsize=12,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.9))

    ax3.tick_params(axis='both', labelsize=12)
    for spine in ax3.spines.values():
        spine.set_linewidth(1.5)

    # Panel D: Schematic / depth bins
    ax4 = fig.add_subplot(2, 2, 4)

    # Bin by depth terciles
    depth_bins = pd.cut(site_data['depth'], bins=3, labels=['Superficial\n(140-266μm)',
                                                             'Middle\n(266-392μm)',
                                                             'Deep\n(392-518μm)'])
    binned = site_data.copy()
    binned['depth_bin'] = depth_bins

    bin_means = binned.groupby('depth_bin').agg({
        'm_s_mean': ['mean', 'std'],
        'm_c_mean': ['mean', 'std']
    })
    bin_means.columns = ['m_s', 'm_s_std', 'm_c', 'm_c_std']
    bin_means = bin_means.reset_index()

    x = np.arange(3)
    width = 0.35

    bars1 = ax4.bar(x - width/2, bin_means['m_s'], width, yerr=bin_means['m_s_std'],
                    label='M_S', color='#2E86AB', edgecolor='black', linewidth=1.5,
                    capsize=5, error_kw={'linewidth': 1.5})
    bars2 = ax4.bar(x + width/2, bin_means['m_c'] * 50, width, yerr=bin_means['m_c_std'] * 50,
                    label='M_C (×50)', color='#E94F37', edgecolor='black', linewidth=1.5,
                    capsize=5, error_kw={'linewidth': 1.5})

    ax4.axhline(y=0, color='black', linewidth=1)
    ax4.set_xticks(x)
    ax4.set_xticklabels(bin_means['depth_bin'], fontsize=11)
    ax4.set_ylabel('Mean Value', fontsize=14, fontweight='bold')
    ax4.set_title('D', fontsize=18, fontweight='bold', loc='left')
    ax4.legend(loc='upper right', fontsize=11)

    ax4.tick_params(axis='both', labelsize=12)
    for spine in ax4.spines.values():
        spine.set_linewidth(1.5)

    plt.tight_layout()
    plt.savefig(out_dir / 'publication_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig(out_dir / 'publication_summary.svg', format='svg', bbox_inches='tight')
    plt.savefig(out_dir / 'publication_summary.pdf', format='pdf', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'publication_summary.png'}")
    print(f"[SAVED] {out_dir / 'publication_summary.pdf'}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("="*60)
    print("XORI ADDITIONAL ANALYSES")
    print("="*60)

    # Load data
    print("\nLoading data...")
    df, depths = load_all_data()
    print(f"Loaded {len(df)} ROIs from {len(depths)} sites")

    # Run analyses
    roi_size_results = run_roi_size_control(df)
    mixed_results = run_mixed_effects(df)
    osi_results = run_osi_interaction(df)
    create_summary_figure(df)

    print("\n" + "="*60)
    print("ADDITIONAL ANALYSES COMPLETE")
    print("="*60)

if __name__ == '__main__':
    main()
