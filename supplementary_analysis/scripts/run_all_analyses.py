#!/usr/bin/env python3
"""
Supplementary Analysis for XORI Project
Runs all additional analyses to strengthen the depth-metric findings.

Analyses included:
1. Subpopulation analysis (high/low OSI, high/low SF)
2. Partial correlations controlling for confounds
3. Within-site variance analysis
4. M_S vs M_C relationship
5. Morphology confound checks
6. Bootstrap confidence intervals
7. ROI-level scatter plots
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set up paths
BASE_DIR = Path("/Users/muditagar/XORI_Analysis/XORI")
RAW_DATA = BASE_DIR / "raw_data" / "bm_data"
METRIC_DATA = BASE_DIR / "metric_data" / "all_roi"
STAT_DIR = BASE_DIR / "stat"
OUTPUT_DIR = BASE_DIR / "supplementary_analysis" / "outputs"

# Results storage
RESULTS = {}

# ============================================================
# DATA LOADING
# ============================================================

def load_site_depths():
    """Load site depth mapping"""
    depths = {}
    with open(RAW_DATA / "site_depth.txt", 'r') as f:
        lines = f.readlines()
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            depths[parts[0]] = float(parts[1])
    return depths

def load_osi_data():
    """Load OSI data for all ROIs"""
    data = []
    with open(RAW_DATA / "roi_osi.txt", 'r') as f:
        lines = f.readlines()
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 8:
            try:
                data.append({
                    'site': int(parts[0]),
                    'roi': int(parts[1]),
                    'bsf': float(parts[2]),
                    'bdr': float(parts[3]),
                    'dsi': float(parts[4]),
                    'osi': float(parts[5]),
                    'dcv': float(parts[6]),
                    'ocv': float(parts[7])
                })
            except ValueError:
                continue
    return pd.DataFrame(data)

def load_roi_stat():
    """Load ROI morphology statistics"""
    data = []
    with open(RAW_DATA / "roi_stat.txt", 'r') as f:
        lines = f.readlines()
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 14:
            try:
                data.append({
                    'site': int(parts[0]),
                    'roi': int(parts[1]),
                    'npix': int(parts[2]),
                    'medx': float(parts[3]),
                    'medy': float(parts[4]),
                    'neuropil_pix': int(parts[5]),
                    'radius': float(parts[6]),
                    'aspect_ratio': float(parts[7]),
                    'skew': float(parts[8]),
                    'std': float(parts[9]),
                    'mrs': float(parts[10]),
                    'mrs0': float(parts[11]),
                    'compact': float(parts[12]),
                    'solidity': float(parts[13])
                })
            except ValueError:
                continue
    return pd.DataFrame(data)

def load_metric_data(site_depths):
    """Load all metric data and merge with site depths"""
    all_data = []

    for metric_file in METRIC_DATA.glob("metrics_site*.txt"):
        site_name = metric_file.stem.replace("metrics_", "")
        site_num = int(site_name.replace("site", ""))

        if site_name not in site_depths:
            continue

        depth = site_depths[site_name]

        with open(metric_file, 'r') as f:
            lines = f.readlines()

        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) >= 8:
                try:
                    all_data.append({
                        'site': site_num,
                        'site_name': site_name,
                        'roi': int(parts[0]),
                        'depth': depth,
                        'm_s': float(parts[1]),
                        'm_c': float(parts[2]),
                        'snr_g': float(parts[3]),
                        'snr_p': float(parts[4]),
                        'm_s_norm': float(parts[5]),
                        'm_s_ratio': float(parts[6]),
                        'm_x': float(parts[7])
                    })
                except ValueError:
                    continue

    return pd.DataFrame(all_data)

def merge_all_data(metric_df, osi_df, stat_df):
    """Merge all dataframes"""
    # Merge metric with OSI
    merged = metric_df.merge(osi_df, on=['site', 'roi'], how='left')
    # Merge with morphology stats
    merged = merged.merge(stat_df, on=['site', 'roi'], how='left')
    return merged

# ============================================================
# ANALYSIS 1: SUBPOPULATION ANALYSIS
# ============================================================

def run_subpopulation_analysis(df):
    """Split by OSI and SF, analyze depth effects separately"""
    print("\n" + "="*60)
    print("ANALYSIS 1: SUBPOPULATION ANALYSIS")
    print("="*60)

    results = {}
    out_dir = OUTPUT_DIR / "subpopulation"
    out_dir.mkdir(parents=True, exist_ok=True)

    # OSI split (median split)
    osi_median = df['osi'].median()
    high_osi = df[df['osi'] >= osi_median]
    low_osi = df[df['osi'] < osi_median]

    print(f"\nOSI median: {osi_median:.3f}")
    print(f"High OSI ROIs: {len(high_osi)}, Low OSI ROIs: {len(low_osi)}")

    # SF split (median split)
    sf_median = df['bsf'].median()
    high_sf = df[df['bsf'] >= sf_median]
    low_sf = df[df['bsf'] < sf_median]

    print(f"\nSF median: {sf_median:.3f} cpd")
    print(f"High SF ROIs: {len(high_sf)}, Low SF ROIs: {len(low_sf)}")

    # Function to compute site-level stats and correlation
    def compute_depth_correlation(subset, label):
        site_means = subset.groupby('site_name').agg({
            'depth': 'first',
            'm_s': 'mean',
            'm_c': 'mean'
        }).reset_index()

        if len(site_means) < 5:
            return None

        r_ms, p_ms = pearsonr(site_means['depth'], site_means['m_s'])
        r_mc, p_mc = pearsonr(site_means['depth'], site_means['m_c'])

        return {
            'label': label,
            'n_rois': len(subset),
            'n_sites': len(site_means),
            'r_ms': r_ms,
            'p_ms': p_ms,
            'r_mc': r_mc,
            'p_mc': p_mc,
            'site_means': site_means
        }

    # Compute for each subpopulation
    results['high_osi'] = compute_depth_correlation(high_osi, 'High OSI')
    results['low_osi'] = compute_depth_correlation(low_osi, 'Low OSI')
    results['high_sf'] = compute_depth_correlation(high_sf, 'High SF')
    results['low_sf'] = compute_depth_correlation(low_sf, 'Low SF')
    results['all'] = compute_depth_correlation(df, 'All ROIs')

    # Print results
    print("\n--- M_S vs Depth by Subpopulation ---")
    for key, res in results.items():
        if res:
            sig = "***" if res['p_ms'] < 0.001 else "**" if res['p_ms'] < 0.01 else "*" if res['p_ms'] < 0.05 else ""
            print(f"{res['label']:12s}: r = {res['r_ms']:+.3f}, p = {res['p_ms']:.2e} {sig} (n={res['n_rois']} ROIs, {res['n_sites']} sites)")

    print("\n--- M_C vs Depth by Subpopulation ---")
    for key, res in results.items():
        if res:
            sig = "***" if res['p_mc'] < 0.001 else "**" if res['p_mc'] < 0.01 else "*" if res['p_mc'] < 0.05 else ""
            print(f"{res['label']:12s}: r = {res['r_mc']:+.3f}, p = {res['p_mc']:.2e} {sig} (n={res['n_rois']} ROIs, {res['n_sites']} sites)")

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Row 1: M_S vs Depth
    for idx, (key, color, marker) in enumerate([('high_osi', 'red', 'o'), ('low_osi', 'blue', 's')]):
        if results[key]:
            sm = results[key]['site_means']
            ax = axes[0, idx]
            ax.scatter(sm['depth'], sm['m_s'], c=color, s=80, alpha=0.7, edgecolors='black', marker=marker)

            # Regression line
            slope, intercept, _, _, _ = linregress(sm['depth'], sm['m_s'])
            x_fit = np.linspace(sm['depth'].min(), sm['depth'].max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, '--', color=color, alpha=0.5, linewidth=2)

            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel('Depth (μm)', fontsize=12)
            ax.set_ylabel('M_S (a.u.)', fontsize=12)
            ax.set_title(f"{results[key]['label']}\nr = {results[key]['r_ms']:.3f}, p = {results[key]['p_ms']:.2e}", fontsize=11)

    # Comparison panel
    ax = axes[0, 2]
    for key, color, marker, label in [('high_osi', 'red', 'o', 'High OSI'), ('low_osi', 'blue', 's', 'Low OSI')]:
        if results[key]:
            sm = results[key]['site_means']
            ax.scatter(sm['depth'], sm['m_s'], c=color, s=60, alpha=0.5, marker=marker, label=label)
            slope, intercept, _, _, _ = linregress(sm['depth'], sm['m_s'])
            x_fit = np.linspace(sm['depth'].min(), sm['depth'].max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, '--', color=color, alpha=0.7, linewidth=2)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('M_S (a.u.)', fontsize=12)
    ax.set_title('OSI Comparison', fontsize=11)
    ax.legend()

    # Row 2: SF splits
    for idx, (key, color, marker) in enumerate([('high_sf', 'green', '^'), ('low_sf', 'purple', 'v')]):
        if results[key]:
            sm = results[key]['site_means']
            ax = axes[1, idx]
            ax.scatter(sm['depth'], sm['m_s'], c=color, s=80, alpha=0.7, edgecolors='black', marker=marker)

            slope, intercept, _, _, _ = linregress(sm['depth'], sm['m_s'])
            x_fit = np.linspace(sm['depth'].min(), sm['depth'].max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, '--', color=color, alpha=0.5, linewidth=2)

            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel('Depth (μm)', fontsize=12)
            ax.set_ylabel('M_S (a.u.)', fontsize=12)
            ax.set_title(f"{results[key]['label']}\nr = {results[key]['r_ms']:.3f}, p = {results[key]['p_ms']:.2e}", fontsize=11)

    # SF Comparison panel
    ax = axes[1, 2]
    for key, color, marker, label in [('high_sf', 'green', '^', 'High SF'), ('low_sf', 'purple', 'v', 'Low SF')]:
        if results[key]:
            sm = results[key]['site_means']
            ax.scatter(sm['depth'], sm['m_s'], c=color, s=60, alpha=0.5, marker=marker, label=label)
            slope, intercept, _, _, _ = linregress(sm['depth'], sm['m_s'])
            x_fit = np.linspace(sm['depth'].min(), sm['depth'].max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, '--', color=color, alpha=0.7, linewidth=2)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('M_S (a.u.)', fontsize=12)
    ax.set_title('SF Comparison', fontsize=11)
    ax.legend()

    plt.tight_layout()
    plt.savefig(out_dir / 'subpopulation_ms_vs_depth.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'subpopulation_ms_vs_depth.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'subpopulation_ms_vs_depth.png'}")

    RESULTS['subpopulation'] = results
    return results

# ============================================================
# ANALYSIS 2: PARTIAL CORRELATIONS
# ============================================================

def partial_correlation(x, y, covariates):
    """Compute partial correlation between x and y controlling for covariates"""
    from numpy.linalg import lstsq

    # Residualize x
    X_cov = np.column_stack([np.ones(len(covariates)), covariates])
    beta_x, _, _, _ = lstsq(X_cov, x, rcond=None)
    x_resid = x - X_cov @ beta_x

    # Residualize y
    beta_y, _, _, _ = lstsq(X_cov, y, rcond=None)
    y_resid = y - X_cov @ beta_y

    # Correlation of residuals
    r, p = pearsonr(x_resid, y_resid)
    return r, p

def run_partial_correlation_analysis(df):
    """Run partial correlations controlling for various confounds"""
    print("\n" + "="*60)
    print("ANALYSIS 2: PARTIAL CORRELATION ANALYSIS")
    print("="*60)

    out_dir = OUTPUT_DIR / "partial_correlation"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Work at site level
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'snr_g': 'mean',
        'osi': 'mean',
        'bsf': 'mean',
        'radius': 'mean',
        'npix': 'mean'
    }).reset_index().dropna()

    results = {}

    # Simple correlations first
    print("\n--- Simple Correlations (no controls) ---")
    r_ms, p_ms = pearsonr(site_data['depth'], site_data['m_s'])
    r_mc, p_mc = pearsonr(site_data['depth'], site_data['m_c'])
    print(f"M_S ~ Depth: r = {r_ms:+.3f}, p = {p_ms:.2e}")
    print(f"M_C ~ Depth: r = {r_mc:+.3f}, p = {p_mc:.2e}")
    results['simple'] = {'r_ms': r_ms, 'p_ms': p_ms, 'r_mc': r_mc, 'p_mc': p_mc}

    # Partial correlations
    print("\n--- Partial Correlations ---")

    # Control for SNR
    covariates_snr = site_data['snr_g'].values.reshape(-1, 1)
    r_ms_snr, p_ms_snr = partial_correlation(site_data['depth'].values, site_data['m_s'].values, covariates_snr)
    r_mc_snr, p_mc_snr = partial_correlation(site_data['depth'].values, site_data['m_c'].values, covariates_snr)
    print(f"M_S ~ Depth | SNR: r = {r_ms_snr:+.3f}, p = {p_ms_snr:.2e}")
    print(f"M_C ~ Depth | SNR: r = {r_mc_snr:+.3f}, p = {p_mc_snr:.2e}")
    results['control_snr'] = {'r_ms': r_ms_snr, 'p_ms': p_ms_snr, 'r_mc': r_mc_snr, 'p_mc': p_mc_snr}

    # Control for OSI
    covariates_osi = site_data['osi'].values.reshape(-1, 1)
    r_ms_osi, p_ms_osi = partial_correlation(site_data['depth'].values, site_data['m_s'].values, covariates_osi)
    r_mc_osi, p_mc_osi = partial_correlation(site_data['depth'].values, site_data['m_c'].values, covariates_osi)
    print(f"M_S ~ Depth | OSI: r = {r_ms_osi:+.3f}, p = {p_ms_osi:.2e}")
    print(f"M_C ~ Depth | OSI: r = {r_mc_osi:+.3f}, p = {p_mc_osi:.2e}")
    results['control_osi'] = {'r_ms': r_ms_osi, 'p_ms': p_ms_osi, 'r_mc': r_mc_osi, 'p_mc': p_mc_osi}

    # Control for SF
    covariates_sf = site_data['bsf'].values.reshape(-1, 1)
    r_ms_sf, p_ms_sf = partial_correlation(site_data['depth'].values, site_data['m_s'].values, covariates_sf)
    r_mc_sf, p_mc_sf = partial_correlation(site_data['depth'].values, site_data['m_c'].values, covariates_sf)
    print(f"M_S ~ Depth | SF: r = {r_ms_sf:+.3f}, p = {p_ms_sf:.2e}")
    print(f"M_C ~ Depth | SF: r = {r_mc_sf:+.3f}, p = {p_mc_sf:.2e}")
    results['control_sf'] = {'r_ms': r_ms_sf, 'p_ms': p_ms_sf, 'r_mc': r_mc_sf, 'p_mc': p_mc_sf}

    # Control for ALL
    covariates_all = np.column_stack([site_data['snr_g'], site_data['osi'], site_data['bsf']])
    r_ms_all, p_ms_all = partial_correlation(site_data['depth'].values, site_data['m_s'].values, covariates_all)
    r_mc_all, p_mc_all = partial_correlation(site_data['depth'].values, site_data['m_c'].values, covariates_all)
    print(f"M_S ~ Depth | SNR+OSI+SF: r = {r_ms_all:+.3f}, p = {p_ms_all:.2e}")
    print(f"M_C ~ Depth | SNR+OSI+SF: r = {r_mc_all:+.3f}, p = {p_mc_all:.2e}")
    results['control_all'] = {'r_ms': r_ms_all, 'p_ms': p_ms_all, 'r_mc': r_mc_all, 'p_mc': p_mc_all}

    # Create forest plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    labels = ['No control', 'Control SNR', 'Control OSI', 'Control SF', 'Control ALL']
    keys = ['simple', 'control_snr', 'control_osi', 'control_sf', 'control_all']

    # M_S plot
    ax = axes[0]
    r_vals = [results[k]['r_ms'] for k in keys]
    p_vals = [results[k]['p_ms'] for k in keys]
    colors = ['green' if p < 0.05 else 'red' for p in p_vals]

    y_pos = np.arange(len(labels))
    ax.barh(y_pos, r_vals, color=colors, alpha=0.7, edgecolor='black')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Correlation (r)')
    ax.set_title('M_S ~ Depth\n(green = p < 0.05)')
    ax.set_xlim(-1, 1)

    for i, (r, p) in enumerate(zip(r_vals, p_vals)):
        ax.text(r + 0.05 if r >= 0 else r - 0.05, i, f'p={p:.3f}',
                va='center', ha='left' if r >= 0 else 'right', fontsize=9)

    # M_C plot
    ax = axes[1]
    r_vals = [results[k]['r_mc'] for k in keys]
    p_vals = [results[k]['p_mc'] for k in keys]
    colors = ['green' if p < 0.05 else 'red' for p in p_vals]

    ax.barh(y_pos, r_vals, color=colors, alpha=0.7, edgecolor='black')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Correlation (r)')
    ax.set_title('M_C ~ Depth\n(green = p < 0.05)')
    ax.set_xlim(-1, 1)

    for i, (r, p) in enumerate(zip(r_vals, p_vals)):
        ax.text(r + 0.05 if r >= 0 else r - 0.05, i, f'p={p:.3f}',
                va='center', ha='left' if r >= 0 else 'right', fontsize=9)

    plt.tight_layout()
    plt.savefig(out_dir / 'partial_correlations_forest.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'partial_correlations_forest.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'partial_correlations_forest.png'}")

    RESULTS['partial_correlation'] = results
    return results

# ============================================================
# ANALYSIS 3: WITHIN-SITE VARIANCE
# ============================================================

def run_within_site_analysis(df):
    """Analyze within-site variance and create error bar plots"""
    print("\n" + "="*60)
    print("ANALYSIS 3: WITHIN-SITE VARIANCE ANALYSIS")
    print("="*60)

    out_dir = OUTPUT_DIR / "within_site"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Compute site-level statistics
    site_stats = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': ['mean', 'std', 'sem', 'count'],
        'm_c': ['mean', 'std', 'sem', 'count']
    }).reset_index()

    site_stats.columns = ['site_name', 'depth', 'm_s_mean', 'm_s_std', 'm_s_sem', 'm_s_n',
                          'm_c_mean', 'm_c_std', 'm_c_sem', 'm_c_n']

    site_stats = site_stats.sort_values('depth')

    print(f"\nSite-level statistics:")
    print(f"Mean ROIs per site: {site_stats['m_s_n'].mean():.1f}")
    print(f"Min ROIs: {site_stats['m_s_n'].min()}, Max ROIs: {site_stats['m_s_n'].max()}")
    print(f"\nMean within-site SD for M_S: {site_stats['m_s_std'].mean():.2f}")
    print(f"Mean within-site SD for M_C: {site_stats['m_c_std'].mean():.3f}")

    # Create figure with error bars
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # M_S with error bars
    ax = axes[0]
    ax.errorbar(site_stats['depth'], site_stats['m_s_mean'],
                yerr=site_stats['m_s_sem'], fmt='o', color='darkorange',
                markersize=10, capsize=4, capthick=1.5, elinewidth=1.5,
                markeredgecolor='black', markeredgewidth=1)

    # Regression
    r, p = pearsonr(site_stats['depth'], site_stats['m_s_mean'])
    slope, intercept, _, _, _ = linregress(site_stats['depth'], site_stats['m_s_mean'])
    x_fit = np.linspace(site_stats['depth'].min(), site_stats['depth'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='black', alpha=0.5, linewidth=2)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (μm)', fontsize=14)
    ax.set_ylabel('M_S (mean ± SEM)', fontsize=14)
    ax.set_title(f'M_S vs Depth with Within-Site Error\nr = {r:.3f}, p = {p:.2e}', fontsize=12)

    # M_C with error bars
    ax = axes[1]
    ax.errorbar(site_stats['depth'], site_stats['m_c_mean'],
                yerr=site_stats['m_c_sem'], fmt='o', color='darkorange',
                markersize=10, capsize=4, capthick=1.5, elinewidth=1.5,
                markeredgecolor='black', markeredgewidth=1)

    r, p = pearsonr(site_stats['depth'], site_stats['m_c_mean'])
    slope, intercept, _, _, _ = linregress(site_stats['depth'], site_stats['m_c_mean'])
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='black', alpha=0.5, linewidth=2)

    ax.set_xlabel('Depth (μm)', fontsize=14)
    ax.set_ylabel('M_C (mean ± SEM)', fontsize=14)
    ax.set_title(f'M_C vs Depth with Within-Site Error\nr = {r:.3f}, p = {p:.2e}', fontsize=12)

    plt.tight_layout()
    plt.savefig(out_dir / 'depth_vs_metric_with_errorbars.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'depth_vs_metric_with_errorbars.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'depth_vs_metric_with_errorbars.png'}")

    RESULTS['within_site'] = site_stats.to_dict()
    return site_stats

# ============================================================
# ANALYSIS 4: M_S vs M_C RELATIONSHIP
# ============================================================

def run_ms_mc_analysis(df):
    """Analyze relationship between M_S and M_C"""
    print("\n" + "="*60)
    print("ANALYSIS 4: M_S vs M_C RELATIONSHIP")
    print("="*60)

    out_dir = OUTPUT_DIR / "ms_mc_relationship"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level analysis
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean'
    }).reset_index()

    # Correlation between M_S and M_C
    r_site, p_site = pearsonr(site_data['m_s'], site_data['m_c'])
    print(f"\nSite-level M_S vs M_C: r = {r_site:+.3f}, p = {p_site:.2e}")

    # ROI-level correlation
    r_roi, p_roi = pearsonr(df['m_s'], df['m_c'])
    print(f"ROI-level M_S vs M_C: r = {r_roi:+.3f}, p = {p_roi:.2e}")

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Site-level M_S vs M_C colored by depth
    ax = axes[0]
    scatter = ax.scatter(site_data['m_s'], site_data['m_c'],
                         c=site_data['depth'], cmap='viridis',
                         s=120, edgecolors='black', linewidths=1.5)
    plt.colorbar(scatter, ax=ax, label='Depth (μm)')

    # Regression line
    slope, intercept, _, _, _ = linregress(site_data['m_s'], site_data['m_c'])
    x_fit = np.linspace(site_data['m_s'].min(), site_data['m_s'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2, alpha=0.7)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('M_S (a.u.)', fontsize=12)
    ax.set_ylabel('M_C (r-value)', fontsize=12)
    ax.set_title(f'Site-level: M_S vs M_C\nr = {r_site:.3f}, p = {p_site:.2e}', fontsize=11)

    # Panel B: Depth bins - how does M_S vs M_C relationship change?
    ax = axes[1]
    depth_bins = pd.qcut(df['depth'], q=3, labels=['Superficial', 'Middle', 'Deep'])
    colors = {'Superficial': 'blue', 'Middle': 'green', 'Deep': 'red'}

    for label in ['Superficial', 'Middle', 'Deep']:
        subset = df[depth_bins == label]
        r, p = pearsonr(subset['m_s'], subset['m_c'])
        ax.scatter(subset['m_s'], subset['m_c'], c=colors[label], alpha=0.1, s=10, label=f'{label} (r={r:.2f})')

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('M_S (a.u.)', fontsize=12)
    ax.set_ylabel('M_C (r-value)', fontsize=12)
    ax.set_title('ROI-level: M_S vs M_C by Depth Tercile', fontsize=11)
    ax.legend(fontsize=9)

    # Panel C: Depth effect visualization
    ax = axes[2]

    # Bin depths and compute means
    depth_bins_cont = pd.cut(df['depth'], bins=6)
    binned = df.groupby(depth_bins_cont).agg({
        'depth': 'mean',
        'm_s': 'mean',
        'm_c': 'mean'
    }).dropna()

    ax.scatter(binned['depth'], binned['m_s'], c='blue', s=150, marker='o', label='M_S', edgecolors='black')
    ax.scatter(binned['depth'], binned['m_c'] * 50, c='red', s=150, marker='s', label='M_C (×50)', edgecolors='black')

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('Mean Value', fontsize=12)
    ax.set_title('M_S and M_C Trends Across Depth', fontsize=11)
    ax.legend()

    plt.tight_layout()
    plt.savefig(out_dir / 'ms_mc_relationship.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'ms_mc_relationship.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'ms_mc_relationship.png'}")

    RESULTS['ms_mc'] = {'r_site': r_site, 'p_site': p_site, 'r_roi': r_roi, 'p_roi': p_roi}
    return RESULTS['ms_mc']

# ============================================================
# ANALYSIS 5: MORPHOLOGY CONFOUNDS
# ============================================================

def run_morphology_analysis(df):
    """Check for morphology-related confounds"""
    print("\n" + "="*60)
    print("ANALYSIS 5: MORPHOLOGY CONFOUND ANALYSIS")
    print("="*60)

    out_dir = OUTPUT_DIR / "morphology_controls"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level morphology
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'radius': 'mean',
        'npix': 'mean',
        'aspect_ratio': 'mean',
        'snr_g': 'mean'
    }).reset_index().dropna()

    results = {}

    # Check correlations
    print("\n--- Morphology vs Depth ---")
    for col in ['radius', 'npix', 'aspect_ratio', 'snr_g']:
        r, p = pearsonr(site_data['depth'], site_data[col])
        sig = "*" if p < 0.05 else ""
        print(f"{col:15s} ~ Depth: r = {r:+.3f}, p = {p:.3f} {sig}")
        results[f'depth_vs_{col}'] = {'r': r, 'p': p}

    print("\n--- Morphology vs M_S ---")
    for col in ['radius', 'npix', 'aspect_ratio', 'snr_g']:
        r, p = pearsonr(site_data['m_s'], site_data[col])
        sig = "*" if p < 0.05 else ""
        print(f"{col:15s} ~ M_S: r = {r:+.3f}, p = {p:.3f} {sig}")
        results[f'ms_vs_{col}'] = {'r': r, 'p': p}

    print("\n--- Morphology vs M_C ---")
    for col in ['radius', 'npix', 'aspect_ratio', 'snr_g']:
        r, p = pearsonr(site_data['m_c'], site_data[col])
        sig = "*" if p < 0.05 else ""
        print(f"{col:15s} ~ M_C: r = {r:+.3f}, p = {p:.3f} {sig}")
        results[f'mc_vs_{col}'] = {'r': r, 'p': p}

    # Create correlation matrix figure
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    pairs = [
        ('depth', 'radius', 'Depth vs ROI Radius'),
        ('depth', 'npix', 'Depth vs ROI Size (pixels)'),
        ('depth', 'snr_g', 'Depth vs SNR'),
        ('m_s', 'radius', 'M_S vs ROI Radius'),
        ('m_s', 'snr_g', 'M_S vs SNR'),
        ('m_c', 'snr_g', 'M_C vs SNR'),
    ]

    for ax, (x_col, y_col, title) in zip(axes.flat, pairs):
        ax.scatter(site_data[x_col], site_data[y_col],
                   c='steelblue', s=80, alpha=0.7, edgecolors='black')
        r, p = pearsonr(site_data[x_col], site_data[y_col])

        # Regression line
        slope, intercept, _, _, _ = linregress(site_data[x_col], site_data[y_col])
        x_fit = np.linspace(site_data[x_col].min(), site_data[x_col].max(), 100)
        ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', alpha=0.5, linewidth=2)

        ax.set_xlabel(x_col, fontsize=11)
        ax.set_ylabel(y_col, fontsize=11)
        color = 'red' if p < 0.05 else 'black'
        ax.set_title(f'{title}\nr = {r:.3f}, p = {p:.3f}', fontsize=10, color=color)

    plt.tight_layout()
    plt.savefig(out_dir / 'morphology_confounds.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'morphology_confounds.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'morphology_confounds.png'}")

    RESULTS['morphology'] = results
    return results

# ============================================================
# ANALYSIS 6: BOOTSTRAP STATISTICS
# ============================================================

def run_bootstrap_analysis(df, n_bootstrap=10000):
    """Bootstrap confidence intervals for depth correlations"""
    print("\n" + "="*60)
    print("ANALYSIS 6: BOOTSTRAP STATISTICS")
    print("="*60)

    out_dir = OUTPUT_DIR / "bootstrap_stats"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean'
    }).reset_index()

    n_sites = len(site_data)

    # Bootstrap
    print(f"\nRunning {n_bootstrap} bootstrap iterations...")

    boot_r_ms = []
    boot_r_mc = []
    boot_slope_ms = []
    boot_slope_mc = []

    np.random.seed(42)
    for i in range(n_bootstrap):
        # Resample sites with replacement
        idx = np.random.choice(n_sites, size=n_sites, replace=True)
        boot_sample = site_data.iloc[idx]

        r_ms, _ = pearsonr(boot_sample['depth'], boot_sample['m_s'])
        r_mc, _ = pearsonr(boot_sample['depth'], boot_sample['m_c'])

        slope_ms, _, _, _, _ = linregress(boot_sample['depth'], boot_sample['m_s'])
        slope_mc, _, _, _, _ = linregress(boot_sample['depth'], boot_sample['m_c'])

        boot_r_ms.append(r_ms)
        boot_r_mc.append(r_mc)
        boot_slope_ms.append(slope_ms)
        boot_slope_mc.append(slope_mc)

    boot_r_ms = np.array(boot_r_ms)
    boot_r_mc = np.array(boot_r_mc)
    boot_slope_ms = np.array(boot_slope_ms)
    boot_slope_mc = np.array(boot_slope_mc)

    # Compute CIs
    ci_95 = lambda x: (np.percentile(x, 2.5), np.percentile(x, 97.5))

    r_ms_ci = ci_95(boot_r_ms)
    r_mc_ci = ci_95(boot_r_mc)
    slope_ms_ci = ci_95(boot_slope_ms)
    slope_mc_ci = ci_95(boot_slope_mc)

    # Original values
    r_ms_orig, _ = pearsonr(site_data['depth'], site_data['m_s'])
    r_mc_orig, _ = pearsonr(site_data['depth'], site_data['m_c'])
    slope_ms_orig, _, _, _, _ = linregress(site_data['depth'], site_data['m_s'])
    slope_mc_orig, _, _, _, _ = linregress(site_data['depth'], site_data['m_c'])

    print("\n--- Bootstrap Results (95% CI) ---")
    print(f"M_S ~ Depth:")
    print(f"  r = {r_ms_orig:.3f} [{r_ms_ci[0]:.3f}, {r_ms_ci[1]:.3f}]")
    print(f"  slope = {slope_ms_orig:.4f} [{slope_ms_ci[0]:.4f}, {slope_ms_ci[1]:.4f}]")
    print(f"M_C ~ Depth:")
    print(f"  r = {r_mc_orig:.3f} [{r_mc_ci[0]:.3f}, {r_mc_ci[1]:.3f}]")
    print(f"  slope = {slope_mc_orig:.5f} [{slope_mc_ci[0]:.5f}, {slope_mc_ci[1]:.5f}]")

    # Check if CIs exclude zero
    print("\n--- Does 95% CI exclude zero? ---")
    print(f"M_S slope: {'YES' if (slope_ms_ci[0] > 0 or slope_ms_ci[1] < 0) else 'NO'}")
    print(f"M_C slope: {'YES' if (slope_mc_ci[0] > 0 or slope_mc_ci[1] < 0) else 'NO'}")

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # M_S correlation distribution
    ax = axes[0, 0]
    ax.hist(boot_r_ms, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(x=r_ms_orig, color='red', linestyle='-', linewidth=2, label=f'Observed: {r_ms_orig:.3f}')
    ax.axvline(x=r_ms_ci[0], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=r_ms_ci[1], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Correlation (r)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'M_S ~ Depth: Bootstrap Distribution\n95% CI: [{r_ms_ci[0]:.3f}, {r_ms_ci[1]:.3f}]', fontsize=11)
    ax.legend()

    # M_S slope distribution
    ax = axes[0, 1]
    ax.hist(boot_slope_ms, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(x=slope_ms_orig, color='red', linestyle='-', linewidth=2, label=f'Observed: {slope_ms_orig:.4f}')
    ax.axvline(x=slope_ms_ci[0], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=slope_ms_ci[1], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Slope', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'M_S Slope: Bootstrap Distribution\n95% CI: [{slope_ms_ci[0]:.4f}, {slope_ms_ci[1]:.4f}]', fontsize=11)
    ax.legend()

    # M_C correlation distribution
    ax = axes[1, 0]
    ax.hist(boot_r_mc, bins=50, color='darkorange', alpha=0.7, edgecolor='black')
    ax.axvline(x=r_mc_orig, color='red', linestyle='-', linewidth=2, label=f'Observed: {r_mc_orig:.3f}')
    ax.axvline(x=r_mc_ci[0], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=r_mc_ci[1], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Correlation (r)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'M_C ~ Depth: Bootstrap Distribution\n95% CI: [{r_mc_ci[0]:.3f}, {r_mc_ci[1]:.3f}]', fontsize=11)
    ax.legend()

    # M_C slope distribution
    ax = axes[1, 1]
    ax.hist(boot_slope_mc, bins=50, color='darkorange', alpha=0.7, edgecolor='black')
    ax.axvline(x=slope_mc_orig, color='red', linestyle='-', linewidth=2, label=f'Observed: {slope_mc_orig:.5f}')
    ax.axvline(x=slope_mc_ci[0], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=slope_mc_ci[1], color='red', linestyle='--', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Slope', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'M_C Slope: Bootstrap Distribution\n95% CI: [{slope_mc_ci[0]:.5f}, {slope_mc_ci[1]:.5f}]', fontsize=11)
    ax.legend()

    plt.tight_layout()
    plt.savefig(out_dir / 'bootstrap_distributions.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'bootstrap_distributions.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'bootstrap_distributions.png'}")

    results = {
        'r_ms': {'value': r_ms_orig, 'ci_low': r_ms_ci[0], 'ci_high': r_ms_ci[1]},
        'r_mc': {'value': r_mc_orig, 'ci_low': r_mc_ci[0], 'ci_high': r_mc_ci[1]},
        'slope_ms': {'value': slope_ms_orig, 'ci_low': slope_ms_ci[0], 'ci_high': slope_ms_ci[1]},
        'slope_mc': {'value': slope_mc_orig, 'ci_low': slope_mc_ci[0], 'ci_high': slope_mc_ci[1]}
    }

    RESULTS['bootstrap'] = results
    return results

# ============================================================
# ANALYSIS 7: ROI-LEVEL SCATTER
# ============================================================

def run_roi_scatter_analysis(df):
    """Create ROI-level scatter plots"""
    print("\n" + "="*60)
    print("ANALYSIS 7: ROI-LEVEL SCATTER PLOTS")
    print("="*60)

    out_dir = OUTPUT_DIR / "roi_scatter"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site means for overlay
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean'
    }).reset_index()

    print(f"\nTotal ROIs: {len(df)}")
    print(f"Total sites: {len(site_data)}")

    # ROI-level correlations
    r_ms_roi, p_ms_roi = pearsonr(df['depth'], df['m_s'])
    r_mc_roi, p_mc_roi = pearsonr(df['depth'], df['m_c'])

    print(f"\nROI-level correlations:")
    print(f"M_S ~ Depth: r = {r_ms_roi:.3f}, p = {p_ms_roi:.2e}")
    print(f"M_C ~ Depth: r = {r_mc_roi:.3f}, p = {p_mc_roi:.2e}")

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # M_S scatter
    ax = axes[0]
    ax.scatter(df['depth'], df['m_s'], c='gray', alpha=0.1, s=5, rasterized=True)
    ax.scatter(site_data['depth'], site_data['m_s'], c='darkorange', s=100,
               edgecolors='black', linewidths=1.5, zorder=10, label='Site means')

    # Regression on site means
    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['m_s'])
    x_fit = np.linspace(df['depth'].min(), df['depth'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='black', linewidth=2, alpha=0.7)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Depth (μm)', fontsize=14)
    ax.set_ylabel('M_S (a.u.)', fontsize=14)
    ax.set_title(f'M_S vs Depth (All ROIs)\nROI-level r = {r_ms_roi:.3f}', fontsize=12)
    ax.legend(loc='upper right')
    ax.set_ylim(np.percentile(df['m_s'], 1), np.percentile(df['m_s'], 99))

    # M_C scatter
    ax = axes[1]
    ax.scatter(df['depth'], df['m_c'], c='gray', alpha=0.1, s=5, rasterized=True)
    ax.scatter(site_data['depth'], site_data['m_c'], c='darkorange', s=100,
               edgecolors='black', linewidths=1.5, zorder=10, label='Site means')

    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['m_c'])
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='black', linewidth=2, alpha=0.7)

    ax.set_xlabel('Depth (μm)', fontsize=14)
    ax.set_ylabel('M_C (r-value)', fontsize=14)
    ax.set_title(f'M_C vs Depth (All ROIs)\nROI-level r = {r_mc_roi:.3f}', fontsize=12)
    ax.legend(loc='upper left')
    ax.set_ylim(np.percentile(df['m_c'], 1), np.percentile(df['m_c'], 99))

    plt.tight_layout()
    plt.savefig(out_dir / 'roi_level_scatter.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'roi_level_scatter.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'roi_level_scatter.png'}")

    RESULTS['roi_scatter'] = {
        'r_ms_roi': r_ms_roi, 'p_ms_roi': p_ms_roi,
        'r_mc_roi': r_mc_roi, 'p_mc_roi': p_mc_roi,
        'n_rois': len(df), 'n_sites': len(site_data)
    }
    return RESULTS['roi_scatter']

# ============================================================
# MAIN
# ============================================================

def main():
    print("="*60)
    print("XORI SUPPLEMENTARY ANALYSIS")
    print("="*60)

    # Load all data
    print("\nLoading data...")
    site_depths = load_site_depths()
    print(f"  Loaded {len(site_depths)} site depths")

    osi_df = load_osi_data()
    print(f"  Loaded {len(osi_df)} OSI records")

    stat_df = load_roi_stat()
    print(f"  Loaded {len(stat_df)} ROI stat records")

    metric_df = load_metric_data(site_depths)
    print(f"  Loaded {len(metric_df)} metric records")

    # Merge all data
    df = merge_all_data(metric_df, osi_df, stat_df)
    print(f"\nMerged dataset: {len(df)} ROIs with complete data")

    # Run all analyses
    run_subpopulation_analysis(df)
    run_partial_correlation_analysis(df)
    run_within_site_analysis(df)
    run_ms_mc_analysis(df)
    run_morphology_analysis(df)
    run_bootstrap_analysis(df)
    run_roi_scatter_analysis(df)

    # Save results summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nAll outputs saved to: {OUTPUT_DIR}")

    return RESULTS

if __name__ == '__main__':
    results = main()
