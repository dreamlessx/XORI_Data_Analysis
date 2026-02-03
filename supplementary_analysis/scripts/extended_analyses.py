#!/usr/bin/env python3
"""
Extended Analyses for XORI Project
- Halfwidth (tuning width) vs depth
- LHI vs depth
- Spatial frequency preference vs depth
- Multi-metric correlation matrix
- Comprehensive depth profile
"""

import os
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
OUTPUT_DIR = BASE_DIR / "supplementary_analysis" / "outputs"

# ============================================================
# DATA LOADING
# ============================================================

def load_site_depths():
    """Load site depth mapping"""
    depths = {}
    with open(RAW_DATA / "site_depth.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                depths[parts[0]] = float(parts[1])
    return depths

def load_halfwidth_data():
    """Load halfwidth (tuning width) data"""
    data = []
    with open(RAW_DATA / "roi_hw_orth.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 5:
                try:
                    data.append({
                        'site': int(float(parts[0])),
                        'roi': int(float(parts[1])),
                        'hw_raw': float(parts[2]),
                        'hw_norm': float(parts[3]),
                        'orth_norm': float(parts[4])
                    })
                except: pass
    return pd.DataFrame(data)

def load_lhi_data():
    """Load LHI data"""
    data = []
    with open(RAW_DATA / "roi_lhi.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    data.append({
                        'site': int(float(parts[0])),
                        'roi': int(float(parts[1])),
                        'lhi2': float(parts[2]),
                        'lhi3': float(parts[3])
                    })
                except: pass
    return pd.DataFrame(data)

def load_sf_data():
    """Load spatial frequency preference data"""
    data = []
    with open(RAW_DATA / "roi_sf.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    data.append({
                        'site': int(float(parts[0])),
                        'roi': int(float(parts[1])),
                        'direction': float(parts[2]),
                        'sf': float(parts[3])
                    })
                except: pass
    return pd.DataFrame(data)

def load_osi_data():
    """Load OSI data"""
    data = []
    with open(RAW_DATA / "roi_osi.txt", 'r') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 8:
                try:
                    data.append({
                        'site': int(float(parts[0])),
                        'roi': int(float(parts[1])),
                        'bsf': float(parts[2]),
                        'bdr': float(parts[3]),
                        'dsi': float(parts[4]),
                        'osi': float(parts[5]),
                        'dcv': float(parts[6]),
                        'ocv': float(parts[7])
                    })
                except: pass
    return pd.DataFrame(data)

def load_metric_data(site_depths):
    """Load metric data"""
    all_data = []
    for mf in METRIC_DATA.glob("metrics_site*.txt"):
        site_name = mf.stem.replace("metrics_", "")
        site_num = int(site_name.replace("site", ""))
        if site_name not in site_depths: continue
        depth = site_depths[site_name]

        with open(mf, 'r') as f:
            for line in f.readlines()[2:]:
                parts = line.strip().split()
                if len(parts) >= 8:
                    try:
                        all_data.append({
                            'site': site_num, 'site_name': site_name,
                            'roi': int(parts[0]), 'depth': depth,
                            'm_s': float(parts[1]), 'm_c': float(parts[2]),
                            'snr_g': float(parts[3]), 'snr_p': float(parts[4]),
                            'm_s_norm': float(parts[5]), 'm_s_ratio': float(parts[6]),
                            'm_x': float(parts[7])
                        })
                    except: pass
    return pd.DataFrame(all_data)

def merge_all_data(site_depths):
    """Load and merge all datasets"""
    metric_df = load_metric_data(site_depths)
    hw_df = load_halfwidth_data()
    lhi_df = load_lhi_data()
    sf_df = load_sf_data()
    osi_df = load_osi_data()

    # Merge
    df = metric_df.merge(hw_df, on=['site', 'roi'], how='left')
    df = df.merge(lhi_df, on=['site', 'roi'], how='left')
    df = df.merge(sf_df, on=['site', 'roi'], how='left')
    df = df.merge(osi_df, on=['site', 'roi'], how='left')

    return df

# ============================================================
# ANALYSIS 1: HALFWIDTH VS DEPTH
# ============================================================

def run_halfwidth_analysis(df):
    """Analyze tuning width (halfwidth) vs depth"""
    print("\n" + "="*60)
    print("ANALYSIS: HALFWIDTH (TUNING WIDTH) VS DEPTH")
    print("="*60)

    out_dir = OUTPUT_DIR / "halfwidth_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level analysis
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'hw_raw': 'mean',
        'hw_norm': 'mean',
        'orth_norm': 'mean',
        'm_s': 'mean',
        'm_c': 'mean'
    }).reset_index().dropna()

    results = {}

    # Halfwidth vs depth
    r_hw_raw, p_hw_raw = pearsonr(site_data['depth'], site_data['hw_raw'])
    r_hw_norm, p_hw_norm = pearsonr(site_data['depth'], site_data['hw_norm'])
    r_orth, p_orth = pearsonr(site_data['depth'], site_data['orth_norm'])

    print(f"\n--- Halfwidth vs Depth ---")
    print(f"HW_raw ~ Depth: r = {r_hw_raw:+.3f}, p = {p_hw_raw:.4f}")
    print(f"HW_norm ~ Depth: r = {r_hw_norm:+.3f}, p = {p_hw_norm:.4f}")
    print(f"Orth_norm ~ Depth: r = {r_orth:+.3f}, p = {p_orth:.4f}")

    results['hw_raw_depth'] = {'r': r_hw_raw, 'p': p_hw_raw}
    results['hw_norm_depth'] = {'r': r_hw_norm, 'p': p_hw_norm}
    results['orth_depth'] = {'r': r_orth, 'p': p_orth}

    # Halfwidth vs M_S/M_C
    r_hw_ms, p_hw_ms = pearsonr(site_data['hw_raw'], site_data['m_s'])
    r_hw_mc, p_hw_mc = pearsonr(site_data['hw_raw'], site_data['m_c'])

    print(f"\n--- Halfwidth vs Metrics ---")
    print(f"HW_raw ~ M_S: r = {r_hw_ms:+.3f}, p = {p_hw_ms:.4f}")
    print(f"HW_raw ~ M_C: r = {r_hw_mc:+.3f}, p = {p_hw_mc:.4f}")

    results['hw_ms'] = {'r': r_hw_ms, 'p': p_hw_ms}
    results['hw_mc'] = {'r': r_hw_mc, 'p': p_hw_mc}

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Row 1: Halfwidth metrics vs Depth
    for idx, (col, label) in enumerate([('hw_raw', 'Raw Halfwidth (°)'),
                                         ('hw_norm', 'Normalized Halfwidth'),
                                         ('orth_norm', 'Orthogonal Norm')]):
        ax = axes[0, idx]
        ax.scatter(site_data['depth'], site_data[col], c='steelblue', s=80,
                   edgecolors='black', alpha=0.7)
        r, p = pearsonr(site_data['depth'], site_data[col])
        slope, intercept, _, _, _ = linregress(site_data['depth'], site_data[col])
        x_fit = np.linspace(site_data['depth'].min(), site_data['depth'].max(), 100)
        ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
        ax.set_xlabel('Depth (μm)', fontsize=12)
        ax.set_ylabel(label, fontsize=12)
        sig = "*" if p < 0.05 else ""
        ax.set_title(f'r = {r:.3f}, p = {p:.3f} {sig}', fontsize=11)

    # Row 2: Halfwidth vs M_S and M_C
    ax = axes[1, 0]
    ax.scatter(site_data['hw_raw'], site_data['m_s'], c='darkorange', s=80,
               edgecolors='black', alpha=0.7)
    r, p = pearsonr(site_data['hw_raw'], site_data['m_s'])
    slope, intercept, _, _, _ = linregress(site_data['hw_raw'], site_data['m_s'])
    x_fit = np.linspace(site_data['hw_raw'].min(), site_data['hw_raw'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Raw Halfwidth (°)', fontsize=12)
    ax.set_ylabel('M_S', fontsize=12)
    sig = "*" if p < 0.05 else ""
    ax.set_title(f'HW vs M_S: r = {r:.3f}, p = {p:.3f} {sig}', fontsize=11)

    ax = axes[1, 1]
    ax.scatter(site_data['hw_raw'], site_data['m_c'], c='darkorange', s=80,
               edgecolors='black', alpha=0.7)
    r, p = pearsonr(site_data['hw_raw'], site_data['m_c'])
    slope, intercept, _, _, _ = linregress(site_data['hw_raw'], site_data['m_c'])
    x_fit = np.linspace(site_data['hw_raw'].min(), site_data['hw_raw'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Raw Halfwidth (°)', fontsize=12)
    ax.set_ylabel('M_C', fontsize=12)
    sig = "*" if p < 0.05 else ""
    ax.set_title(f'HW vs M_C: r = {r:.3f}, p = {p:.3f} {sig}', fontsize=11)

    # Summary text
    ax = axes[1, 2]
    ax.axis('off')
    summary_text = """
    HALFWIDTH SUMMARY

    Halfwidth vs Depth:
    • Raw HW:  r = {:.3f}, p = {:.4f}
    • Norm HW: r = {:.3f}, p = {:.4f}
    • Orth:    r = {:.3f}, p = {:.4f}

    Halfwidth vs Metrics:
    • HW ~ M_S: r = {:.3f}, p = {:.4f}
    • HW ~ M_C: r = {:.3f}, p = {:.4f}
    """.format(r_hw_raw, p_hw_raw, r_hw_norm, p_hw_norm, r_orth, p_orth,
               r_hw_ms, p_hw_ms, r_hw_mc, p_hw_mc)
    ax.text(0.1, 0.5, summary_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='center', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(out_dir / 'halfwidth_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'halfwidth_analysis.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'halfwidth_analysis.png'}")

    return results

# ============================================================
# ANALYSIS 2: LHI VS DEPTH
# ============================================================

def run_lhi_analysis(df):
    """Analyze LHI vs depth"""
    print("\n" + "="*60)
    print("ANALYSIS: LHI VS DEPTH")
    print("="*60)

    out_dir = OUTPUT_DIR / "lhi_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level analysis
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'lhi2': 'mean',
        'lhi3': 'mean',
        'm_s': 'mean',
        'm_c': 'mean'
    }).reset_index().dropna()

    results = {}

    # LHI vs depth
    r_lhi2, p_lhi2 = pearsonr(site_data['depth'], site_data['lhi2'])
    r_lhi3, p_lhi3 = pearsonr(site_data['depth'], site_data['lhi3'])

    print(f"\n--- LHI vs Depth ---")
    print(f"LHI2 ~ Depth: r = {r_lhi2:+.3f}, p = {p_lhi2:.4f}")
    print(f"LHI3 ~ Depth: r = {r_lhi3:+.3f}, p = {p_lhi3:.4f}")

    results['lhi2_depth'] = {'r': r_lhi2, 'p': p_lhi2}
    results['lhi3_depth'] = {'r': r_lhi3, 'p': p_lhi3}

    # LHI vs M_S/M_C
    r_lhi2_ms, p_lhi2_ms = pearsonr(site_data['lhi2'], site_data['m_s'])
    r_lhi2_mc, p_lhi2_mc = pearsonr(site_data['lhi2'], site_data['m_c'])

    print(f"\n--- LHI vs Metrics ---")
    print(f"LHI2 ~ M_S: r = {r_lhi2_ms:+.3f}, p = {p_lhi2_ms:.4f}")
    print(f"LHI2 ~ M_C: r = {r_lhi2_mc:+.3f}, p = {p_lhi2_mc:.4f}")

    results['lhi2_ms'] = {'r': r_lhi2_ms, 'p': p_lhi2_ms}
    results['lhi2_mc'] = {'r': r_lhi2_mc, 'p': p_lhi2_mc}

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # LHI2 vs Depth
    ax = axes[0, 0]
    ax.scatter(site_data['depth'], site_data['lhi2'], c='purple', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['lhi2'])
    x_fit = np.linspace(site_data['depth'].min(), site_data['depth'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('LHI2', fontsize=12)
    sig = "*" if p_lhi2 < 0.05 else ""
    ax.set_title(f'LHI2 vs Depth: r = {r_lhi2:.3f}, p = {p_lhi2:.3f} {sig}', fontsize=11)

    # LHI3 vs Depth
    ax = axes[0, 1]
    ax.scatter(site_data['depth'], site_data['lhi3'], c='purple', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['lhi3'])
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('LHI3', fontsize=12)
    sig = "*" if p_lhi3 < 0.05 else ""
    ax.set_title(f'LHI3 vs Depth: r = {r_lhi3:.3f}, p = {p_lhi3:.3f} {sig}', fontsize=11)

    # LHI2 vs M_S
    ax = axes[1, 0]
    ax.scatter(site_data['lhi2'], site_data['m_s'], c='green', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['lhi2'], site_data['m_s'])
    x_fit = np.linspace(site_data['lhi2'].min(), site_data['lhi2'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('LHI2', fontsize=12)
    ax.set_ylabel('M_S', fontsize=12)
    sig = "*" if p_lhi2_ms < 0.05 else ""
    ax.set_title(f'LHI2 vs M_S: r = {r_lhi2_ms:.3f}, p = {p_lhi2_ms:.3f} {sig}', fontsize=11)

    # LHI2 vs M_C
    ax = axes[1, 1]
    ax.scatter(site_data['lhi2'], site_data['m_c'], c='green', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['lhi2'], site_data['m_c'])
    x_fit = np.linspace(site_data['lhi2'].min(), site_data['lhi2'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('LHI2', fontsize=12)
    ax.set_ylabel('M_C', fontsize=12)
    sig = "*" if p_lhi2_mc < 0.05 else ""
    ax.set_title(f'LHI2 vs M_C: r = {r_lhi2_mc:.3f}, p = {p_lhi2_mc:.3f} {sig}', fontsize=11)

    plt.tight_layout()
    plt.savefig(out_dir / 'lhi_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'lhi_analysis.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'lhi_analysis.png'}")

    return results

# ============================================================
# ANALYSIS 3: SPATIAL FREQUENCY VS DEPTH
# ============================================================

def run_sf_analysis(df):
    """Analyze spatial frequency preference vs depth"""
    print("\n" + "="*60)
    print("ANALYSIS: SPATIAL FREQUENCY VS DEPTH")
    print("="*60)

    out_dir = OUTPUT_DIR / "sf_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level analysis
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'sf': 'mean',
        'bsf': 'mean',  # from OSI data
        'm_s': 'mean',
        'm_c': 'mean',
        'osi': 'mean'
    }).reset_index().dropna()

    results = {}

    # SF vs depth
    r_sf, p_sf = pearsonr(site_data['depth'], site_data['sf'])
    r_bsf, p_bsf = pearsonr(site_data['depth'], site_data['bsf'])

    print(f"\n--- Spatial Frequency vs Depth ---")
    print(f"SF ~ Depth: r = {r_sf:+.3f}, p = {p_sf:.4f}")
    print(f"BSF ~ Depth: r = {r_bsf:+.3f}, p = {p_bsf:.4f}")

    results['sf_depth'] = {'r': r_sf, 'p': p_sf}
    results['bsf_depth'] = {'r': r_bsf, 'p': p_bsf}

    # SF vs M_S/M_C
    r_sf_ms, p_sf_ms = pearsonr(site_data['sf'], site_data['m_s'])
    r_sf_mc, p_sf_mc = pearsonr(site_data['sf'], site_data['m_c'])

    print(f"\n--- SF vs Metrics ---")
    print(f"SF ~ M_S: r = {r_sf_ms:+.3f}, p = {p_sf_ms:.4f}")
    print(f"SF ~ M_C: r = {r_sf_mc:+.3f}, p = {p_sf_mc:.4f}")

    results['sf_ms'] = {'r': r_sf_ms, 'p': p_sf_ms}
    results['sf_mc'] = {'r': r_sf_mc, 'p': p_sf_mc}

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # SF vs Depth
    ax = axes[0, 0]
    ax.scatter(site_data['depth'], site_data['sf'], c='teal', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['sf'])
    x_fit = np.linspace(site_data['depth'].min(), site_data['depth'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('Spatial Frequency (cpd)', fontsize=12)
    sig = "*" if p_sf < 0.05 else ""
    ax.set_title(f'SF vs Depth: r = {r_sf:.3f}, p = {p_sf:.3f} {sig}', fontsize=11)

    # OSI vs Depth
    r_osi, p_osi = pearsonr(site_data['depth'], site_data['osi'])
    ax = axes[0, 1]
    ax.scatter(site_data['depth'], site_data['osi'], c='coral', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['depth'], site_data['osi'])
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Depth (μm)', fontsize=12)
    ax.set_ylabel('OSI', fontsize=12)
    sig = "*" if p_osi < 0.05 else ""
    ax.set_title(f'OSI vs Depth: r = {r_osi:.3f}, p = {p_osi:.3f} {sig}', fontsize=11)
    results['osi_depth'] = {'r': r_osi, 'p': p_osi}

    # SF vs M_S
    ax = axes[1, 0]
    ax.scatter(site_data['sf'], site_data['m_s'], c='orange', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['sf'], site_data['m_s'])
    x_fit = np.linspace(site_data['sf'].min(), site_data['sf'].max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Spatial Frequency (cpd)', fontsize=12)
    ax.set_ylabel('M_S', fontsize=12)
    sig = "*" if p_sf_ms < 0.05 else ""
    ax.set_title(f'SF vs M_S: r = {r_sf_ms:.3f}, p = {p_sf_ms:.3f} {sig}', fontsize=11)

    # SF vs M_C
    ax = axes[1, 1]
    ax.scatter(site_data['sf'], site_data['m_c'], c='orange', s=80,
               edgecolors='black', alpha=0.7)
    slope, intercept, _, _, _ = linregress(site_data['sf'], site_data['m_c'])
    ax.plot(x_fit, slope * x_fit + intercept, '--', color='red', linewidth=2)
    ax.set_xlabel('Spatial Frequency (cpd)', fontsize=12)
    ax.set_ylabel('M_C', fontsize=12)
    sig = "*" if p_sf_mc < 0.05 else ""
    ax.set_title(f'SF vs M_C: r = {r_sf_mc:.3f}, p = {p_sf_mc:.3f} {sig}', fontsize=11)

    plt.tight_layout()
    plt.savefig(out_dir / 'sf_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'sf_analysis.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'sf_analysis.png'}")

    return results

# ============================================================
# ANALYSIS 4: MULTI-METRIC CORRELATION MATRIX
# ============================================================

def run_correlation_matrix(df):
    """Create comprehensive correlation matrix"""
    print("\n" + "="*60)
    print("ANALYSIS: MULTI-METRIC CORRELATION MATRIX")
    print("="*60)

    out_dir = OUTPUT_DIR / "correlation_matrix"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'm_x': 'mean',
        'snr_g': 'mean',
        'hw_raw': 'mean',
        'lhi2': 'mean',
        'sf': 'mean',
        'osi': 'mean'
    }).reset_index().dropna()

    # Select columns for correlation
    cols = ['depth', 'm_s', 'm_c', 'm_x', 'snr_g', 'hw_raw', 'lhi2', 'sf', 'osi']
    labels = ['Depth', 'M_S', 'M_C', 'M_X', 'SNR', 'Halfwidth', 'LHI', 'SF', 'OSI']

    # Compute correlation matrix
    n = len(cols)
    corr_matrix = np.zeros((n, n))
    pval_matrix = np.zeros((n, n))

    for i, col1 in enumerate(cols):
        for j, col2 in enumerate(cols):
            if i == j:
                corr_matrix[i, j] = 1.0
                pval_matrix[i, j] = 0.0
            else:
                r, p = pearsonr(site_data[col1], site_data[col2])
                corr_matrix[i, j] = r
                pval_matrix[i, j] = p

    # Print significant correlations
    print("\n--- Significant Correlations (p < 0.05) ---")
    for i in range(n):
        for j in range(i+1, n):
            if pval_matrix[i, j] < 0.05:
                stars = "***" if pval_matrix[i, j] < 0.001 else "**" if pval_matrix[i, j] < 0.01 else "*"
                print(f"{labels[i]:12s} ~ {labels[j]:12s}: r = {corr_matrix[i,j]:+.3f}, p = {pval_matrix[i,j]:.4f} {stars}")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Create heatmap
    im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Pearson r', fontsize=12)

    # Set ticks
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(labels, fontsize=11, rotation=45, ha='right')
    ax.set_yticklabels(labels, fontsize=11)

    # Add correlation values as text
    for i in range(n):
        for j in range(n):
            r = corr_matrix[i, j]
            p = pval_matrix[i, j]
            color = 'white' if abs(r) > 0.5 else 'black'
            text = f'{r:.2f}'
            if i != j and p < 0.05:
                text += '*'
            if i != j and p < 0.01:
                text += '*'
            if i != j and p < 0.001:
                text += '*'
            ax.text(j, i, text, ha='center', va='center', color=color, fontsize=9)

    ax.set_title('Multi-Metric Correlation Matrix\n(Site-Level, * p<.05, ** p<.01, *** p<.001)', fontsize=14)

    plt.tight_layout()
    plt.savefig(out_dir / 'correlation_matrix.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'correlation_matrix.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'correlation_matrix.png'}")

    return {'corr': corr_matrix, 'pval': pval_matrix, 'labels': labels}

# ============================================================
# ANALYSIS 5: COMPREHENSIVE DEPTH PROFILE
# ============================================================

def run_depth_profile(df):
    """Create comprehensive depth profile showing all metrics"""
    print("\n" + "="*60)
    print("ANALYSIS: COMPREHENSIVE DEPTH PROFILE")
    print("="*60)

    out_dir = OUTPUT_DIR / "depth_profile"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'm_x': 'mean',
        'snr_g': 'mean',
        'hw_raw': 'mean',
        'lhi2': 'mean',
        'sf': 'mean',
        'osi': 'mean'
    }).reset_index().dropna()

    site_data = site_data.sort_values('depth')

    # Z-score normalize for comparison
    for col in ['m_s', 'm_c', 'm_x', 'snr_g', 'hw_raw', 'lhi2', 'sf', 'osi']:
        site_data[f'{col}_z'] = (site_data[col] - site_data[col].mean()) / site_data[col].std()

    # Compute all correlations with depth
    results = {}
    metrics = ['m_s', 'm_c', 'm_x', 'snr_g', 'hw_raw', 'lhi2', 'sf', 'osi']
    labels = ['M_S', 'M_C', 'M_X', 'SNR', 'Halfwidth', 'LHI2', 'SF', 'OSI']

    print("\n--- All Metrics vs Depth ---")
    for metric, label in zip(metrics, labels):
        r, p = pearsonr(site_data['depth'], site_data[metric])
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"{label:12s} ~ Depth: r = {r:+.3f}, p = {p:.4f} {sig}")
        results[metric] = {'r': r, 'p': p}

    # Create figure
    fig = plt.figure(figsize=(16, 12))

    # Panel A: All metrics vs depth (z-scored)
    ax1 = fig.add_subplot(2, 2, 1)
    colors = plt.cm.tab10(np.linspace(0, 1, len(metrics)))

    for i, (metric, label) in enumerate(zip(metrics, labels)):
        r, p = results[metric]['r'], results[metric]['p']
        style = '-' if p < 0.05 else '--'
        alpha = 1.0 if p < 0.05 else 0.5
        ax1.plot(site_data['depth'], site_data[f'{metric}_z'], style,
                 color=colors[i], linewidth=2, alpha=alpha,
                 label=f'{label} (r={r:.2f}{"*" if p<0.05 else ""})')

    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax1.set_xlabel('Depth (μm)', fontsize=12)
    ax1.set_ylabel('Z-score', fontsize=12)
    ax1.set_title('A. All Metrics vs Depth (Z-scored)', fontsize=14)
    ax1.legend(loc='upper right', fontsize=9, ncol=2)

    # Panel B: Forest plot of correlations
    ax2 = fig.add_subplot(2, 2, 2)
    r_vals = [results[m]['r'] for m in metrics]
    p_vals = [results[m]['p'] for m in metrics]
    colors_bar = ['forestgreen' if p < 0.05 else 'indianred' for p in p_vals]

    y_pos = np.arange(len(metrics))
    ax2.barh(y_pos, r_vals, color=colors_bar, alpha=0.8, edgecolor='black')
    ax2.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=11)
    ax2.set_xlabel('Correlation with Depth (r)', fontsize=12)
    ax2.set_title('B. Depth Correlations\n(green = p < 0.05)', fontsize=14)
    ax2.set_xlim(-1, 1)

    for i, (r, p) in enumerate(zip(r_vals, p_vals)):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        ax2.text(0.95 if r < 0 else -0.95, i, f'{r:.2f}{sig}',
                 va='center', ha='left' if r < 0 else 'right', fontsize=10)

    # Panel C: Depth bins comparison
    ax3 = fig.add_subplot(2, 2, 3)

    # Bin into superficial/middle/deep
    bins = pd.cut(site_data['depth'], bins=3, labels=['Superficial', 'Middle', 'Deep'])
    site_data['depth_bin'] = bins

    bin_means = site_data.groupby('depth_bin').agg({
        'm_s': 'mean', 'm_c': 'mean', 'hw_raw': 'mean', 'sf': 'mean', 'osi': 'mean'
    }).reset_index()

    x = np.arange(3)
    width = 0.15

    metrics_plot = ['m_s', 'm_c', 'hw_raw', 'sf', 'osi']
    labels_plot = ['M_S', 'M_C', 'HW', 'SF', 'OSI']
    colors_plot = ['#2E86AB', '#E94F37', '#8338EC', '#06D6A0', '#FF9F1C']

    # Normalize for comparison
    for i, (m, l, c) in enumerate(zip(metrics_plot, labels_plot, colors_plot)):
        vals = bin_means[m].values
        vals_norm = (vals - vals.mean()) / (vals.std() + 1e-9)
        ax3.bar(x + i*width - 2*width, vals_norm, width, label=l, color=c, edgecolor='black')

    ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax3.set_xticks(x)
    ax3.set_xticklabels(['Superficial\n(140-266μm)', 'Middle\n(266-392μm)', 'Deep\n(392-518μm)'], fontsize=10)
    ax3.set_ylabel('Normalized Value', fontsize=12)
    ax3.set_title('C. Metrics by Depth Bin', fontsize=14)
    ax3.legend(loc='upper right', fontsize=9)

    # Panel D: Summary statistics table
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    # Create summary table
    table_data = []
    for metric, label in zip(metrics, labels):
        r, p = results[metric]['r'], results[metric]['p']
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        table_data.append([label, f'{r:+.3f}', f'{p:.4f}', sig])

    table = ax4.table(cellText=table_data,
                      colLabels=['Metric', 'r', 'p-value', 'Sig'],
                      loc='center', cellLoc='center',
                      colWidths=[0.25, 0.2, 0.25, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)

    # Color code significant rows
    for i, p in enumerate(p_vals):
        color = '#c8e6c9' if p < 0.05 else '#ffcdd2'
        for j in range(4):
            table[(i+1, j)].set_facecolor(color)

    ax4.set_title('D. Summary Statistics', fontsize=14, pad=20)

    plt.tight_layout()
    plt.savefig(out_dir / 'depth_profile.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'depth_profile.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'depth_profile.png'}")

    return results

# ============================================================
# ANALYSIS 6: MEDIATION ANALYSIS
# ============================================================

def run_mediation_analysis(df):
    """Test if other metrics mediate the M_S ~ Depth relationship"""
    print("\n" + "="*60)
    print("ANALYSIS: MEDIATION ANALYSIS")
    print("="*60)

    out_dir = OUTPUT_DIR / "mediation_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Site-level data
    site_data = df.groupby('site_name').agg({
        'depth': 'first',
        'm_s': 'mean',
        'm_c': 'mean',
        'hw_raw': 'mean',
        'lhi2': 'mean',
        'sf': 'mean',
        'osi': 'mean',
        'snr_g': 'mean'
    }).reset_index().dropna()

    # Simple M_S ~ Depth
    r_simple, p_simple = pearsonr(site_data['depth'], site_data['m_s'])
    print(f"\nSimple: M_S ~ Depth: r = {r_simple:.3f}, p = {p_simple:.4f}")

    # Partial correlation function
    def partial_corr(x, y, covariates):
        from numpy.linalg import lstsq
        X_cov = np.column_stack([np.ones(len(covariates)), covariates])
        beta_x, _, _, _ = lstsq(X_cov, x, rcond=None)
        x_resid = x - X_cov @ beta_x
        beta_y, _, _, _ = lstsq(X_cov, y, rcond=None)
        y_resid = y - X_cov @ beta_y
        return pearsonr(x_resid, y_resid)

    # Test each potential mediator
    mediators = ['hw_raw', 'lhi2', 'sf', 'osi', 'snr_g']
    labels = ['Halfwidth', 'LHI', 'SF', 'OSI', 'SNR']

    results = {'simple': {'r': r_simple, 'p': p_simple}}

    print("\n--- Controlling for potential mediators ---")
    for med, label in zip(mediators, labels):
        cov = site_data[med].values.reshape(-1, 1)
        r_part, p_part = partial_corr(site_data['depth'].values, site_data['m_s'].values, cov)
        change = ((r_simple - r_part) / r_simple) * 100
        print(f"M_S ~ Depth | {label:10s}: r = {r_part:+.3f}, p = {p_part:.4f} ({change:+.1f}% change)")
        results[med] = {'r': r_part, 'p': p_part, 'change': change}

    # Create visualization
    fig, ax = plt.subplots(figsize=(10, 6))

    labels_plot = ['No control'] + labels
    r_vals = [r_simple] + [results[m]['r'] for m in mediators]
    p_vals = [p_simple] + [results[m]['p'] for m in mediators]
    changes = [0] + [results[m]['change'] for m in mediators]

    colors = ['forestgreen' if p < 0.05 else 'indianred' for p in p_vals]

    y_pos = np.arange(len(labels_plot))
    bars = ax.barh(y_pos, r_vals, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels_plot, fontsize=12)
    ax.set_xlabel('Partial Correlation (r)', fontsize=14)
    ax.set_title('M_S ~ Depth: Mediation Analysis\n(green = p < 0.05)', fontsize=14)
    ax.set_xlim(-1, 0.2)

    for i, (r, p, c) in enumerate(zip(r_vals, p_vals, changes)):
        text = f'r={r:.2f}'
        if i > 0:
            text += f' ({c:+.0f}%)'
        ax.text(0.05, i, text, va='center', ha='left', fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_dir / 'mediation_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig(out_dir / 'mediation_analysis.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"\n[SAVED] {out_dir / 'mediation_analysis.png'}")

    return results

# ============================================================
# MAIN
# ============================================================

def main():
    print("="*60)
    print("XORI EXTENDED ANALYSES")
    print("="*60)

    # Load data
    print("\nLoading all data...")
    site_depths = load_site_depths()
    df = merge_all_data(site_depths)
    print(f"Loaded {len(df)} ROIs with all metrics")

    # Run all analyses
    hw_results = run_halfwidth_analysis(df)
    lhi_results = run_lhi_analysis(df)
    sf_results = run_sf_analysis(df)
    corr_results = run_correlation_matrix(df)
    depth_results = run_depth_profile(df)
    med_results = run_mediation_analysis(df)

    print("\n" + "="*60)
    print("EXTENDED ANALYSES COMPLETE")
    print("="*60)

    # Summary
    print("\n" + "="*60)
    print("KEY FINDINGS SUMMARY")
    print("="*60)

    print("\nMetrics that correlate with Depth (p < 0.05):")
    for metric, label in [('m_s', 'M_S'), ('m_c', 'M_C'), ('snr_g', 'SNR')]:
        if metric in depth_results and depth_results[metric]['p'] < 0.05:
            print(f"  ✓ {label}: r = {depth_results[metric]['r']:.3f}")

    print("\nMetrics that do NOT correlate with Depth:")
    for metric, label in [('hw_raw', 'Halfwidth'), ('lhi2', 'LHI'), ('sf', 'SF'), ('osi', 'OSI')]:
        if metric in depth_results and depth_results[metric]['p'] >= 0.05:
            print(f"  ✗ {label}: r = {depth_results[metric]['r']:.3f}, p = {depth_results[metric]['p']:.3f}")

if __name__ == '__main__':
    main()
