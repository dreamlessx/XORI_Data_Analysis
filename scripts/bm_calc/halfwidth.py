#!/usr/bin/env python3
"""
halfwidth_analysis.py
Generates halfwidth analysis plots for all_roi metrics (using raw halfwidth values).

Structure:
data_halfwidth/
└── raw/
    └── all_roi/
        ├── metric_site/ (M_S, M_S_ratio, M_S_log, M_C, M_X vs halfwidth per site)
        ├── pearson_site/ (R-values vs depth)
        └── depth_site/ (halfwidth vs depth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
HW_ORTH_FILE = Path('./raw_data/bm_data/roi_hw_orth.txt')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_ALL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_halfwidth')

# ========= Load site depths =========
def load_site_depths(depth_file):
    """
    Read site_depth.txt with format:
    site        depth
    --------------------
    site038      140
    ...
    Returns: {site_name: depth_value}
    """
    depths = {}
    with open(depth_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header and separator
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            site = parts[0].strip()
            depth = float(parts[1])
            depths[site] = depth
    
    return depths

# ========= Load halfwidth data =========
def load_halfwidth_data(hw_orth_file):
    """
    Read roi_hw_orth.txt and return halfwidth data organized by site.
    Returns: {site_name: {roi_idx: hw_norm}}
    Only includes ROIs with valid (non-null) hw_norm values.
    """
    df = pd.read_csv(hw_orth_file, sep=r'\s+', skiprows=2, header=None,
                     names=['site', 'roi', 'hw_raw', 'hw_norm', 'orth_norm'])
    
    # Convert numeric columns, replacing 'null' with NaN
    for col in ['hw_raw', 'hw_norm', 'orth_norm']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Only keep rows with valid hw_norm
    df_clean = df.dropna(subset=['hw_norm'])
    
    # Organize by site
    site_hw_data = {}
    for _, row in df_clean.iterrows():
        site_num = int(row['site'])
        site_name = f"site{site_num:03d}"
        roi_idx = int(row['roi'])
        hw_norm = row['hw_norm']
        
        if site_name not in site_hw_data:
            site_hw_data[site_name] = {}
        
        site_hw_data[site_name][roi_idx] = hw_norm
    
    return site_hw_data

# ========= Load metrics from metric files =========
def load_metrics_from_file(metric_file):
    """
    Read metrics file and return dict: {roi: {M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_S_log, M_X}}
    M_S_log is log2(M_S_ratio) for positive values, None for negative values
    M_X is peak suppression metric
    """
    metrics = {}
    
    with open(metric_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header and separator (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 8:
            try:
                roi = int(parts[0])
                m_s_ratio = float(parts[6])
                
                # Calculate M_S_log: log2 of M_S_ratio, only if positive
                if m_s_ratio > 0:
                    m_s_log = np.log2(m_s_ratio)
                else:
                    m_s_log = None  # Mark as invalid for negative values
                
                metrics[roi] = {
                    'M_S': float(parts[1]),
                    'M_C': float(parts[2]),
                    'SNR_g': float(parts[3]),
                    'SNR_p': float(parts[4]),
                    'M_S_norm': float(parts[5]),
                    'M_S_ratio': m_s_ratio,
                    'M_S_log': m_s_log,
                    'M_X': float(parts[7])  # Peak suppression metric
                }
            except (ValueError, IndexError):
                continue
    
    return metrics

# ========= Plot metric vs halfwidth for one site =========
def plot_metric_vs_halfwidth_site(site, hw_data, metrics, metric_key, metric_label, 
                                   out_file):
    """
    Create scatter plot of metric vs halfwidth for one site.
    hw_data: dict of {roi_idx: hw_norm}
    metrics: dict of {roi: {metric_key: value}}
    """
    x_vals = []
    y_vals = []
    
    for roi_idx in sorted(metrics.keys()):
        if roi_idx in hw_data:
            hw = hw_data[roi_idx]
            y_val = metrics[roi_idx][metric_key]
            
            # Skip if M_S_log is None (negative M_S_ratio)
            if metric_key == 'M_S_log' and y_val is None:
                continue
            
            if hw > 0:
                x_vals.append(hw)
                y_vals.append(y_val)
    
    if len(x_vals) < 2:
        print(f"[SKIP] {site}: not enough data points for {metric_key}")
        return
    
    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(x_vals, y_vals, s=32, alpha=0.7)
    
    # Fit line
    slope, intercept, *_ = linregress(x_vals, y_vals)
    xx = np.linspace(np.min(x_vals), np.max(x_vals), 200)
    yy = intercept + slope * xx
    ax.plot(xx, yy, '--', color='black', linewidth=2)
    
    # Stats
    r_val, p_val = pearsonr(x_vals, y_vals)
    n = len(x_vals)
    txt = f"y = {slope:+.2f}x {intercept:+.2f}\nPearson r = {r_val:+.2f}, p = {p_val:.2g}, n = {n}"
    ax.text(0.05, 0.95, txt, transform=ax.transAxes, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), fontsize=12)
    
    # Labels
    x_label = "Half-width (°)"
    y_label = metric_label
    
    ax.set_title(f"{site}: {x_label} vs {y_label}", fontsize=14)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Add padding to axes
    x_range = np.max(x_vals) - np.min(x_vals)
    y_range = np.max(y_vals) - np.min(y_vals)
    x_padding = x_range * 0.05
    y_padding = y_range * 0.1
    ax.set_xlim(np.min(x_vals) - x_padding, np.max(x_vals) + x_padding)
    ax.set_ylim(np.min(y_vals) - y_padding, np.max(y_vals) + y_padding)
    
    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=150)
    plt.close()

# ========= Calculate R-values per site =========
def calculate_r_values_per_site(site_hw_data, metric_data_dir):
    """
    For each site, calculate Pearson r between halfwidth and each metric.
    Returns: {metric_key: {site: (r, p)}}
    """
    results = {
        'M_S': {},
        'M_S_norm': {},
        'M_S_ratio': {},
        'M_S_log': {},
        'M_C': {},
        'M_X': {}
    }
    
    for metric_file in sorted(metric_data_dir.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        if site_name not in site_hw_data:
            continue
        
        hw_data = site_hw_data[site_name]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            
            for roi_idx in sorted(metrics.keys()):
                if roi_idx in hw_data:
                    hw = hw_data[roi_idx]
                    y_val = metrics[roi_idx][metric_key]
                    
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    
                    if hw > 0:
                        x_vals.append(hw)
                        y_vals.append(y_val)
            
            if len(x_vals) >= 2:
                r, p = pearsonr(x_vals, y_vals)
                results[metric_key][site_name] = (r, p)
    
    return results

# ========= Plot R-values vs depth =========
def plot_r_vs_depth(r_values, site_depths, metric_key, metric_label, out_file):
    """
    Plot Pearson r-values vs depth for one metric.
    r_values: {site: (r, p)}
    """
    depths = []
    r_vals = []
    p_vals = []
    
    for site, (r, p) in r_values.items():
        if site in site_depths:
            depths.append(site_depths[site])
            r_vals.append(r)
            p_vals.append(p)
    
    if len(depths) == 0:
        print(f"[SKIP] No data for {metric_key} R vs depth")
        return
    
    depths = np.array(depths)
    r_vals = np.array(r_vals)
    p_vals = np.array(p_vals)
    
    sig_mask = p_vals < 0.05
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Significant points (filled)
    if np.any(sig_mask):
        ax.scatter(depths[sig_mask], r_vals[sig_mask],
                  s=150, c='blue', marker='o',
                  label='p < 0.05', edgecolors='black', linewidths=2)
    
    # Non-significant points (hollow)
    if np.any(~sig_mask):
        ax.scatter(depths[~sig_mask], r_vals[~sig_mask],
                  s=150, facecolors='none', edgecolors='blue',
                  marker='o', linewidths=2, label='p ≥ 0.05')
    
    # Zero line
    ax.axhline(y=0, color='black', linestyle='--', linewidth=2, label='r = 0')
    
    # Labels and title based on metric
    if metric_key == 'M_S':
        ylabel = 'S vs. Half-width (r-value)'
        title = f'Metric S (Half-width): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_C':
        ylabel = 'C vs. Half-width (r-value)'
        title = f'Metric C (Half-width): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_S_ratio':
        ylabel = 'S vs. Half-width (r-value)'
        title = f'Metric S (Half-width): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_S_log':
        ylabel = 'S vs. Half-width (r-value)'
        title = f'Metric S (Half-width): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_X':
        ylabel = 'X vs. Half-width (r-value)'
        title = f'Metric X (Half-width): r-value vs. Depth (N={len(depths)})'
    else:
        ylabel = f'{metric_label} vs. Half-width (r-value)'
        title = f'{metric_label} (Half-width): r-value vs. Depth (N={len(depths)})'
    
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_ylabel(ylabel, fontsize=20, fontweight='bold', labelpad=15)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Legend: explicitly not bold
    legend = ax.legend(loc='best', fontsize=12, frameon=True)
    for text in legend.get_texts():
        text.set_fontweight('normal')
    
    # Set y-axis limits to -0.5 to 0.4 for all metrics
    ax.set_ylim(-0.5, 0.4)
    ax.set_yticks(np.arange(-0.5, 0.41, 0.1))
    
    # Add padding to x-axis
    if len(depths) > 0:
        x_min, x_max = depths.min(), depths.max()
        x_range = x_max - x_min
        x_padding = x_range * 0.05
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
    
    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=150)
    plt.savefig(out_file.with_suffix('.svg'), format='svg')
    plt.close()
    
    print(f"[SAVED] {out_file}")

# ========= Plot halfwidth vs depth =========
def plot_halfwidth_vs_depth(site_hw_data, site_depths, out_file):
    """
    Create publication-quality halfwidth vs depth plot (raw values).
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site, hw_dict in site_hw_data.items():
        if site not in site_depths:
            continue
        
        depth = site_depths[site]
        
        # Get all halfwidth values for this site
        hw_values = [hw for hw in hw_dict.values() if hw > 0]
        
        if not hw_values:
            continue
        
        # Site mean and SEM (raw values)
        mean_hw = np.mean(hw_values)
        n = len(hw_values)
        if n > 1:
            sd_hw = np.std(hw_values, ddof=1)
            sem_hw = sd_hw / np.sqrt(n)
        else:
            sem_hw = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_hw)
        site_sems.append(sem_hw)
    
    site_means_depths = np.array(site_means_depths)
    site_means = np.array(site_means)
    site_sems = np.array(site_sems)
    
    # Regression on site means
    if len(site_means_depths) >= 2:
        slope, intercept, *_ = linregress(site_means_depths, site_means)
        r_val, p_val = pearsonr(site_means_depths, site_means)
        xx = np.linspace(site_means_depths.min(), site_means_depths.max(), 200)
        yy = intercept + slope * xx
    else:
        slope = intercept = r_val = p_val = None
        xx = yy = None
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Site means with SEM error bars (orange, prominent)
    ax.errorbar(site_means_depths, site_means, yerr=site_sems,
                fmt='o', markersize=8, color='darkorange',
                markeredgecolor='black', markeredgewidth=1.5,
                ecolor='black', elinewidth=1.5, capsize=4, capthick=1.5,
                alpha=1.0, zorder=3)
    
    # Regression line (half contrast)
    if xx is not None:
        ax.plot(xx, yy, linestyle='--', color='black', linewidth=2, alpha=0.4, zorder=2)
    
    # Stats annotation - TOP LEFT
    if slope is not None and r_val is not None:
        eq_line = f"y = {slope:.4f}x + {intercept:.2f}"
        stats_line = f"r = {r_val:+.3f}, p = {p_val:.3g}, N = {len(site_means_depths)}"
        stats_txt = f"{eq_line}\n{stats_line}"
        ax.text(0.05, 0.95, stats_txt,
                transform=ax.transAxes,
                fontsize=16, fontweight='normal',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='black', alpha=0.9, linewidth=1.5))
    
    # Labels - axis labels 20pt bold, no title, extra spacing with labelpad
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_ylabel('Half-width (°)', fontsize=20, fontweight='bold', labelpad=15)
    
    # Add padding to both axes
    if len(site_means_depths) > 0:
        x_min, x_max = site_means_depths.min(), site_means_depths.max()
        x_range = x_max - x_min
        x_padding = x_range * 0.05
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
    
    if len(site_means) > 0:
        y_min, y_max = site_means.min(), site_means.max()
        y_range = y_max - y_min
        y_padding = y_range * 0.1
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
    
    # Tick styling - increased font size
    ax.tick_params(axis='both', which='major', labelsize=16, width=1.5, length=6)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Grid for cleaner look
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.savefig(out_file.with_suffix('.svg'), format='svg', bbox_inches='tight')
    plt.close()
    
    print(f"[SAVED] {out_file}")

# ========= Main processing function =========
def process_halfwidth_analysis():
    """
    Process halfwidth analysis for all_roi data (raw values).
    """
    print(f"\n{'='*60}")
    print(f"Processing all_roi with raw halfwidth")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Load halfwidth data
    site_hw_data = load_halfwidth_data(HW_ORTH_FILE)
    print(f"[INFO] Loaded halfwidth data for {len(site_hw_data)} sites")
    
    # Create output directories
    output_dir = OUTPUT_BASE / 'raw' / 'all_roi'
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',
              pearson_site_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs halfwidth plots
    print(f"\n[INFO] Generating metric vs halfwidth plots for each site...")
    for metric_file in sorted(METRIC_DATA_ALL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        if site_name not in site_hw_data:
            print(f"[SKIP] {site_name}: no halfwidth data")
            continue
        
        hw_data = site_hw_data[site_name]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_S', 'M_S',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_halfwidth.png'
        )
        
        # M_S_norm (percentage of baseline)
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_S_norm', 'M_S_norm',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_halfwidth.png'
        )
        
        # M_S_ratio
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_S_ratio', 'M_S_ratio',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_halfwidth.png'
        )
        
        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_S_log', 'M_S_log',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_halfwidth.png'
        )
        
        # M_C
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_C', 'M_C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_halfwidth.png'
        )
        
        # M_X (peak suppression)
        plot_metric_vs_halfwidth_site(
            site_name, hw_data, metrics, 'M_X', 'M_X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_halfwidth.png'
        )
    
    # Calculate R-values and plot vs depth
    print(f"\n[INFO] Calculating R-values and plotting vs depth...")
    r_values = calculate_r_values_per_site(site_hw_data, METRIC_DATA_ALL)
    
    plot_r_vs_depth(r_values['M_S'], site_depths, 'M_S', 'M_S',
                   pearson_site_dir / 'm_s.png')
    plot_r_vs_depth(r_values['M_S_norm'], site_depths, 'M_S_norm', 'M_S_norm',
                   pearson_site_dir / 'm_s_norm.png')
    plot_r_vs_depth(r_values['M_S_ratio'], site_depths, 'M_S_ratio', 'M_S_ratio',
                   pearson_site_dir / 'm_s_r.png')
    plot_r_vs_depth(r_values['M_S_log'], site_depths, 'M_S_log', 'M_S_log',
                   pearson_site_dir / 'm_s_l.png')
    plot_r_vs_depth(r_values['M_C'], site_depths, 'M_C', 'M_C',
                   pearson_site_dir / 'm_c.png')
    plot_r_vs_depth(r_values['M_X'], site_depths, 'M_X', 'M_X',
                   pearson_site_dir / 'm_x.png')
    
    # Plot halfwidth vs depth
    print(f"\n[INFO] Plotting halfwidth vs depth...")
    plot_halfwidth_vs_depth(site_hw_data, site_depths,
                           depth_site_dir / 'halfwidth_vs_depth.png')

# ========= Main =========
def main():
    print("Starting halfwidth analysis...")
    
    # Process halfwidth analysis for all_roi only
    print("\n" + "="*60)
    print("RAW HALFWIDTH ANALYSIS - ALL ROI")
    print("="*60)
    
    process_halfwidth_analysis()
    
    print(f"\n{'='*60}")
    print("Halfwidth analysis complete!")
    print(f"Output directory: {OUTPUT_BASE / 'raw' / 'all_roi'}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()