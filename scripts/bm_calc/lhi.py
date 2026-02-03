#!/usr/bin/env python3
"""
lhi_analysis.py
Generates Local Homogeneity Index analysis plots for all_roi metrics.

Structure:
data_lhi/
├── 2d/          # LHI2 (2D local homogeneity)
│   ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_S_log, M_C, M_X vs LHI2 per site)
│   ├── pearson_site/ (R-values vs depth)
│   └── depth_site/ (LHI2 vs depth)
└── 3d/          # LHI3 (3D local homogeneity)
    ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_S_log, M_C, M_X vs LHI3 per site)
    ├── pearson_site/ (R-values vs depth)
    └── depth_site/ (LHI3 vs depth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
LHI_FILE = Path('./raw_data/bm_data/roi_lhi.txt')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_ALL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_lhi')

# ========= Load LHI data from roi_lhi.txt =========
def load_lhi_data(lhi_file):
    """
    Read roi_lhi.txt and return dicts for both LHI methods
    Format: site | ROI | LHI2 | LHI3
    Returns: 
        lhi2_data: {site: {roi: lhi2_value}} - 2D local homogeneity
        lhi3_data: {site: {roi: lhi3_value}} - 3D local homogeneity
    """
    with open(lhi_file, 'r') as f:
        lines = f.readlines()
    
    lhi2_data = {}
    lhi3_data = {}
    
    # Skip header and separator (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        
        try:
            site = int(float(parts[0]))
            roi = int(float(parts[1]))
            lhi2 = float(parts[2])
            lhi3 = float(parts[3])
            
            # LHI2 (2D)
            if site not in lhi2_data:
                lhi2_data[site] = {}
            lhi2_data[site][roi] = lhi2
            
            # LHI3 (3D)
            if site not in lhi3_data:
                lhi3_data[site] = {}
            lhi3_data[site][roi] = lhi3
            
        except (ValueError, IndexError):
            continue
    
    return lhi2_data, lhi3_data

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

# ========= Convert site number to site name =========
def site_num_to_name(site_num):
    """Convert site number (e.g., 2) to site name (e.g., 'site002')"""
    return f"site{int(site_num):03d}"

# ========= Plot metric vs LHI for one site =========
def plot_metric_vs_lhi_site(site_name, lhi_values, metrics, metric_key, 
                            metric_label, out_file, lhi_type='2d'):
    """
    Create scatter plot of metric vs LHI for one site.
    lhi_values: dict of {roi: lhi_value}
    metrics: dict of {roi: {metric_key: value}}
    lhi_type: '2d' or '3d'
    """
    # Match LHI to metrics by ROI
    x_vals = []
    y_vals = []
    
    for roi_idx in sorted(metrics.keys()):
        if roi_idx in lhi_values:
            y_val = metrics[roi_idx][metric_key]
            
            # Skip if M_S_log is None (negative M_S_ratio)
            if metric_key == 'M_S_log' and y_val is None:
                continue
            
            x_vals.append(lhi_values[roi_idx])
            y_vals.append(y_val)
    
    if len(x_vals) < 2:
        print(f"[SKIP] {site_name}: not enough data points for {metric_key}")
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
    lhi_label = "LHI2 (2D)" if lhi_type == '2d' else "LHI3 (3D)"
    
    if metric_key == 'M_S_norm':
        y_label = "M_S (% of baseline)"
    else:
        y_label = metric_label
    
    ax.set_title(f"{site_name}: {lhi_label} vs {y_label}", fontsize=14)
    ax.set_xlabel(lhi_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Add padding to axes
    x_range = np.max(x_vals) - np.min(x_vals)
    y_range = np.max(y_vals) - np.min(y_vals)
    x_padding = x_range * 0.05 if x_range > 0 else 0.05
    y_padding = y_range * 0.1 if y_range > 0 else 0.1
    ax.set_xlim(np.min(x_vals) - x_padding, np.max(x_vals) + x_padding)
    ax.set_ylim(np.min(y_vals) - y_padding, np.max(y_vals) + y_padding)
    
    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=150)
    plt.close()

# ========= Calculate R-values per site =========
def calculate_r_values_per_site(lhi_data, metric_data_dir):
    """
    For each site, calculate Pearson r between LHI and each metric.
    Returns: {metric_key: {site_name: (r, p)}}
    lhi_data: {site_num: {roi: lhi_value}}
    """
    results = {
        'M_S': {},
        'M_S_norm': {},
        'M_S_ratio': {},
        'M_S_log': {},
        'M_C': {},
        'M_X': {}  # Added M_X
    }
    
    for metric_file in sorted(metric_data_dir.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        # Extract site number from site_name (e.g., 'site002' -> 2)
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in lhi_data:
            continue
        
        lhi_values = lhi_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            
            for roi_idx in sorted(metrics.keys()):
                if roi_idx in lhi_values:
                    y_val = metrics[roi_idx][metric_key]
                    
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    
                    x_vals.append(lhi_values[roi_idx])
                    y_vals.append(y_val)
            
            if len(x_vals) >= 2:
                r, p = pearsonr(x_vals, y_vals)
                results[metric_key][site_name] = (r, p)
    
    return results

# ========= Plot R-values vs depth =========
def plot_r_vs_depth(r_values, site_depths, metric_key, metric_label, out_file, lhi_type='2d'):
    """
    Plot Pearson r-values vs depth for one metric.
    r_values: {site_name: (r, p)}
    lhi_type: '2d' or '3d'
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
    
    # Set y-axis bounds
    # For M_S, M_C, M_S_ratio, M_S_log, and M_X: use consistent bounds for easy comparison
    if metric_key in ['M_S', 'M_C', 'M_S_ratio', 'M_S_log', 'M_X']:
        # Use consistent bounds across these metrics
        ax.set_ylim(-0.7, 0.5)
        y_ticks = np.arange(-0.7, 0.51, 0.1)
        ax.set_yticks(y_ticks)
        # Add horizontal line at r=0
        ax.axhline(y=0, color='black', linestyle='--', linewidth=2, label='r = 0')
    else:
        # Auto-scale for other metrics
        y_min, y_max = np.min(r_vals), np.max(r_vals)
        y_range = y_max - y_min
        if y_range < 0.1:  # minimum range
            y_center = (y_min + y_max) / 2
            y_min, y_max = y_center - 0.05, y_center + 0.05
            y_range = 0.1
        
        # Add padding
        y_pad = y_range * 0.15
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        
        # Set nice tick intervals
        y_span = (y_max + y_pad) - (y_min - y_pad)
        tick_interval = 0.1 if y_span <= 0.8 else 0.2
        y_ticks = np.arange(
            np.floor((y_min - y_pad) / tick_interval) * tick_interval,
            np.ceil((y_max + y_pad) / tick_interval) * tick_interval + tick_interval/2,
            tick_interval
        )
        ax.set_yticks(y_ticks)
        
        # Add horizontal line at r=0 if it's in range
        y_min_plot, y_max_plot = ax.get_ylim()
        if y_min_plot <= 0 <= y_max_plot:
            ax.axhline(y=0, color='black', linestyle='--', linewidth=2, label='r = 0')
    
    # Labels and title based on metric
    lhi_label = r"$\mathregular{LHI_{2D}}$" if lhi_type == '2d' else r"$\mathregular{LHI_{3D}}$"
    
    if metric_key == 'M_S':
        title = f'Metric S ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_C':
        title = f'Metric C ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_S_norm':
        title = f'Metric S_norm ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'best'
    elif metric_key == 'M_S_ratio':
        # Remove "_ratio" from title but keep in filename
        title = f'Metric S ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_S_log':
        # Rename from S_log to just S
        title = f'Metric S ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_X':
        title = f'Metric X ({lhi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    else:
        title = f'{metric_label} vs. {lhi_label}: r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    
    # Font sizes: axis labels 20pt bold, title 22pt bold, extra spacing with labelpad
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_ylabel('r-value', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Legend: explicitly not bold, regular font size
    legend = ax.legend(loc=legend_loc, fontsize=12, frameon=True)
    for text in legend.get_texts():
        text.set_fontweight('normal')
    
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

# ========= Plot LHI vs depth =========
def plot_lhi_vs_depth(lhi_data, site_depths, out_file, lhi_type='2d'):
    """
    Create publication-quality LHI vs depth plot.
    lhi_data: {site_num: {roi: lhi_value}}
    lhi_type: '2d' or '3d'
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site_num, lhi_values in lhi_data.items():
        site_name = site_num_to_name(site_num)
        
        if site_name not in site_depths:
            continue
        
        depth = site_depths[site_name]
        lhi_list = list(lhi_values.values())
        
        if not lhi_list:
            continue
        
        # Site mean and SEM
        mean_lhi = np.mean(lhi_list)
        n = len(lhi_list)
        if n > 1:
            sd_lhi = np.std(lhi_list, ddof=1)
            sem_lhi = sd_lhi / np.sqrt(n)
        else:
            sem_lhi = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_lhi)
        site_sems.append(sem_lhi)
    
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
                ecolor='black', elinewidth=2.0, capsize=5, capthick=2.0,
                alpha=1.0, zorder=3)
    
    # Regression line (half contrast)
    if xx is not None:
        ax.plot(xx, yy, linestyle='--', color='black', linewidth=2, alpha=0.4, zorder=2)
    
    # Stats annotation - BOTTOM LEFT (with left-aligned text)
    if slope is not None and r_val is not None:
        eq_line = f"y = {slope:.6f}x + {intercept:.3f}"
        stats_line = f"r = {r_val:+.3f}, p = {p_val:.3g}, N = {len(site_means_depths)}"
        stats_txt = f"{eq_line}\n{stats_line}"
        ax.text(0.05, 0.05, stats_txt,
                transform=ax.transAxes,
                fontsize=16, fontweight='normal',
                verticalalignment='bottom',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='black', alpha=0.9, linewidth=1.5))
    
    # Labels - axis labels 20pt bold, no title, extra spacing with labelpad
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    if lhi_type == '2d':
        ax.set_ylabel(r'$\mathregular{LHI_{2D}}$ (Local Heterogeneity)', fontsize=20, fontweight='bold', labelpad=15)
    else:
        ax.set_ylabel(r'$\mathregular{LHI_{3D}}$ (Local Heterogeneity)', fontsize=20, fontweight='bold', labelpad=15)
    
    # Auto-scale y-axis with expanded padding to show full error bars
    if len(site_means) > 0:
        # Include error bars in range calculation
        y_min_with_err = np.min(site_means - site_sems)
        y_max_with_err = np.max(site_means + site_sems)
        y_range = y_max_with_err - y_min_with_err
        y_padding = y_range * 0.15  # Increased from 0.1 to 0.15 for better visibility
        
        # Set upper bound to 0.350
        y_upper = 0.350
        ax.set_ylim(y_min_with_err - y_padding, y_upper)
    
    # Add padding to x-axis only
    if len(site_means_depths) > 0:
        x_min, x_max = site_means_depths.min(), site_means_depths.max()
        x_range = x_max - x_min
        x_padding = x_range * 0.05
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
    
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
def process_lhi_analysis(lhi_data, output_dir, lhi_type='2d'):
    """
    Process LHI analysis for all_roi data.
    lhi_type: '2d' or '3d'
    """
    type_label = "LHI2 (2D Local Homogeneity)" if lhi_type == '2d' else "LHI3 (3D Local Homogeneity)"
    
    print(f"\n{'='*60}")
    print(f"Processing all_roi with {type_label}")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Create output directories
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',  # Added M_X
              pearson_site_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs LHI plots
    print(f"\n[INFO] Generating metric vs LHI plots for each site...")
    for metric_file in sorted(METRIC_DATA_ALL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        # Extract site number
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in lhi_data:
            print(f"[SKIP] {site_name}: no LHI data")
            continue
        
        lhi_values = lhi_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_S', 'M_S',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_lhi.png',
            lhi_type=lhi_type
        )
        
        # M_S_norm (percentage of baseline)
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_S_norm', 'M_S_norm',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_lhi.png',
            lhi_type=lhi_type
        )
        
        # M_S_ratio
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_S_ratio', 'M_S_ratio',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_lhi.png',
            lhi_type=lhi_type
        )
        
        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_S_log', 'M_S_log',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_lhi.png',
            lhi_type=lhi_type
        )
        
        # M_C
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_C', 'M_C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_lhi.png',
            lhi_type=lhi_type
        )
        
        # M_X (peak suppression)
        plot_metric_vs_lhi_site(
            site_name, lhi_values, metrics, 'M_X', 'M_X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_lhi.png',
            lhi_type=lhi_type
        )
    
    # Calculate R-values and plot vs depth
    print(f"\n[INFO] Calculating R-values and plotting vs depth...")
    r_values = calculate_r_values_per_site(lhi_data, METRIC_DATA_ALL)
    
    plot_r_vs_depth(r_values['M_S'], site_depths, 'M_S', 'M_S',
                   pearson_site_dir / 'm_s.png', lhi_type=lhi_type)
    plot_r_vs_depth(r_values['M_S_norm'], site_depths, 'M_S_norm', 'M_S_norm',
                   pearson_site_dir / 'm_s_norm.png', lhi_type=lhi_type)
    plot_r_vs_depth(r_values['M_S_ratio'], site_depths, 'M_S_ratio', 'M_S_ratio',
                   pearson_site_dir / 'm_s_r.png', lhi_type=lhi_type)
    plot_r_vs_depth(r_values['M_S_log'], site_depths, 'M_S_log', 'M_S_log',
                   pearson_site_dir / 'm_s_l.png', lhi_type=lhi_type)
    plot_r_vs_depth(r_values['M_C'], site_depths, 'M_C', 'M_C',
                   pearson_site_dir / 'm_c.png', lhi_type=lhi_type)
    plot_r_vs_depth(r_values['M_X'], site_depths, 'M_X', 'M_X',
                   pearson_site_dir / 'm_x.png', lhi_type=lhi_type)
    
    # Plot LHI vs depth
    print(f"\n[INFO] Plotting LHI vs depth...")
    plot_lhi_vs_depth(lhi_data, site_depths,
                     depth_site_dir / 'lhi_vs_depth.png',
                     lhi_type=lhi_type)

# ========= Main =========
def main():
    print("Starting LHI analysis...")
    
    # Load LHI data (both methods)
    print(f"\n[INFO] Loading LHI data from {LHI_FILE}")
    lhi2_data, lhi3_data = load_lhi_data(LHI_FILE)
    print(f"[INFO] Loaded data for {len(lhi2_data)} sites")
    
    # Process LHI2 (2D) analysis
    print("\n" + "="*60)
    print("LHI ANALYSIS - 2D LOCAL HOMOGENEITY (LHI2)")
    print("="*60)
    
    process_lhi_analysis(
        lhi2_data,
        OUTPUT_BASE / '2d',
        lhi_type='2d'
    )
    
    # Process LHI3 (3D) analysis
    print("\n" + "="*60)
    print("LHI ANALYSIS - 3D LOCAL HOMOGENEITY (LHI3)")
    print("="*60)
    
    process_lhi_analysis(
        lhi3_data,
        OUTPUT_BASE / '3d',
        lhi_type='3d'
    )
    
    print(f"\n{'='*60}")
    print("LHI analysis complete!")
    print(f"2D (LHI2) output: {OUTPUT_BASE / '2d'}")
    print(f"3D (LHI3) output: {OUTPUT_BASE / '3d'}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()