#!/usr/bin/env python3
"""
osi_analysis.py
Generates orientation selectivity index analysis plots for all_roi metrics.

Structure:
data_osi/
├── variance/    # Using raw circular variance (OCV - lower = more selective)
│   ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_S_log, M_C, M_X vs OCV per site)
│   ├── pearson_site/ (R-values vs depth, inverted y-axis)
│   └── depth_site/ (OCV vs depth)
└── osi/         # Using Gaku's method (depth of modulation - higher = more selective)
    ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_S_log, M_C, M_X vs OSI per site)
    ├── pearson_site/ (R-values vs depth)
    └── depth_site/ (OSI vs depth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
OSI_FILE = Path('./raw_data/bm_data/roi_osi.txt')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_ALL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_osi')

# ========= Load OSI data from roi_osi.txt =========
def load_osi_data(osi_file):
    """
    Read roi_osi.txt and return dicts for both OSI methods
    Format: site | roi | bsf | bdr | dsi | osi | dcv | ocv
    Returns: 
        osi_gaku: {site: {roi: osi_value}} - Gaku's method (higher = more selective)
        osi_cv: {site: {roi: ocv_value}} - Raw circular variance (lower = more selective)
    """
    with open(osi_file, 'r') as f:
        lines = f.readlines()
    
    osi_gaku = {}
    osi_cv = {}
    
    # Skip header and separator (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) < 8:
            continue
        
        try:
            site = int(float(parts[0]))
            roi = int(float(parts[1]))
            osi_val = float(parts[5])  # Gaku's OSI
            ocv = float(parts[7])      # Orientation circular variance
            
            # Gaku's OSI
            if site not in osi_gaku:
                osi_gaku[site] = {}
            osi_gaku[site][roi] = osi_val
            
            # Circular variance (raw OCV values)
            if site not in osi_cv:
                osi_cv[site] = {}
            osi_cv[site][roi] = ocv
            
        except (ValueError, IndexError):
            continue
    
    return osi_gaku, osi_cv

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
    M_X is peak suppression at preferred orientation
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
                    'M_X': float(parts[7])
                }
            except ValueError:
                continue
    
    return metrics

# ========= Convert site number to site name =========
def site_num_to_name(site_num):
    """Convert site number (e.g., 2) to site name (e.g., 'site002')"""
    return f"site{int(site_num):03d}"

# ========= Plot metric vs OSI for one site =========
def plot_metric_vs_osi_site(site_name, osi_values, metrics, metric_key, 
                            metric_label, out_file, osi_type='gaku'):
    """
    Create scatter plot of metric vs OSI for one site.
    osi_values: dict of {roi: osi_value}
    metrics: dict of {roi: {metric_key: value}}
    osi_type: 'gaku' or 'cv'
    """
    # Match OSI to metrics by ROI
    x_vals = []
    y_vals = []
    
    for roi_idx in sorted(metrics.keys()):
        if roi_idx in osi_values:
            y_val = metrics[roi_idx][metric_key]
            
            # Skip if M_S_log is None (negative M_S_ratio)
            if metric_key == 'M_S_log' and y_val is None:
                continue
            
            x_vals.append(osi_values[roi_idx])
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
    osi_label = "Circular Variance" if osi_type == 'cv' else "OSI(ratio)"
    
    if metric_key == 'M_S_norm':
        y_label = "P_diff (% baseline)"
    else:
        y_label = metric_label

    ax.set_title(f"{site_name}: {osi_label} vs {y_label}", fontsize=14)
    ax.set_xlabel(osi_label, fontsize=12)
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
def calculate_r_values_per_site(osi_data, metric_data_dir):
    """
    For each site, calculate Pearson r between OSI and each metric.
    Returns: {metric_key: {site_name: (r, p)}}
    osi_data: {site_num: {roi: osi_value}}
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
        
        # Extract site number from site_name (e.g., 'site002' -> 2)
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in osi_data:
            continue
        
        osi_values = osi_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            
            for roi_idx in sorted(metrics.keys()):
                if roi_idx in osi_values:
                    y_val = metrics[roi_idx][metric_key]
                    
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    
                    x_vals.append(osi_values[roi_idx])
                    y_vals.append(y_val)
            
            if len(x_vals) >= 2:
                r, p = pearsonr(x_vals, y_vals)
                results[metric_key][site_name] = (r, p)
    
    return results

# ========= Plot R-values vs depth =========
def plot_r_vs_depth(r_values, site_depths, metric_key, metric_label, out_file, osi_type='gaku'):
    """
    Plot Pearson r-values vs depth for one metric.
    r_values: {site_name: (r, p)}
    osi_type: 'gaku' or 'cv'
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
        # Use different bounds for CV vs OSI to emphasize expected correlation direction
        if osi_type == 'cv':
            # For circular variance, emphasize negative correlations
            ax.set_ylim(-0.7, 0.5)
            y_ticks = np.arange(-0.7, 0.51, 0.1)
            ax.set_yticks(y_ticks)
        else:
            # For OSI, emphasize positive correlations
            ax.set_ylim(-0.5, 0.7)
            y_ticks = np.arange(-0.5, 0.71, 0.1)
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
    osi_label = "Circular Variance" if osi_type == 'cv' else "OSI(ratio)"
    
    if metric_key == 'M_S':
        title = f'P_diff ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_C':
        title = f'C ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_X':
        title = f'X ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_S_norm':
        title = f'P_diff (% baseline) ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'best'
    elif metric_key == 'M_S_ratio':
        title = f'P ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    elif metric_key == 'M_S_log':
        title = f'log\u2082(P) ({osi_label}): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper left'
    else:
        title = f'{metric_label} vs. {osi_label}: r-value vs. Depth (N={len(depths)})'
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

# ========= Plot OSI vs depth =========
def plot_osi_vs_depth(osi_data, site_depths, out_file, osi_type='gaku'):
    """
    Create publication-quality OSI vs depth plot.
    osi_data: {site_num: {roi: osi_value}}
    osi_type: 'gaku' or 'cv'
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site_num, osi_values in osi_data.items():
        site_name = site_num_to_name(site_num)
        
        if site_name not in site_depths:
            continue
        
        depth = site_depths[site_name]
        osi_list = list(osi_values.values())
        
        if not osi_list:
            continue
        
        # Site mean and SEM
        mean_osi = np.mean(osi_list)
        n = len(osi_list)
        if n > 1:
            sd_osi = np.std(osi_list, ddof=1)
            sem_osi = sd_osi / np.sqrt(n)
        else:
            sem_osi = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_osi)
        site_sems.append(sem_osi)
    
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
    
    # Stats annotation - position based on osi_type
    if slope is not None and r_val is not None:
        eq_line = f"y = {slope:.6f}x + {intercept:.3f}"
        stats_line = f"r = {r_val:+.3f}, p = {p_val:.3g}, N = {len(site_means_depths)}"
        stats_txt = f"{eq_line}\n{stats_line}"
        
        if osi_type == 'cv':
            # Top left for circular variance
            ax.text(0.05, 0.95, stats_txt,
                    transform=ax.transAxes,
                    fontsize=16, fontweight='normal',
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox=dict(boxstyle='round', facecolor='white',
                             edgecolor='black', alpha=0.9, linewidth=1.5))
        else:
            # Bottom left for Gaku's method
            ax.text(0.05, 0.05, stats_txt,
                    transform=ax.transAxes,
                    fontsize=16, fontweight='normal',
                    verticalalignment='bottom',
                    horizontalalignment='left',
                    bbox=dict(boxstyle='round', facecolor='white',
                             edgecolor='black', alpha=0.9, linewidth=1.5))
    
    # Labels - axis labels 20pt bold, no title, extra spacing with labelpad
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    if osi_type == 'cv':
        ax.set_ylabel('Circular Variance', fontsize=20, fontweight='bold', labelpad=15)
    else:
        ax.set_ylabel('OSI(ratio)', fontsize=20, fontweight='bold', labelpad=15)
    
    # Auto-scale y-axis with expanded padding to show full error bars
    if len(site_means) > 0:
        # Include error bars in range calculation
        y_min_with_err = np.min(site_means - site_sems)
        y_max_with_err = np.max(site_means + site_sems)
        y_range = y_max_with_err - y_min_with_err
        y_padding = y_range * 0.15  # Increased from 0.1 to 0.15 for better visibility
        ax.set_ylim(y_min_with_err - y_padding, y_max_with_err + y_padding)
    
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
def process_osi_analysis(osi_data, output_dir, osi_type='gaku'):
    """
    Process OSI analysis for all_roi data.
    osi_type: 'gaku' or 'cv'
    """
    type_label = "Raw Circular Variance (OCV)" if osi_type == 'cv' else "Gaku's Method"
    
    print(f"\n{'='*60}")
    print(f"Processing all_roi with OSI - {type_label}")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Create output directories
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',
              pearson_site_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs OSI plots
    print(f"\n[INFO] Generating metric vs OSI plots for each site...")
    for metric_file in sorted(METRIC_DATA_ALL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        # Extract site number
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in osi_data:
            print(f"[SKIP] {site_name}: no OSI data")
            continue
        
        osi_values = osi_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_S', 'P_diff',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_osi.png',
            osi_type=osi_type
        )

        # M_S_norm (percentage of baseline)
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_S_norm', 'P_diff (% baseline)',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_osi.png',
            osi_type=osi_type
        )

        # M_S_ratio
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_S_ratio', 'P',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_osi.png',
            osi_type=osi_type
        )

        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_S_log', 'log\u2082(P)',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_osi.png',
            osi_type=osi_type
        )

        # M_C
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_C', 'C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_osi.png',
            osi_type=osi_type
        )

        # M_X (peak suppression)
        plot_metric_vs_osi_site(
            site_name, osi_values, metrics, 'M_X', 'X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_osi.png',
            osi_type=osi_type
        )
    
    # Calculate R-values and plot vs depth
    print(f"\n[INFO] Calculating R-values and plotting vs depth...")
    r_values = calculate_r_values_per_site(osi_data, METRIC_DATA_ALL)
    
    plot_r_vs_depth(r_values['M_S'], site_depths, 'M_S', 'P_diff',
                   pearson_site_dir / 'm_s.png', osi_type=osi_type)
    plot_r_vs_depth(r_values['M_S_norm'], site_depths, 'M_S_norm', 'P_diff (% baseline)',
                   pearson_site_dir / 'm_s_norm.png', osi_type=osi_type)
    plot_r_vs_depth(r_values['M_S_ratio'], site_depths, 'M_S_ratio', 'P',
                   pearson_site_dir / 'm_s_r.png', osi_type=osi_type)
    plot_r_vs_depth(r_values['M_S_log'], site_depths, 'M_S_log', 'log\u2082(P)',
                   pearson_site_dir / 'm_s_l.png', osi_type=osi_type)
    plot_r_vs_depth(r_values['M_C'], site_depths, 'M_C', 'C',
                   pearson_site_dir / 'm_c.png', osi_type=osi_type)
    plot_r_vs_depth(r_values['M_X'], site_depths, 'M_X', 'X',
                   pearson_site_dir / 'm_x.png', osi_type=osi_type)
    
    # Plot OSI vs depth
    print(f"\n[INFO] Plotting OSI vs depth...")
    plot_osi_vs_depth(osi_data, site_depths,
                     depth_site_dir / 'osi_vs_depth.png',
                     osi_type=osi_type)

# ========= Main =========
def main():
    print("Starting OSI analysis...")
    
    # Load OSI data (both methods)
    print(f"\n[INFO] Loading OSI data from {OSI_FILE}")
    osi_gaku_data, osi_cv_data = load_osi_data(OSI_FILE)
    print(f"[INFO] Loaded data for {len(osi_gaku_data)} sites")
    
    # Process Gaku's OSI analysis
    print("\n" + "="*60)
    print("OSI ANALYSIS - GAKU'S METHOD")
    print("="*60)
    
    process_osi_analysis(
        osi_gaku_data,
        OUTPUT_BASE / 'osi',
        osi_type='gaku'
    )
    
    # Process Circular Variance OSI analysis
    print("\n" + "="*60)
    print("OSI ANALYSIS - CIRCULAR VARIANCE METHOD")
    print("="*60)
    
    process_osi_analysis(
        osi_cv_data,
        OUTPUT_BASE / 'variance',
        osi_type='cv'
    )
    
    print(f"\n{'='*60}")
    print("OSI analysis complete!")
    print(f"Gaku's method output: {OUTPUT_BASE / 'osi'}")
    print(f"Circular variance output: {OUTPUT_BASE / 'variance'}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()