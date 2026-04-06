#!/usr/bin/env python3
"""
spatial_analysis.py
Generates spatial frequency analysis plots for all_roi metrics.

Structure:
data_spatial/
├── all_roi/
│   ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_C, M_X vs SF per site)
│   ├── pearson_site/ (R-values vs depth)
│   └── depth_site/ (SF vs depth)
└── all_roi_log10/
    ├── metric_site/ (metrics vs log10 SF per site)
    ├── pearson_site/ (R-values vs depth for log10 SF correlations)
    └── depth_site/ (log10 SF vs depth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
SF_FILE = Path('./raw_data/bm_data/roi_sf.txt')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_ALL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_spatial')

# ========= Load spatial frequency data from roi_sf.txt =========
def load_sf_data(sf_file):
    """
    Read roi_sf.txt and return dict of {site: {roi: sf_value}}
    Format: site | roi | direction | sf_value
    Includes all non-null SF values
    """
    with open(sf_file, 'r') as f:
        lines = f.readlines()
    
    sf_data = {}
    
    # Skip header and separator (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        
        try:
            site = int(float(parts[0]))
            roi = int(float(parts[1]))
            sf = float(parts[3])
            
            if site not in sf_data:
                sf_data[site] = {}
            
            sf_data[site][roi] = sf
            
        except (ValueError, IndexError):
            continue
    
    return sf_data

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

# ========= Plot metric vs SF for one site =========
def plot_metric_vs_sf_site(site_name, sf_values, metrics, metric_key,
                           metric_label, out_file, use_log=False):
    """
    Create scatter plot of metric vs spatial frequency for one site.
    sf_values: dict of {roi: sf_value}
    metrics: dict of {roi: {metric_key: value}}
    use_log: if True, use log-scaled x-axis with linear SF tick labels
    """
    from matplotlib.ticker import FixedLocator, FuncFormatter

    # Match SF to metrics by ROI
    raw_sf = []
    y_vals = []

    for roi_idx in sorted(metrics.keys()):
        if roi_idx in sf_values:
            x_val = sf_values[roi_idx]
            y_val = metrics[roi_idx][metric_key]

            if metric_key == 'M_S_log' and y_val is None:
                continue

            if x_val > 0:
                raw_sf.append(x_val)
                y_vals.append(y_val)

    if len(raw_sf) < 2:
        print(f"[SKIP] {site_name}: not enough data points for {metric_key}")
        return

    raw_sf = np.array(raw_sf)
    y_vals = np.array(y_vals)

    # Stats computed on log10(SF) regardless of display
    log_sf = np.log10(raw_sf)
    r_val, p_val = pearsonr(log_sf, y_vals)
    slope, intercept, *_ = linregress(log_sf, y_vals)
    n = len(raw_sf)

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))

    if use_log:
        # Plot raw SF values on a log-scaled axis
        ax.scatter(raw_sf, y_vals, s=32, alpha=0.7)
        ax.set_xscale('log')

        # Regression line in log space, displayed in data coords
        xx_log = np.linspace(log_sf.min(), log_sf.max(), 200)
        yy = intercept + slope * xx_log
        ax.plot(10**xx_log, yy, '--', color='black', linewidth=2)

        # Tick labels as actual SF values
        sf_unique = np.unique(raw_sf)
        tick_candidates = [0.5, 0.64, 1, 1.5, 2, 3, 4, 5, 6, 7, 8]
        sf_min, sf_max = raw_sf.min() * 0.85, raw_sf.max() * 1.15
        ticks = [t for t in tick_candidates if sf_min <= t <= sf_max]
        if ticks:
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            ax.xaxis.set_major_formatter(FuncFormatter(
                lambda x, _: f'{x:g}'))

        x_label = "Spatial Frequency [cyc/deg]"
        title = f"{site_name}: Spatial Frequency vs {metric_label}"
    else:
        ax.scatter(raw_sf, y_vals, s=32, alpha=0.7)

        # Regression line on raw SF
        slope_raw, intercept_raw, *_ = linregress(raw_sf, y_vals)
        xx = np.linspace(raw_sf.min(), raw_sf.max(), 200)
        yy = intercept_raw + slope_raw * xx
        ax.plot(xx, yy, '--', color='black', linewidth=2)

        x_label = "Spatial Frequency (cyc/deg)"
        title = f"{site_name}: Spatial Frequency vs {metric_label}"

    # Stats annotation (always based on log10 correlation)
    txt = f"y = {slope:+.2f}x {intercept:+.2f}\nPearson r = {r_val:+.2f}, p = {p_val:.2g}, n = {n}"
    ax.text(0.05, 0.95, txt, transform=ax.transAxes, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), fontsize=12)

    if metric_key == 'M_S_norm':
        y_label = "M_S (% of baseline)"
    else:
        y_label = metric_label

    ax.set_title(title, fontsize=14)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True, alpha=0.3)

    # Y-axis padding
    y_range = y_vals.max() - y_vals.min()
    y_padding = y_range * 0.1
    ax.set_ylim(y_vals.min() - y_padding, y_vals.max() + y_padding)

    if not use_log:
        x_range = raw_sf.max() - raw_sf.min()
        x_padding = x_range * 0.05
        ax.set_xlim(raw_sf.min() - x_padding, raw_sf.max() + x_padding)

    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=150)
    plt.close()

# ========= Calculate R-values per site =========
def calculate_r_values_per_site(sf_data, metric_data_dir, use_log=False):
    """
    For each site, calculate Pearson r between SF and each metric.
    Returns: {metric_key: {site_name: (r, p)}}
    sf_data: {site_num: {roi: sf_value}}
    use_log: if True, apply log10 to SF only, not metrics
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
        
        if site_num not in sf_data:
            continue
        
        sf_values = sf_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            
            for roi_idx in sorted(metrics.keys()):
                if roi_idx in sf_values:
                    x_val = sf_values[roi_idx]
                    y_val = metrics[roi_idx][metric_key]
                    
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    
                    # Apply log10 to SF only if requested, skip non-positive SF values
                    if use_log:
                        if x_val > 0:
                            x_vals.append(np.log10(x_val))
                            y_vals.append(y_val)
                    else:
                        x_vals.append(x_val)
                        y_vals.append(y_val)
            
            if len(x_vals) >= 2:
                r, p = pearsonr(x_vals, y_vals)
                results[metric_key][site_name] = (r, p)
    
    return results

# ========= Plot R-values vs depth =========
def plot_r_vs_depth(r_values, site_depths, metric_key, metric_label, out_file, use_log=False):
    """
    Plot Pearson r-values vs depth for one metric.
    r_values: {site_name: (r, p)}
    use_log: if True, indicates this is for log10-transformed data
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
    
    # Set fixed y-axis bounds
    ax.set_ylim(-0.3, 0.6)
    
    # Set tick intervals at 0.1 increments
    y_ticks = np.arange(-0.3, 0.7, 0.1)
    ax.set_yticks(y_ticks)
    
    # Add horizontal line at r=0
    ax.axhline(y=0, color='black', linestyle='--', linewidth=2, label='r = 0')
    
    # Labels and title based on metric
    metric_display = {
        'M_S': 'P_diff',
        'M_C': 'C',
        'M_X': 'X',
        'M_S_norm': 'P_diff_norm',
        'M_S_ratio': 'P',
        'M_S_log': 'P',
    }
    display_name = metric_display.get(metric_key, metric_key)
    title = f'{display_name} (SF): r-value vs. Depth (N={len(depths)})'
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

# ========= Plot SF vs depth =========
def plot_sf_vs_depth(sf_data, site_depths, out_file, use_log=False):
    """
    Create publication-quality spatial frequency vs depth plot.
    sf_data: {site_num: {roi: sf_value}}
    use_log: if True, plot log10(SF) vs depth with exponential y-axis formatting
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site_num, sf_values in sf_data.items():
        site_name = site_num_to_name(site_num)
        
        if site_name not in site_depths:
            continue
        
        depth = site_depths[site_name]
        sf_list = list(sf_values.values())
        
        if not sf_list:
            continue
        
        # Apply log10 if requested
        if use_log:
            sf_list = [np.log10(sf) for sf in sf_list if sf > 0]
            if not sf_list:
                continue
        
        # Site mean and SEM (Standard Error of Mean)
        mean_sf = np.mean(sf_list)
        n = len(sf_list)
        if n > 1:
            sd_sf = np.std(sf_list, ddof=1)
            sem_sf = sd_sf / np.sqrt(n)
        else:
            sem_sf = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_sf)
        site_sems.append(sem_sf)
    
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
    
    if use_log:
        # For log scale, convert log10 values back to actual SF values
        # Data is stored as log10(SF), so convert back with 10^x
        site_means_actual = 10 ** site_means
        site_sems_actual_lower = 10 ** (site_means - site_sems)
        site_sems_actual_upper = 10 ** (site_means + site_sems)
        
        # Calculate asymmetric error bars in linear space
        yerr_lower = site_means_actual - site_sems_actual_lower
        yerr_upper = site_sems_actual_upper - site_means_actual
        yerr = [yerr_lower, yerr_upper]
        
        # Plot with actual SF values
        ax.errorbar(site_means_depths, site_means_actual, yerr=yerr,
                    fmt='o', markersize=8, color='darkorange',
                    markeredgecolor='black', markeredgewidth=1.5,
                    ecolor='black', elinewidth=1.5, capsize=4, capthick=1.5,
                    alpha=1.0, zorder=3)
        
        # Regression line also needs to be converted
        if xx is not None:
            yy_actual = 10 ** yy
            ax.plot(xx, yy_actual, linestyle='--', color='black', linewidth=2, alpha=0.4, zorder=2)
    else:
        # Site means with SEM error bars (orange, prominent)
        ax.errorbar(site_means_depths, site_means, yerr=site_sems,
                    fmt='o', markersize=8, color='darkorange',
                    markeredgecolor='black', markeredgewidth=1.5,
                    ecolor='black', elinewidth=1.5, capsize=4, capthick=1.5,
                    alpha=1.0, zorder=3)
        
        # Regression line (half contrast)
        if xx is not None:
            ax.plot(xx, yy, linestyle='--', color='black', linewidth=2, alpha=0.4, zorder=2)
    
    # Stats annotation - TOP RIGHT
    if slope is not None and r_val is not None:
        eq_line = f"y = {slope:.4f}x + {intercept:.2f}"
        stats_line = f"r = {r_val:+.3f}, p = {p_val:.3g}, N = {len(site_means_depths)}"
        stats_txt = f"{eq_line}\n{stats_line}"
        ax.text(0.95, 0.95, stats_txt,
                transform=ax.transAxes,
                fontsize=16, fontweight='normal',
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='black', alpha=0.9, linewidth=1.5))
    
    # Labels - axis labels 20pt bold, no title, extra spacing with labelpad
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    
    if use_log:
        # For log scale, use plain label and log scale axis
        y_label = 'Spatial Frequency [cyc/deg]'
        ax.set_ylabel(y_label, fontsize=20, fontweight='bold', labelpad=15)
        
        # Use logarithmic scale
        ax.set_yscale('log')
        
        # Set round number ticks - choose based on data range
        # Common spatial frequencies: 2, 3, 4, 5, 6, 7
        from matplotlib.ticker import FixedLocator
        tick_values = [2, 3, 4, 5, 6, 7]
        ax.yaxis.set_major_locator(FixedLocator(tick_values))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x)}'))
        
        # Set y-limits to show the data nicely
        if len(site_means) > 0:
            # Data is already converted to actual SF values
            y_min_actual = 10 ** (site_means.min() - site_sems[site_means.argmin()])
            y_max_actual = 10 ** (site_means.max() + site_sems[site_means.argmax()])
            ax.set_ylim(y_min_actual * 0.9, y_max_actual * 1.1)
    else:
        y_label = 'Spatial Frequency (cyc/deg)'
        ax.set_ylabel(y_label, fontsize=20, fontweight='bold', labelpad=15)
        # Fixed y-axis bounds: 1.75 to 7
        ax.set_ylim(1.75, 7)
    
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
def process_sf_analysis(sf_data, use_log=False):
    """
    Process spatial frequency analysis for all_roi data.
    use_log: if True, apply log10 transformation to SF only, not metrics
    """
    suffix = "_log10" if use_log else ""
    log_desc = "LOG10 " if use_log else ""
    
    print(f"\n{'='*60}")
    print(f"Processing all_roi{suffix}")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Create output directories
    output_dir = OUTPUT_BASE / f'all_roi{suffix}'
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',
              pearson_site_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs SF plots
    print(f"\n[INFO] Generating {log_desc}metric vs SF plots for each site...")
    for metric_file in sorted(METRIC_DATA_ALL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        # Extract site number
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in sf_data:
            print(f"[SKIP] {site_name}: no SF data")
            continue
        
        sf_values = sf_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_S', 'M_S',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_sf.png',
            use_log=use_log
        )
        
        # M_S_norm (percentage of baseline)
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_S_norm', 'M_S_norm',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_sf.png',
            use_log=use_log
        )
        
        # M_S_ratio
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_S_ratio', 'M_S_ratio',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_sf.png',
            use_log=use_log
        )
        
        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_S_log', 'M_S_log',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_sf.png',
            use_log=use_log
        )
        
        # M_C
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_C', 'M_C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_sf.png',
            use_log=use_log
        )
        
        # M_X (peak suppression)
        plot_metric_vs_sf_site(
            site_name, sf_values, metrics, 'M_X', 'M_X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_sf.png',
            use_log=use_log
        )
    
    # Calculate R-values and plot vs depth
    print(f"\n[INFO] Calculating {log_desc}R-values and plotting vs depth...")
    r_values = calculate_r_values_per_site(sf_data, METRIC_DATA_ALL, use_log=use_log)
    
    plot_r_vs_depth(r_values['M_S'], site_depths, 'M_S', 'M_S',
                   pearson_site_dir / 'm_s.png', use_log=use_log)
    plot_r_vs_depth(r_values['M_S_norm'], site_depths, 'M_S_norm', 'M_S_norm',
                   pearson_site_dir / 'm_s_norm.png', use_log=use_log)
    plot_r_vs_depth(r_values['M_S_ratio'], site_depths, 'M_S_ratio', 'M_S_ratio',
                   pearson_site_dir / 'm_s_r.png', use_log=use_log)
    plot_r_vs_depth(r_values['M_S_log'], site_depths, 'M_S_log', 'M_S_log',
                   pearson_site_dir / 'm_s_l.png', use_log=use_log)
    plot_r_vs_depth(r_values['M_C'], site_depths, 'M_C', 'M_C',
                   pearson_site_dir / 'm_c.png', use_log=use_log)
    plot_r_vs_depth(r_values['M_X'], site_depths, 'M_X', 'M_X',
                   pearson_site_dir / 'm_x.png', use_log=use_log)
    
    # Plot SF vs depth
    print(f"\n[INFO] Plotting {log_desc}spatial frequency vs depth...")
    plot_sf_vs_depth(sf_data, site_depths,
                    depth_site_dir / 'sf_vs_depth.png', use_log=use_log)

# ========= Main =========
def main():
    print("Starting spatial frequency analysis...")
    
    # Load SF data
    print(f"\n[INFO] Loading SF data from {SF_FILE}")
    sf_all_data = load_sf_data(SF_FILE)
    print(f"[INFO] Loaded data for {len(sf_all_data)} sites")
    
    # Process standard all_roi analysis
    print("\n" + "="*60)
    print("ALL_ROI SPATIAL FREQUENCY ANALYSIS")
    print("="*60)
    
    process_sf_analysis(sf_all_data, use_log=False)
    
    # Process log10 all_roi analysis
    print("\n" + "="*60)
    print("ALL_ROI LOG10 SPATIAL FREQUENCY ANALYSIS")
    print("="*60)
    
    process_sf_analysis(sf_all_data, use_log=True)
    
    print(f"\n{'='*60}")
    print("Spatial frequency analysis complete!")
    print(f"Standard output: {OUTPUT_BASE / 'all_roi'}")
    print(f"Log10 output: {OUTPUT_BASE / 'all_roi_log10'}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()