#!/usr/bin/env python3
"""
size_analysis.py
Generates ROI size (npix) analysis plots for null_roi metrics.
Structure:
data_size/
└── null_roi/
    ├── metric_site/ (M_S, M_S_norm, M_S_ratio, M_S_log, M_C, M_X vs npix per site)
    ├── pearson_site/ (R-values vs depth)
    ├── pearson_comparison/ (Separate plots: Size vs Bandwidth, SF, Baseline)
    └── depth_site/ (npix vs depth)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
SIZE_FILE = Path('./raw_data/bm_data/roi_stat.txt')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_NULL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_size')

# Paths to other data sources
BANDWIDTH_FILE = Path('./raw_data/bm_data/roi_hw_orth.txt')  # Use orth_norm column
SF_FILE = Path('./raw_data/bm_data/roi_sf.txt')  # Use sf column
TC_DATA_DIR = Path('./raw_data/tc_data')  # For baseline extraction

# ========= Baseline reading function =========
def read_baseline_from_xplot(path):
    """
    Read baseline from xplot file (line 30, second value).
    Fallback: search from bottom for line with >= 2 floats.
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Preferred: line 31 (index 30), second value
    try:
        toks = lines[30].split()
        return float(toks[1])
    except Exception:
        # Fallback: search upward
        for ln in reversed(lines):
            toks = ln.split()
            if len(toks) >= 2:
                try:
                    return float(toks[1])
                except Exception:
                    continue
        raise ValueError(f"Could not parse baseline in {path}")

# ========= Collect baselines for all ROIs in all sites =========
def collect_all_baselines(tc_data_dir):
    """
    Returns: {site_num: [baseline_values]}
    """
    site_baselines = {}
    
    for site_folder in sorted(os.listdir(tc_data_dir)):
        site_path = tc_data_dir / site_folder
        if not site_path.is_dir():
            continue
        
        # Extract site number from folder name (e.g., 'site002' -> 2)
        try:
            site_num = int(site_folder.replace('site', ''))
        except ValueError:
            continue
        
        baselines = []
        for roi_file in sorted(os.listdir(site_path)):
            roi_path = site_path / roi_file
            if not roi_path.is_file():
                continue
            try:
                baseline = read_baseline_from_xplot(roi_path)
                baselines.append(baseline)
            except Exception as e:
                print(f"[WARN] Could not read baseline from {site_folder}/{roi_file}: {e}")
        
        if baselines:
            site_baselines[site_num] = baselines
    
    return site_baselines

# ========= Load ROI size data from roi_stat.txt =========
def load_size_data(size_file):
    """
    Read roi_stat.txt and extract npix (3rd value on each line)
    Format expected: site | roi | npix | [other values...]
    Returns: {site: {roi: npix}}
    """
    with open(size_file, 'r') as f:
        lines = f.readlines()
    
    size_data = {}
    # Skip header if present (first 2 lines if they contain header/separator)
    start_idx = 0
    if len(lines) > 0 and ('site' in lines[0].lower() or '---' in lines[0]):
        start_idx = 2
    
    for line in lines[start_idx:]:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        
        try:
            site = int(float(parts[0]))
            roi = int(float(parts[1]))
            npix = float(parts[2])
            
            # Skip invalid sizes
            if npix <= 0:
                continue
            
            # Store in dictionary
            if site not in size_data:
                size_data[site] = {}
            size_data[site][roi] = npix
        except (ValueError, IndexError):
            continue
    
    return size_data

# ========= Load bandwidth data (orth_norm) =========
def load_bandwidth_data(bandwidth_file):
    """
    Read roi_hw_orth.txt and extract orth_norm values (column 4)
    Returns: {site: {roi: bandwidth}}
    """
    bandwidth_data = {}
    with open(bandwidth_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 5:  # Need at least 5 columns
            try:
                site = int(float(parts[0]))
                roi = int(float(parts[1]))
                orth_norm = float(parts[4])  # Column 4 is orth_norm (hw)
                
                if site not in bandwidth_data:
                    bandwidth_data[site] = {}
                bandwidth_data[site][roi] = orth_norm
            except (ValueError, IndexError):
                continue
    
    return bandwidth_data

# ========= Load spatial frequency data =========
def load_sf_data(sf_file):
    """
    Read roi_sf.txt and extract sf values (column 3)
    Returns: {site: {roi: sf}}
    """
    sf_data = {}
    with open(sf_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            try:
                site = int(float(parts[0]))
                roi = int(float(parts[1]))
                sf = float(parts[3])  # Column 3 is sf
                
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

# ========= Calculate correlations between size and other properties =========
def calculate_size_property_correlations(size_data, bandwidth_data, sf_data, baseline_data):
    """
    For each site, calculate Pearson r between size and each property.
    Returns: {
        'bandwidth': {site_name: (r, p)},
        'sf': {site_name: (r, p)},
        'baseline': {site_name: (r, p)}
    }
    """
    results = {
        'bandwidth': {},
        'sf': {},
        'baseline': {}
    }
    
    # Process each site
    for site_num in size_data.keys():
        site_name = site_num_to_name(site_num)
        size_values = size_data[site_num]
        
        # Size vs Bandwidth
        if site_num in bandwidth_data:
            bw_values = bandwidth_data[site_num]
            size_list = []
            bw_list = []
            for roi in size_values.keys():
                if roi in bw_values:
                    size_list.append(size_values[roi])
                    bw_list.append(bw_values[roi])
            
            if len(size_list) >= 2:
                r, p = pearsonr(size_list, bw_list)
                results['bandwidth'][site_name] = (r, p)
        
        # Size vs SF
        if site_num in sf_data:
            sf_values = sf_data[site_num]
            size_list = []
            sf_list = []
            for roi in size_values.keys():
                if roi in sf_values:
                    size_list.append(size_values[roi])
                    sf_list.append(sf_values[roi])
            
            if len(size_list) >= 2:
                r, p = pearsonr(size_list, sf_list)
                results['sf'][site_name] = (r, p)
        
        # Size vs Baseline
        if site_num in baseline_data:
            baselines = baseline_data[site_num]
            size_list = []
            baseline_list = []
            for roi in size_values.keys():
                if roi < len(baselines):
                    baseline = baselines[roi]
                    if baseline > 0:  # Only use valid baselines
                        size_list.append(size_values[roi])
                        baseline_list.append(baseline)
            
            if len(size_list) >= 2:
                r, p = pearsonr(size_list, baseline_list)
                results['baseline'][site_name] = (r, p)
    
    return results

# ========= Plot single property r-values vs depth =========
def plot_single_property_r_vs_depth(r_vals_dict, site_depths, property_name, 
                                     property_color, marker, out_file):
    """
    Plot r-values from size vs single property.
    """
    if not r_vals_dict:
        print(f"[SKIP] No data for Size vs {property_name}")
        return
    
    # Extract depths and r-values
    depths = []
    r_vals = []
    p_vals = []
    
    for site, (r, p) in r_vals_dict.items():
        if site in site_depths:
            depths.append(site_depths[site])
            r_vals.append(r)
            p_vals.append(p)
    
    if len(depths) == 0:
        print(f"[SKIP] No depth data for Size vs {property_name}")
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
                  s=150, c=property_color, marker=marker,
                  label='p < 0.05', edgecolors='black', linewidths=2, zorder=3)
    
    # Non-significant points (hollow)
    if np.any(~sig_mask):
        ax.scatter(depths[~sig_mask], r_vals[~sig_mask],
                  s=150, facecolors='none', edgecolors=property_color,
                  marker=marker, linewidths=2, label='p ≥ 0.05', zorder=3)
    
    # Set fixed y-axis bounds from -0.4 to 0.7
    ax.set_ylim(-0.4, 0.7)
    
    # Set tick intervals at 0.1 increments
    y_ticks = np.arange(-0.4, 0.71, 0.1)
    ax.set_yticks(y_ticks)
    
    # Add horizontal line at r=0
    ax.axhline(y=0, color='black', linestyle='--', linewidth=2, alpha=0.5, zorder=1)
    
    # Title with parentheses
    title = f'ROI Size ({property_name}): r-value vs. Depth'
    
    # Font sizes: axis labels 20pt bold, title 22pt bold
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_ylabel('r-value', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Legend: explicitly not bold, regular font size
    legend = ax.legend(loc='upper right', fontsize=12, frameon=True)
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
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.savefig(out_file.with_suffix('.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"[SAVED] {out_file}")

# ========= Plot metric vs size for one site =========
def plot_metric_vs_size_site(site_name, size_values, metrics, metric_key, 
                              metric_label, out_file):
    """
    Create scatter plot of metric vs ROI size (npix) for one site.
    size_values: dict of {roi: npix}
    metrics: dict of {roi: {metric_key: value}}
    """
    # Match sizes to metrics by ROI
    x_vals = []
    y_vals = []
    for roi_idx in sorted(metrics.keys()):
        if roi_idx in size_values:
            y_val = metrics[roi_idx][metric_key]
            # Skip if M_S_log is None (negative M_S_ratio)
            if metric_key == 'M_S_log' and y_val is None:
                continue
            x_vals.append(size_values[roi_idx])
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
    x_label = "ROI Size (pixels)"
    if metric_key == 'M_S_norm':
        y_label = "M_S (% of baseline)"
    else:
        y_label = metric_label
    
    ax.set_title(f"{site_name}: {x_label} vs {y_label}", fontsize=14)
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
def calculate_r_values_per_site(size_data, metric_data_dir):
    """
    For each site, calculate Pearson r between ROI size and each metric.
    Returns: {metric_key: {site_name: (r, p)}}
    size_data: {site_num: {roi: npix}}
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
        
        if site_num not in size_data:
            continue
        
        size_values = size_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            for roi_idx in sorted(metrics.keys()):
                if roi_idx in size_values:
                    y_val = metrics[roi_idx][metric_key]
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    x_vals.append(size_values[roi_idx])
                    y_vals.append(y_val)
            
            if len(x_vals) >= 2:
                r, p = pearsonr(x_vals, y_vals)
                results[metric_key][site_name] = (r, p)
    
    return results

# ========= Plot R-values vs depth =========
def plot_r_vs_depth(r_values, site_depths, metric_key, metric_label, out_file):
    """
    Plot Pearson r-values vs depth for one metric.
    r_values: {site_name: (r, p)}
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
    ax.set_ylim(-0.5, 0.4)
    
    # Set tick intervals at 0.1 increments
    y_ticks = np.arange(-0.5, 0.41, 0.1)
    ax.set_yticks(y_ticks)
    
    # Add horizontal line at r=0
    ax.axhline(y=0, color='black', linestyle='--', linewidth=2, label='r = 0')
    
    # Labels and title based on metric
    if metric_key == 'M_S':
        title = f'Metric S (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    elif metric_key == 'M_C':
        title = f'Metric C (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    elif metric_key == 'M_X':
        title = f'Metric X (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    elif metric_key == 'M_S_norm':
        title = f'Metric S_norm (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    elif metric_key == 'M_S_ratio':
        # Remove "_ratio" from title but keep in filename
        title = f'Metric S (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    elif metric_key == 'M_S_log':
        # Rename from S_log to just S
        title = f'Metric S (Size): r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    else:
        title = f'{metric_label} vs. Size: r-value vs. Depth (N={len(depths)})'
        legend_loc = 'upper right'
    
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

# ========= Plot ROI size vs depth =========
def plot_size_vs_depth(size_data, site_depths, out_file):
    """
    Create publication-quality ROI size vs depth plot.
    size_data: {site_num: {roi: npix}}
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site_num, size_values in size_data.items():
        site_name = site_num_to_name(site_num)
        if site_name not in site_depths:
            continue
        
        depth = site_depths[site_name]
        size_list = list(size_values.values())
        
        if not size_list:
            continue
        
        # Site mean and SEM
        mean_size = np.mean(size_list)
        n = len(size_list)
        if n > 1:
            sd_size = np.std(size_list, ddof=1)
            sem_size = sd_size / np.sqrt(n)
        else:
            sem_size = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_size)
        site_sems.append(sem_size)
    
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
    ax.set_ylabel('ROI Size (pixels)', fontsize=20, fontweight='bold', labelpad=15)
    
    # Set y-axis bounds: auto bottom with padding, fixed top at 170
    if len(site_means) > 0:
        y_min = site_means.min()
        y_range = site_means.max() - y_min
        y_padding = y_range * 0.1
        ax.set_ylim(y_min - y_padding, 170)
    
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
def process_size_analysis(size_data, bandwidth_data, sf_data, baseline_data):
    """
    Process ROI size (npix) analysis for null_roi data.
    """
    print(f"\n{'='*60}")
    print(f"Processing null_roi with ROI size (npix)")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Create output directories
    output_dir = OUTPUT_BASE / 'null_roi'
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    pearson_comparison_dir = output_dir / 'pearson_comparison'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',
              pearson_site_dir, pearson_comparison_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs size plots
    print(f"\n[INFO] Generating metric vs size plots for each site...")
    for metric_file in sorted(METRIC_DATA_NULL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        # Extract site number
        try:
            site_num = int(site_name.replace('site', ''))
        except ValueError:
            continue
        
        if site_num not in size_data:
            print(f"[SKIP] {site_name}: no size data")
            continue
        
        size_values = size_data[site_num]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_S', 'M_S',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_size.png'
        )
        
        # M_S_norm (percentage of baseline)
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_S_norm', 'M_S_norm',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_size.png'
        )
        
        # M_S_ratio
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_S_ratio', 'M_S_ratio',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_size.png'
        )
        
        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_S_log', 'M_S_log',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_size.png'
        )
        
        # M_C
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_C', 'M_C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_size.png'
        )
        
        # M_X (peak suppression)
        plot_metric_vs_size_site(
            site_name, size_values, metrics, 'M_X', 'M_X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_size.png'
        )
    
    # Calculate R-values and plot vs depth (for metrics)
    print(f"\n[INFO] Calculating R-values and plotting vs depth (metrics)...")
    r_values = calculate_r_values_per_site(size_data, METRIC_DATA_NULL)
    
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
    
    # Calculate correlations between size and other properties
    print(f"\n[INFO] Calculating correlations between size and other properties...")
    property_correlations = calculate_size_property_correlations(
        size_data, bandwidth_data, sf_data, baseline_data
    )
    
    # Generate separate comparison plots
    print(f"\n[INFO] Generating separate comparison plots...")
    
    # Bandwidth plot (circles)
    plot_single_property_r_vs_depth(
        property_correlations['bandwidth'],
        site_depths,
        'Bandwidth',
        '#1f77b4',
        'o',
        pearson_comparison_dir / 'size_bandwidth.png'
    )
    
    # Spatial Frequency plot (squares)
    plot_single_property_r_vs_depth(
        property_correlations['sf'],
        site_depths,
        'Spatial Frequency',
        '#ff7f0e',
        's',
        pearson_comparison_dir / 'size_sf.png'
    )
    
    # Baseline plot (triangles)
    plot_single_property_r_vs_depth(
        property_correlations['baseline'],
        site_depths,
        'Baseline',
        '#2ca02c',
        '^',
        pearson_comparison_dir / 'size_baseline.png'
    )
    
    # Plot size vs depth
    print(f"\n[INFO] Plotting ROI size vs depth...")
    plot_size_vs_depth(size_data, site_depths,
                      depth_site_dir / 'size_vs_depth.png')

# ========= Main =========
def main():
    print("Starting ROI size analysis...")
    
    # Load size data
    print(f"\n[INFO] Loading size data from {SIZE_FILE}")
    size_data = load_size_data(SIZE_FILE)
    print(f"[INFO] Loaded data for {len(size_data)} sites")
    
    # Load bandwidth data (orth_norm)
    print(f"\n[INFO] Loading bandwidth data from {BANDWIDTH_FILE}")
    bandwidth_data = load_bandwidth_data(BANDWIDTH_FILE)
    print(f"[INFO] Loaded bandwidth data for {len(bandwidth_data)} sites")
    
    # Load SF data
    print(f"\n[INFO] Loading SF data from {SF_FILE}")
    sf_data = load_sf_data(SF_FILE)
    print(f"[INFO] Loaded SF data for {len(sf_data)} sites")
    
    # Load baseline data from tc_data
    print(f"\n[INFO] Loading baseline data from {TC_DATA_DIR}")
    baseline_data = collect_all_baselines(TC_DATA_DIR)
    print(f"[INFO] Loaded baseline data for {len(baseline_data)} sites")
    
    # Process size analysis
    print("\n" + "="*60)
    print("ROI SIZE (NPIX) ANALYSIS")
    print("="*60)
    process_size_analysis(size_data, bandwidth_data, sf_data, baseline_data)
    
    print(f"\n{'='*60}")
    print("ROI size analysis complete!")
    print(f"Output directory: {OUTPUT_BASE / 'null_roi'}")
    print(f"  - pearson_comparison/ contains 3 separate plots:")
    print(f"    * size_bandwidth.png (circles)")
    print(f"    * size_sf.png (squares)")
    print(f"    * size_baseline.png (triangles)")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()