#!/usr/bin/env python3
"""
baseline_analysis.py
Generates log10-transformed baseline analysis plots for all_roi metrics.

Structure:
data_baseline/
└── log10/
    └── all_roi/
        ├── metric_site/ (M_S, M_S_ratio, M_S_log, M_C, M_X vs log10(baseline) per site)
        ├── pearson_site/ (R-values vs depth)
        └── depth_site/ (baseline vs depth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ========= Configuration =========
TC_DATA_DIR = Path('./raw_data/tc_data')
SITE_DEPTH_FILE = Path('./raw_data/bm_data/site_depth.txt')
METRIC_DATA_ALL = Path('./metric_data/all_roi')
OUTPUT_BASE = Path('./data_baseline')

N_TUNING = 12

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

# ========= Collect baselines for all ROIs in all sites =========
def collect_all_baselines(tc_data_dir):
    """
    Returns: {site: [baseline_values]}
    """
    site_baselines = {}
    
    for site_folder in sorted(os.listdir(tc_data_dir)):
        site_path = tc_data_dir / site_folder
        if not site_path.is_dir():
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
            site_baselines[site_folder] = baselines
    
    return site_baselines

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
        if len(parts) >= 8:  # Now expect 8 columns instead of 7
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

# ========= Plot metric vs baseline for one site =========
def plot_metric_vs_baseline_site(site, baselines, metrics, metric_key, metric_label, 
                                 out_file):
    """
    Create scatter plot of metric vs log10(baseline) for one site.
    baselines: list of baseline values (one per ROI)
    metrics: dict of {roi: {metric_key: value}}
    """
    # Match baselines to metrics by ROI index
    x_vals = []
    y_vals = []
    
    for roi_idx in sorted(metrics.keys()):
        if roi_idx < len(baselines):
            baseline = baselines[roi_idx]
            y_val = metrics[roi_idx][metric_key]
            
            # Skip if M_S_log is None (negative M_S_ratio)
            if metric_key == 'M_S_log' and y_val is None:
                continue
            
            if baseline > 0:
                x_vals.append(np.log10(baseline))
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
    x_label = "log10(Baseline (a.u.))"
    
    # For M_S, show it as percentage of baseline
    if metric_key == 'M_S_norm':
        y_label = "P_diff (% baseline)"
    else:
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
def calculate_r_values_per_site(site_baselines, metric_data_dir):
    """
    For each site, calculate Pearson r between log10(baseline) and each metric.
    Returns: {metric_key: {site: (r, p)}}
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
        
        if site_name not in site_baselines:
            continue
        
        baselines = site_baselines[site_name]
        metrics = load_metrics_from_file(metric_file)
        
        for metric_key in ['M_S', 'M_S_norm', 'M_S_ratio', 'M_S_log', 'M_C', 'M_X']:
            x_vals = []
            y_vals = []
            
            for roi_idx in sorted(metrics.keys()):
                if roi_idx < len(baselines):
                    baseline = baselines[roi_idx]
                    y_val = metrics[roi_idx][metric_key]
                    
                    # Skip if M_S_log is None (negative M_S_ratio)
                    if metric_key == 'M_S_log' and y_val is None:
                        continue
                    
                    if baseline > 0:
                        x_vals.append(np.log10(baseline))
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
        ylabel = 'P_diff vs. Baseline (r-value)'
        title = f'P_diff (Baseline): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_C':
        ylabel = 'C vs. Baseline (r-value)'
        title = f'C (Baseline): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_S_ratio':
        ylabel = 'P vs. Baseline (r-value)'
        title = f'P (Baseline): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_S_log':
        ylabel = 'log\u2082(P) vs. Baseline (r-value)'
        title = f'log\u2082(P) (Baseline): r-value vs. Depth (N={len(depths)})'
    elif metric_key == 'M_X':
        ylabel = 'X vs. Baseline (r-value)'
        title = f'X (Baseline): r-value vs. Depth (N={len(depths)})'
    else:
        ylabel = f'{metric_label} vs. Baseline (r-value)'
        title = f'{metric_label} (Baseline): r-value vs. Depth (N={len(depths)})'
    
    ax.set_xlabel('Depth (μm)', fontsize=20, fontweight='bold', labelpad=15)
    ax.set_ylabel(ylabel, fontsize=20, fontweight='bold', labelpad=15)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Legend: explicitly not bold
    legend = ax.legend(loc='best', fontsize=12, frameon=True)
    for text in legend.get_texts():
        text.set_fontweight('normal')
    
    # Set y-axis limits to -0.4 to 0.5 for all metrics
    ax.set_ylim(-0.4, 0.5)
    ax.set_yticks(np.arange(-0.4, 0.51, 0.1))
    
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

# ========= Plot baseline vs depth =========
def plot_baseline_vs_depth(site_baselines, site_depths, out_file):
    """
    Create publication-quality log10(baseline) vs depth plot.
    """
    # Collect site means and SEMs
    site_means_depths = []
    site_means = []
    site_sems = []
    
    for site, baselines in site_baselines.items():
        if site not in site_depths:
            continue
        
        depth = site_depths[site]
        
        # Filter valid baselines (positive values for log10)
        valid_baselines = [b for b in baselines if b > 0]
        
        if not valid_baselines:
            continue
        
        # Convert to log10
        log_baselines = [np.log10(b) for b in valid_baselines]
        
        # Site mean and SEM
        mean_baseline = np.mean(log_baselines)
        n = len(log_baselines)
        if n > 1:
            sd_baseline = np.std(log_baselines, ddof=1)
            sem_baseline = sd_baseline / np.sqrt(n)
        else:
            sem_baseline = 0
        
        site_means_depths.append(depth)
        site_means.append(mean_baseline)
        site_sems.append(sem_baseline)
    
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
    ax.set_ylabel(r'$\mathregular{F_{0}}$ Baseline (a.u.)', fontsize=20, fontweight='bold', labelpad=15)
    
    # Set y-axis limits from 100 to 1000 (in log10 space)
    ax.set_ylim(np.log10(100), np.log10(1000))
    
    # Set y-axis ticks at round numbers (in log10 space)
    round_values = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    tick_positions = [np.log10(v) for v in round_values]
    tick_labels = [str(v) for v in round_values]
    
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
    
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
def process_baseline_analysis():
    """
    Process log10 baseline analysis for all_roi data.
    """
    print(f"\n{'='*60}")
    print(f"Processing all_roi with log10 baseline")
    print(f"{'='*60}")
    
    # Load site depths
    site_depths = load_site_depths(SITE_DEPTH_FILE)
    
    # Load baselines
    site_baselines = collect_all_baselines(TC_DATA_DIR)
    print(f"[INFO] Using all ROI baselines")
    
    # Create output directories
    output_dir = OUTPUT_BASE / 'log10' / 'all_roi'
    metric_site_dir = output_dir / 'metric_site'
    pearson_site_dir = output_dir / 'pearson_site'
    depth_site_dir = output_dir / 'depth_site'
    
    for d in [metric_site_dir / 'm_s', metric_site_dir / 'm_s_norm',
              metric_site_dir / 'm_s_r', metric_site_dir / 'm_s_l',
              metric_site_dir / 'm_c', metric_site_dir / 'm_x',  # Added M_X
              pearson_site_dir, depth_site_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Process each site - metric vs baseline plots
    print(f"\n[INFO] Generating metric vs baseline plots for each site...")
    for metric_file in sorted(METRIC_DATA_ALL.glob('metrics_*.txt')):
        site_name = metric_file.stem.replace('metrics_', '')
        
        if site_name not in site_baselines:
            print(f"[SKIP] {site_name}: no baseline data")
            continue
        
        baselines = site_baselines[site_name]
        metrics = load_metrics_from_file(metric_file)
        
        # M_S (raw)
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_S', 'P_diff',
            metric_site_dir / 'm_s' / f'{site_name}_m_s_vs_baseline.png'
        )

        # M_S_norm (percentage of baseline)
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_S_norm', 'P_diff (% baseline)',
            metric_site_dir / 'm_s_norm' / f'{site_name}_m_s_norm_vs_baseline.png'
        )

        # M_S_ratio
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_S_ratio', 'P',
            metric_site_dir / 'm_s_r' / f'{site_name}_m_s_ratio_vs_baseline.png'
        )

        # M_S_log (log2 of positive M_S_ratio values)
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_S_log', 'log\u2082(P)',
            metric_site_dir / 'm_s_l' / f'{site_name}_m_s_log_vs_baseline.png'
        )

        # M_C
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_C', 'C',
            metric_site_dir / 'm_c' / f'{site_name}_m_c_vs_baseline.png'
        )

        # M_X (peak suppression)
        plot_metric_vs_baseline_site(
            site_name, baselines, metrics, 'M_X', 'X',
            metric_site_dir / 'm_x' / f'{site_name}_m_x_vs_baseline.png'
        )
    
    # Calculate R-values and plot vs depth
    print(f"\n[INFO] Calculating R-values and plotting vs depth...")
    r_values = calculate_r_values_per_site(site_baselines, METRIC_DATA_ALL)
    
    plot_r_vs_depth(r_values['M_S'], site_depths, 'M_S', 'P_diff',
                   pearson_site_dir / 'm_s.png')
    plot_r_vs_depth(r_values['M_S_norm'], site_depths, 'M_S_norm', 'P_diff (% baseline)',
                   pearson_site_dir / 'm_s_norm.png')
    plot_r_vs_depth(r_values['M_S_ratio'], site_depths, 'M_S_ratio', 'P',
                   pearson_site_dir / 'm_s_r.png')
    plot_r_vs_depth(r_values['M_S_log'], site_depths, 'M_S_log', 'log\u2082(P)',
                   pearson_site_dir / 'm_s_l.png')
    plot_r_vs_depth(r_values['M_C'], site_depths, 'M_C', 'C',
                   pearson_site_dir / 'm_c.png')
    plot_r_vs_depth(r_values['M_X'], site_depths, 'M_X', 'X',
                   pearson_site_dir / 'm_x.png')
    
    # Plot baseline vs depth
    print(f"\n[INFO] Plotting baseline vs depth...")
    plot_baseline_vs_depth(site_baselines, site_depths,
                          depth_site_dir / 'baseline_vs_depth.png')

# ========= Main =========
def main():
    print("Starting baseline analysis...")
    
    # Process log10 baseline analysis for all_roi only
    print("\n" + "="*60)
    print("LOG10 BASELINE ANALYSIS - ALL ROI")
    print("="*60)
    
    process_baseline_analysis()
    
    print(f"\n{'='*60}")
    print("Baseline analysis complete!")
    print(f"Output directory: {OUTPUT_BASE / 'log10' / 'all_roi'}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()