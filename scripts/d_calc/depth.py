#!/usr/bin/env python3
"""
depth_data_analysis.py

Generates comprehensive depth analysis with ROI maps and depth vs metric plots.

Structure:
depth_data/
├── all_roi/
├── cull_roi/
│   ├── per_cull/ (top_70, top_80, top_90)
│   └── thr_cull/ (above_0_5, above_1_0, above_1_5)
├── null_roi/
└── r_cull_roi/
    ├── per_r_cull/ (bottom_10, bottom_20, bottom_30)
    └── thr_r_cull/ (below_0_5, below_1_0, below_1_5)

Each endpoint contains:
├── roi_map/ (m_s, m_c, m_s_norm, m_s_r, m_x subfolders with per-site maps)
└── depth_metric/ (5 depth vs metric plots)
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, pearsonr
from pathlib import Path

# ========= Configuration =========
SITE_DEPTH_FILE = Path("./raw_data/bm_data/site_depth.txt")
METRICS_DIR = "metric_data/all_roi"  # Changed from data_site_roi
STAT_DIR = "stat"
NULL_ROI_METRICS_DIR = "metric_data/null_roi"
OUTPUT_BASE = Path("depth_data")

# Metric column indices and names
METRICS = {
    'm_s': {'col': 1, 'label': r'Metric $\mathregular{S_{\Delta}}$ (a.u.)', 'short': 'M_S'},
    'm_c': {'col': 2, 'label': 'Metric C (r-value)', 'short': 'M_C'},
    'm_s_norm': {'col': 5, 'label': r'Metric $\mathregular{S_{\Delta}}$ (% baseline)', 'short': 'M_S_norm'},
    'm_s_r': {'col': 6, 'label': 'Metric S (ratio)', 'short': 'M_S_ratio'},
    'm_x': {'col': 7, 'label': 'Metric X (peak suppression)', 'short': 'M_X'}
}

# ========= Helper Functions =========

def read_site_depths(depth_file):
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
    
    # Skip header and separator (first 2 lines)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            site = parts[0].strip()
            depth = float(parts[1])
            depths[site] = depth
    
    return depths

def _norm(name):
    """Normalize header names for robust matching"""
    return re.sub(r'[^a-z0-9]+', '', name.lower())

def read_per_roi_metric_with_snr(site_name, col_idx, snr_col_idx=3):
    """Read metric values and SNR values for a site, returns values and their original ROI indices"""
    path = os.path.join(METRICS_DIR, f"metrics_{site_name}.txt")
    if not os.path.isfile(path):
        print(f"[DEBUG] File not found: {path}")
        return None, None, None
    
    vals = []
    snrs = []
    roi_indices = []
    
    with open(path, 'r') as f:
        lines = f.readlines()
        print(f"[DEBUG] {site_name}: Read {len(lines)} lines from {path}")
        for i, ln in enumerate(lines[2:]):  # Skip header
            parts = ln.strip().split()
            if len(parts) > max(col_idx, snr_col_idx):
                try:
                    roi_idx = int(parts[0])  # First column is ROI index
                    vals.append(float(parts[col_idx]))
                    snrs.append(float(parts[snr_col_idx]))
                    roi_indices.append(roi_idx)
                except ValueError as e:
                    print(f"[DEBUG] Line {i+2}: ValueError - {e}, parts: {parts}")
                    pass
    
    if not vals:
        print(f"[DEBUG] {site_name}: No values parsed (checked column {col_idx})")
        return None, None, None
    
    print(f"[DEBUG] {site_name}: Successfully read {len(vals)} ROIs")
    return np.array(vals, dtype=float), np.array(snrs, dtype=float), np.array(roi_indices, dtype=int)

def get_null_roi_set():
    """Get set of (site, roi) tuples from null_roi dataset"""
    valid_rois = set()
    
    if not Path(NULL_ROI_METRICS_DIR).exists():
        return valid_rois
    
    for metric_file in Path(NULL_ROI_METRICS_DIR).glob('metrics_*.txt'):
        site_name = metric_file.stem.replace('metrics_', '')
        
        with open(metric_file, 'r') as f:
            for line in f.readlines()[2:]:
                parts = line.strip().split()
                if len(parts) >= 1:
                    try:
                        roi = int(parts[0])
                        valid_rois.add((site_name, roi))
                    except ValueError:
                        continue
    
    return valid_rois

def filter_rois_null(roi_metrics, roi_snrs, roi_indices, site_name, null_roi_set):
    """Filter to only null ROIs"""
    if roi_metrics is None or roi_snrs is None or roi_indices is None:
        return None, None, None
    
    filtered_metrics = []
    filtered_snrs = []
    filtered_indices = []
    
    for i, roi_idx in enumerate(roi_indices):
        if (site_name, roi_idx) in null_roi_set:
            filtered_metrics.append(roi_metrics[i])
            filtered_snrs.append(roi_snrs[i])
            filtered_indices.append(roi_idx)
    
    if not filtered_metrics:
        return None, None, None
    
    return np.array(filtered_metrics), np.array(filtered_snrs), np.array(filtered_indices)

def filter_rois_top_percent(roi_metrics, roi_snrs, roi_indices, percentage):
    """Keep top percentage of ROIs by SNR"""
    if roi_metrics is None or roi_snrs is None or roi_indices is None or len(roi_snrs) == 0:
        return None, None, None
    
    n = len(roi_snrs)
    num_to_take = max(1, int(np.ceil((percentage / 100.0) * n)))
    sorted_idx = np.argsort(-roi_snrs)  # Descending
    selected = sorted_idx[:num_to_take]
    
    return roi_metrics[selected], roi_snrs[selected], roi_indices[selected]

def filter_rois_bottom_percent(roi_metrics, roi_snrs, roi_indices, percentage):
    """Keep bottom percentage of ROIs by SNR"""
    if roi_metrics is None or roi_snrs is None or roi_indices is None or len(roi_snrs) == 0:
        return None, None, None
    
    n = len(roi_snrs)
    num_to_take = max(1, int(np.ceil((percentage / 100.0) * n)))
    sorted_idx = np.argsort(roi_snrs)  # Ascending (lowest SNR first)
    selected = sorted_idx[:num_to_take]
    
    return roi_metrics[selected], roi_snrs[selected], roi_indices[selected]

def filter_rois_above_threshold(roi_metrics, roi_snrs, roi_indices, threshold):
    """Keep only ROIs with SNR >= threshold"""
    if roi_metrics is None or roi_snrs is None or roi_indices is None:
        return None, None, None
    
    idx = np.where(roi_snrs >= threshold)[0]
    if len(idx) == 0:
        return None, None, None
    
    return roi_metrics[idx], roi_snrs[idx], roi_indices[idx]

def filter_rois_below_threshold(roi_metrics, roi_snrs, roi_indices, threshold):
    """Keep only ROIs with SNR < threshold"""
    if roi_metrics is None or roi_snrs is None or roi_indices is None:
        return None, None, None
    
    idx = np.where(roi_snrs < threshold)[0]
    if len(idx) == 0:
        return None, None, None
    
    return roi_metrics[idx], roi_snrs[idx], roi_indices[idx]

# ========= ROI Map Functions =========

def map_roi(xn, yn, bgval, stat, r, g, b):
    """Map RGB values to ROI pixels"""
    sr = np.full((yn, xn), bgval, dtype='float32')
    sg = np.full((yn, xn), bgval, dtype='float32')
    sb = np.full((yn, xn), bgval, dtype='float32')
    nroi = len(stat)
    
    for i in range(nroi):
        x_array = stat[i]['xpix']
        y_array = stat[i]['ypix']
        for j in range(len(x_array)):
            xx = x_array[j]
            yy = y_array[j]
            sr[yy, xx] = r[i]
            sg[yy, xx] = g[i]
            sb[yy, xx] = b[i]
    
    return sr, sg, sb

def plot_roi_metric_map(stat_file, xn, yn, metric, roi_indices, nsd, out_path, title):
    """
    Create ROI map colored by metric values
    roi_indices: list of ROI indices that correspond to the metric values
    """
    if not os.path.isfile(stat_file):
        print(f"[SKIP] stat file not found: {stat_file}")
        return
    
    bgval = 0.0
    stat = np.load(stat_file, allow_pickle=True)
    nroi_total = len(stat)
    
    # Create a mapping from original ROI index to metric value
    roi_to_metric = {}
    for idx, roi_idx in enumerate(roi_indices):
        if roi_idx < len(metric):
            roi_to_metric[roi_idx] = metric[idx]
    
    # Calculate std from available metrics
    sd = metric.std() if len(metric) > 1 else 1e-9
    
    r = np.full(nroi_total, 0.1, dtype='float32')
    g = np.full(nroi_total, 0.1, dtype='float32')
    b = np.full(nroi_total, 0.1, dtype='float32')
    
    for i in range(nroi_total):
        if i in roi_to_metric:
            v = roi_to_metric[i] / (nsd * sd)
            redval = bluval = 0.0
            if v > 0.0:
                redval = min(v, 1.0)
            elif v < 0.0:
                bluval = min(-v, 1.0)
            r[i] = redval + 0.1
            g[i] = 0.1
            b[i] = bluval + 0.1
            if r[i] > 1.0:
                r[i] = 1.0
            if b[i] > 1.0:
                b[i] = 1.0
    
    im = np.empty((3, xn, yn), dtype='float32')
    im[0], im[1], im[2] = map_roi(xn, yn, bgval, stat, r, g, b)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(np.transpose(im, (2, 1, 0)))
    ax.set_title(title, fontsize=14)
    ax.axis('off')
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[SAVED] {out_path}")

# ========= Depth vs Metric Plot Functions =========

def plot_depth_vs_metric(site_data, metric_key, metric_info, out_path, use_log2=False):
    """
    Create depth vs metric plot
    site_data: list of (site_name, depth, filtered_metrics, filtered_snrs)
    """
    if not site_data:
        print(f"[SKIP] No data for {metric_key}")
        return
    
    # Calculate site means
    depths = []
    means = []
    
    for site_name, depth, roi_metrics, roi_snrs in site_data:
        if roi_metrics is not None and len(roi_metrics) > 0:
            if use_log2 and metric_key == 'm_s_r':
                # Filter positive values for log2
                positive = roi_metrics[roi_metrics > 0]
                if len(positive) > 0:
                    mean_val = np.mean(np.log2(positive))
                    depths.append(depth)
                    means.append(mean_val)
            else:
                depths.append(depth)
                means.append(np.mean(roi_metrics))
    
    if len(depths) == 0:
        print(f"[SKIP] No valid data for {metric_key}")
        return
    
    depths = np.array(depths)
    means = np.array(means)
    n = len(depths)
    
    # Linear regression
    lr = linregress(depths, means)
    slope, intercept = lr.slope, lr.intercept
    r_val, p_val = pearsonr(depths, means)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot site means
    ax.scatter(depths, means,
               s=100, color='darkorange',
               edgecolors='black', linewidths=1.5,
               label=f"Mean {metric_info['short']} at each depth",
               zorder=3)
    
    # Regression line
    x_fit = np.linspace(depths.min(), depths.max(), 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit,
            '--', color='black', linewidth=2.5, alpha=0.4,
            label=f"Fit: y = {slope:.4f}x + {intercept:.2f}",
            zorder=2)
    
    # Reference line at y=0 for certain metrics
    if metric_key in ['m_s', 'm_s_r', 'm_x']:
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.6, zorder=1)
    
    # Stats annotation
    eq_line = f"y = {slope:.4f}x + {intercept:.2f}"
    stats_line = f"r = {r_val:+.3f}, p = {p_val:.3g}, N = {n}"
    stats_txt = f"{eq_line}\n{stats_line}"
    
    # Position stats box based on metric type
    if metric_key == 'm_c':
        # Bottom-right area for Metric C (positioned to fit without overflow)
        ax.text(0.50, 0.05, stats_txt,
                transform=ax.transAxes,
                fontsize=16,
                verticalalignment='bottom',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='black', alpha=0.9, linewidth=1.5))
    else:
        # Bottom-left for other metrics
        ax.text(0.05, 0.05, stats_txt,
                transform=ax.transAxes,
                fontsize=16,
                verticalalignment='bottom',
                horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='black', alpha=0.9, linewidth=1.5))
    
    # Labels
    ax.set_xlabel('Depth (μm)', fontsize=33, fontweight='bold')
    ax.set_ylabel(metric_info['label'], fontsize=33, fontweight='bold')
    
    # Auto-scale y-axis with padding
    y_min, y_max = means.min(), means.max()
    y_range = y_max - y_min
    
    if metric_key == 'm_s':
        min_range = 10
    elif metric_key in ['m_c', 'm_x']:
        min_range = 0.1
    else:
        min_range = 0.2
    
    if y_range < min_range:
        y_center = (y_min + y_max) / 2
        y_min, y_max = y_center - min_range/2, y_center + min_range/2
        y_range = min_range
    
    y_pad = y_range * 0.15
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    
    # Custom tick formatting for m_s_r in log2 space
    if use_log2 and metric_key == 'm_s_r':
        y_ticks_log2 = ax.get_yticks()
        y_tick_labels = []
        for tick in y_ticks_log2:
            ratio_val = 2**tick
            if ratio_val >= 1:
                if ratio_val == int(ratio_val):
                    y_tick_labels.append(f'{int(ratio_val)}')
                else:
                    y_tick_labels.append(f'{ratio_val:.1f}')
            else:
                denominator = int(round(1 / ratio_val))
                if abs(ratio_val * denominator - 1) < 0.01:
                    y_tick_labels.append(f'1/{denominator}')
                else:
                    y_tick_labels.append(f'{ratio_val:.2f}')
        ax.set_yticklabels(y_tick_labels)
    
    # Tick styling
    ax.tick_params(axis='both', which='major', labelsize=16, width=1.5, length=6)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Legend
    legend_loc = 'upper left' if metric_key in ['m_c', 'm_x'] else 'best'
    ax.legend(loc=legend_loc, fontsize=14, frameon=True)
    
    plt.tight_layout()
    
    # Save
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.savefig(out_path.with_suffix('.svg'), format='svg', bbox_inches='tight')
    plt.close()
    
    print(f"[SAVED] {out_path}")

# ========= Main Processing Function =========

def process_dataset(dataset_name, filter_func, site_depths, null_roi_set=None):
    """
    Process one dataset (e.g., all_roi, top_70, etc.)
    
    filter_func: function(roi_metrics, roi_snrs, site_name, null_roi_set) -> (filtered_metrics, filtered_snrs)
    """
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")
    
    base_path = OUTPUT_BASE / dataset_name
    roi_map_base = base_path / "roi_map"
    depth_metric_base = base_path / "depth_metric"
    
    # Gather data for all sites and all metrics
    all_site_data = {}  # {metric_key: [(site, depth, roi_metrics, roi_snrs), ...]}
    
    for metric_key in METRICS.keys():
        all_site_data[metric_key] = []
    
    # Sort sites by depth
    sorted_sites = sorted(site_depths.items(), key=lambda x: x[1])
    
    for order, (site_name, depth) in enumerate(sorted_sites, start=1):
        # Read metrics for each metric type
        for metric_key, metric_info in METRICS.items():
            roi_metrics, roi_snrs, roi_indices = read_per_roi_metric_with_snr(site_name, metric_info['col'])
            
            if roi_metrics is None or roi_snrs is None or roi_indices is None:
                continue
            
            # Apply filter
            filtered_metrics, filtered_snrs, filtered_indices = filter_func(roi_metrics, roi_snrs, roi_indices, site_name, null_roi_set)
            
            if filtered_metrics is None or len(filtered_metrics) == 0:
                continue
            
            # Store for depth plots
            all_site_data[metric_key].append((site_name, depth, filtered_metrics, filtered_snrs))
            
            # Create ROI map
            stat_file = os.path.join(STAT_DIR, f"stat_{site_name}.npy")
            map_dir = roi_map_base / metric_key
            map_path = map_dir / f"{order:02d}_roi_map_{site_name}.png"
            
            title = f"{site_name} - {metric_info['short']} (Depth = {depth} µm, Avg = {filtered_metrics.mean():.3f})"
            plot_roi_metric_map(stat_file, 512, 512, filtered_metrics, filtered_indices, 2.0, map_path, title)
    
    # Create depth vs metric plots
    print(f"\n[INFO] Creating depth vs metric plots for {dataset_name}...")
    for metric_key, metric_info in METRICS.items():
        out_path = depth_metric_base / f"{metric_key}_vs_depth"
        use_log2 = (metric_key == 'm_s_r')
        plot_depth_vs_metric(all_site_data[metric_key], metric_key, metric_info, out_path, use_log2)

# ========= Main =========

def main():
    print("Starting depth_data analysis...")
    
    # Load site depths
    site_depths = read_site_depths(SITE_DEPTH_FILE)
    print(f"[INFO] Loaded {len(site_depths)} site depths")
    
    # Load null ROI set
    null_roi_set = get_null_roi_set()
    print(f"[INFO] Loaded {len(null_roi_set)} null ROIs")
    
    # Define all datasets and their filter functions
    # Filter signature: (metrics, snrs, indices, site_name, null_set) -> (filtered_metrics, filtered_snrs, filtered_indices)
    datasets = {
        # All ROIs
        'all_roi': lambda m, s, idx, site, null: (m, s, idx),
        
        # Null ROIs
        'null_roi': lambda m, s, idx, site, null: filter_rois_null(m, s, idx, site, null),
        
        # Cull ROI - Percentage (top performers)
        'cull_roi/per_cull/top_70': lambda m, s, idx, site, null: filter_rois_top_percent(m, s, idx, 70),
        'cull_roi/per_cull/top_80': lambda m, s, idx, site, null: filter_rois_top_percent(m, s, idx, 80),
        'cull_roi/per_cull/top_90': lambda m, s, idx, site, null: filter_rois_top_percent(m, s, idx, 90),
        
        # Cull ROI - Threshold (above SNR threshold)
        'cull_roi/thr_cull/above_0_5': lambda m, s, idx, site, null: filter_rois_above_threshold(m, s, idx, 0.5),
        'cull_roi/thr_cull/above_1_0': lambda m, s, idx, site, null: filter_rois_above_threshold(m, s, idx, 1.0),
        'cull_roi/thr_cull/above_1_5': lambda m, s, idx, site, null: filter_rois_above_threshold(m, s, idx, 1.5),
        
        # Reverse Cull ROI - Percentage (bottom performers)
        'r_cull_roi/per_r_cull/bottom_10': lambda m, s, idx, site, null: filter_rois_bottom_percent(m, s, idx, 10),
        'r_cull_roi/per_r_cull/bottom_20': lambda m, s, idx, site, null: filter_rois_bottom_percent(m, s, idx, 20),
        'r_cull_roi/per_r_cull/bottom_30': lambda m, s, idx, site, null: filter_rois_bottom_percent(m, s, idx, 30),
        
        # Reverse Cull ROI - Threshold (below SNR threshold)
        'r_cull_roi/thr_r_cull/below_0_5': lambda m, s, idx, site, null: filter_rois_below_threshold(m, s, idx, 0.5),
        'r_cull_roi/thr_r_cull/below_1_0': lambda m, s, idx, site, null: filter_rois_below_threshold(m, s, idx, 1.0),
        'r_cull_roi/thr_r_cull/below_1_5': lambda m, s, idx, site, null: filter_rois_below_threshold(m, s, idx, 1.5),
    }
    
    # Process each dataset
    for dataset_name, filter_func in datasets.items():
        process_dataset(dataset_name, filter_func, site_depths, null_roi_set)
    
    print(f"\n{'='*60}")
    print("Depth data analysis complete!")
    print(f"Output directory: {OUTPUT_BASE}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()