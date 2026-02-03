#!/usr/bin/env python3
"""
metric_calc_reverse_culled.py
Calculates metrics from tuning curve data with reverse SNR-based culling (low SNR).
Creates two types of reverse culled outputs:
1. per_r_cull: Filters by bottom percentage of SNR (bottom_10, bottom_20, bottom_30)
2. thr_r_cull: Filters by SNR threshold (below_0.5, below_1.0, below_1.5)

Output columns: ROI | M_S | M_C | SNR_g | SNR_p | M_S_norm | M_S_ratio | M_X
"""

import os
import math
import numpy as np
from scipy.stats import pearsonr

def xplot_read_xori(infile):
    """
    Reads an xplot file containing two tuning curves (grat and plaid) and a baseline.
    """
    with open(infile, 'r') as fin:
        lines = fin.readlines()
    
    n = 12
    grat = np.empty((3, n), dtype='float32')
    plaid = np.empty((3, n), dtype='float32')
    
    for i in range(n):
        a = lines[2 + i].split()
        grat[0][i] = float(a[0])
        grat[1][i] = float(a[1])
        grat[2][i] = float(a[2])
    
    for i in range(n):
        a = lines[16 + i].split()
        plaid[0][i] = float(a[0])
        plaid[1][i] = float(a[1])
        plaid[2][i] = float(a[2])
    
    base = float(lines[30].split()[1])
    
    return grat, plaid, base

def metric_suppression(plaid, grat, baseline):
    tc = grat[1] - baseline
    shift_steps = 3
    n = len(tc)
    tc_copy = np.copy(tc)
    for i in range(n):
        tc[(i + shift_steps) % n] += tc_copy[i]
    diff = (plaid[1] - baseline) - tc
    dmean = np.mean(diff)
    return dmean

def metric_suppression_normalized(plaid, grat, baseline):
    dmean = metric_suppression(plaid, grat, baseline)
    return (dmean / baseline) * 100.0

def metric_suppression_ratio(plaid, grat, baseline):
    tc = grat[1] - baseline
    shift_steps = 3
    n = len(tc)
    tc_copy = np.copy(tc)
    for i in range(n):
        tc[(i + shift_steps) % n] += tc_copy[i]
    tm = np.mean(tc)
    diff = plaid[1] - baseline
    pm = np.mean(diff)
    ratio = pm / tm
    return ratio

def metric_x(plaid, grat, baseline):
    """
    Computes M_X (Metric X) - peak suppression:
    Find the peak (max) response in grating, then compute:
    (grat_peak - plaid_at_peak) / grat_peak
    
    Returns: suppression ratio at preferred orientation
    """
    # Subtract baseline from both curves
    grat_response = grat[1] - baseline
    plaid_response = plaid[1] - baseline
    
    # Find index of peak grating response
    peak_idx = np.argmax(grat_response)
    
    # Get values at peak orientation
    grat_peak = grat_response[peak_idx]
    plaid_at_peak = plaid_response[peak_idx]
    
    # Calculate suppression ratio
    # Avoid division by zero
    if grat_peak == 0.0:
        return 0.0
    
    suppression = (grat_peak - plaid_at_peak) / grat_peak
    
    return suppression

def metric_correlation(plaid, grat, baseline):
    tc = grat[1] - baseline
    shift_steps = 3
    n = len(tc)
    tc_copy = np.copy(tc)
    for i in range(n):
        tc[(i - shift_steps) % n] += tc_copy[i]
    r, _ = pearsonr(tc, plaid[1] - baseline)
    return r

def signal_to_noise(tuning_curve):
    dmu = tuning_curve[1]
    dsd = tuning_curve[2]
    m = dmu.shape[0]
    fn = 6
    sh2 = np.mean(dsd**2)
    ybar = np.mean(dmu)
    ssq = np.sum((dmu - ybar)**2)
    der2 = (ssq - ((m - 1) * sh2 / fn)) / m
    if sh2 == 0.0:
        return 0.0
    return der2 / sh2

def write_metrics_file(output_file, roi_metric_values):
    """Write metrics to file in standard format"""
    with open(output_file, 'w') as fout:
        roi_col_w = 8
        ms_col_w = 12
        mc_col_w = 12
        snrg_col_w = 12
        snrp_col_w = 12
        ms_norm_col_w = 12
        ms_ratio_col_w = 12
        mx_col_w = 12
        
        header_fmt = (f"{{:<{roi_col_w}}}"
                     f"{{:>{ms_col_w}}}"
                     f"{{:>{mc_col_w}}}"
                     f"{{:>{snrg_col_w}}}"
                     f"{{:>{snrp_col_w}}}"
                     f"{{:>{ms_norm_col_w}}}"
                     f"{{:>{ms_ratio_col_w}}}"
                     f"{{:>{mx_col_w}}}\n")
        
        data_fmt = (f"{{:<{roi_col_w}}}"
                   f"{{:>{ms_col_w}.5g}}"
                   f"{{:>{mc_col_w}.5g}}"
                   f"{{:>{snrg_col_w}.5g}}"
                   f"{{:>{snrp_col_w}.5g}}"
                   f"{{:>{ms_norm_col_w}.5g}}"
                   f"{{:>{ms_ratio_col_w}.5g}}"
                   f"{{:>{mx_col_w}.5g}}\n")
        
        fout.write(header_fmt.format("ROI", "M_S", "M_C", "SNR_g", "SNR_p", "M_S_norm", "M_S_ratio", "M_X"))
        
        sep_len = roi_col_w + ms_col_w + mc_col_w + snrg_col_w + snrp_col_w + ms_norm_col_w + ms_ratio_col_w + mx_col_w
        fout.write("-" * sep_len + "\n")
        
        for roi_idx, M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X in roi_metric_values:
            fout.write(data_fmt.format(roi_idx, M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X))

def main():
    root_directory = './raw_data/tc_data'
    
    # Define reverse culling criteria (low SNR)
    percentages = [10, 20, 30]  # bottom_10, bottom_20, bottom_30
    thresholds = [0.5, 1.0, 1.5]  # below_0.5, below_1.0, below_1.5
    
    # Create output directories
    base_output = 'metric_data/r_cull_roi'
    per_r_cull_dir = os.path.join(base_output, 'per_r_cull')
    thr_r_cull_dir = os.path.join(base_output, 'thr_r_cull')
    
    for perc in percentages:
        os.makedirs(os.path.join(per_r_cull_dir, f'bottom_{perc}'), exist_ok=True)
    
    for thr in thresholds:
        thr_str = str(thr).replace('.', '_')
        os.makedirs(os.path.join(thr_r_cull_dir, f'below_{thr_str}'), exist_ok=True)
    
    # Process each site
    site_folders = sorted(os.listdir(root_directory))
    
    for site_name in site_folders:
        site_path = os.path.join(root_directory, site_name)
        if not os.path.isdir(site_path):
            continue
        
        roi_files = sorted(os.listdir(site_path))
        
        # Collect all ROI data with indices
        all_roi_data = []
        
        for idx, filename in enumerate(roi_files):
            file_path = os.path.join(site_path, filename)
            grat, plaid, base = xplot_read_xori(file_path)
            
            M_S = metric_suppression(plaid, grat, base)
            M_C = metric_correlation(plaid, grat, base)
            SNR_g = signal_to_noise(grat)
            SNR_p = signal_to_noise(plaid)
            M_S_norm = metric_suppression_normalized(plaid, grat, base)
            M_S_ratio = metric_suppression_ratio(plaid, grat, base)
            M_X = metric_x(plaid, grat, base)
            
            all_roi_data.append((idx, M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X))
        
        if len(all_roi_data) == 0:
            continue
        
        print(f"\nProcessing {site_name} ({len(all_roi_data)} ROIs)")
        
        # Process percentage-based reverse culling (BOTTOM percentages)
        for perc in percentages:
            # Calculate how many ROIs to keep
            n_total = len(all_roi_data)
            n_keep = max(1, int(math.ceil((perc / 100.0) * n_total)))
            
            # Sort by SNR_g (ascending) and take bottom N
            sorted_data = sorted(all_roi_data, key=lambda x: x[3])  # ascending
            bottom_rois = sorted_data[:n_keep]
            
            # Re-sort by original ROI index for output
            bottom_rois_sorted = sorted(bottom_rois, key=lambda x: x[0])
            
            # Prepare output data
            output_data = [(roi[0], roi[1], roi[2], roi[3], roi[4], roi[5], roi[6], roi[7]) 
                          for roi in bottom_rois_sorted]
            
            # Write to file
            output_file = os.path.join(per_r_cull_dir, f'bottom_{perc}', f'metrics_{site_name}.txt')
            write_metrics_file(output_file, output_data)
            print(f"  bottom_{perc}: kept {len(bottom_rois)} ROIs")
        
        # Process threshold-based reverse culling (BELOW threshold)
        for thr in thresholds:
            # Filter ROIs by SNR below threshold
            filtered_rois = [roi for roi in all_roi_data if roi[3] < thr]
            
            if len(filtered_rois) == 0:
                print(f"  below_{thr}: no ROIs below threshold (skipped)")
                continue
            
            # Sort by original ROI index for output
            filtered_rois_sorted = sorted(filtered_rois, key=lambda x: x[0])
            
            # Prepare output data
            output_data = [(roi[0], roi[1], roi[2], roi[3], roi[4], roi[5], roi[6], roi[7]) 
                          for roi in filtered_rois_sorted]
            
            # Write to file
            thr_str = str(thr).replace('.', '_')
            output_file = os.path.join(thr_r_cull_dir, f'below_{thr_str}', f'metrics_{site_name}.txt')
            write_metrics_file(output_file, output_data)
            print(f"  below_{thr}: kept {len(filtered_rois)} ROIs")

if __name__ == '__main__':
    main()