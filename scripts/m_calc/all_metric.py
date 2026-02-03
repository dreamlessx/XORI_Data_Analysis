# -*- coding: utf-8 -*-
"""
metric_calc.py
Calculates metrics from tuning curve data and outputs tables only (no plots).
Output columns: M_S | M_C | SNR_g | SNR_p | M_S_norm | M_S_ratio | M_X
"""
from __future__ import unicode_literals
import os
import numpy as np
from scipy.stats import pearsonr

def xplot_read_xori(infile):
    """
    Reads an xplot file containing two tuning curves (grat and plaid) and a baseline.
    Expected file structure:
      - Lines 2-13: first tuning curve ("grat") with 3 values per line.
      - Lines 16-27: second tuning curve ("plaid") with 3 values per line.
      - Line 30: baseline info (the second value on this line is used as the baseline).
    """
    with open(infile, 'r') as fin:
        lines = fin.readlines()
    
    n = 12  # Number of points per tuning curve
    grat = np.empty((3, n), dtype='float32')
    plaid = np.empty((3, n), dtype='float32')
    
    # Read the "grat" tuning curve
    for i in range(n):
        a = lines[2 + i].split()
        grat[0][i] = float(a[0])
        grat[1][i] = float(a[1])
        grat[2][i] = float(a[2])
    
    # Read the "plaid" tuning curve
    for i in range(n):
        a = lines[16 + i].split()
        plaid[0][i] = float(a[0])
        plaid[1][i] = float(a[1])
        plaid[2][i] = float(a[2])
    
    # Extract the baseline from line 30 (the second value)
    base = float(lines[30].split()[1])
    
    return grat, plaid, base

def metric_suppression(plaid, grat, baseline):
    """
    Computes RAW M_S (Metric Suppression) - raw difference.
    Previously called: additive_curve
    Returns: dmean (raw difference)
    """
    tc = grat[1] - baseline
    shift_steps = 3  # for a 90° shift if each step is 30°
    n = len(tc)
    tc_copy = np.copy(tc)
    for i in range(n):
        tc[(i + shift_steps) % n] += tc_copy[i]
    diff = (plaid[1] - baseline) - tc
    dmean = np.mean(diff)
    return dmean

def metric_suppression_normalized(plaid, grat, baseline):
    """
    Computes NORMALIZED M_S as percentage of baseline.
    Previously called: additive_curve_normalized
    Returns: (dmean / baseline) * 100
    """
    dmean = metric_suppression(plaid, grat, baseline)
    return (dmean / baseline) * 100.0

def metric_suppression_ratio(plaid, grat, baseline):
    """
    Computes M_S_ratio.
    Previously called: additive_curve_ratio
    Returns: ratio of plaid mean to shifted grat mean
    """
    tc = grat[1] - baseline
    shift_steps = 3  # for a 90° shift if each step is 30°
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
    """
    Computes M_C (Metric Correlation) - Pearson r between shifted grat and plaid.
    Previously called: curve_association
    """
    tc = grat[1] - baseline
    shift_steps = 3  # for a 90° shift if each step is 30°
    n = len(tc)
    tc_copy = np.copy(tc)
    for i in range(n):
        tc[(i - shift_steps) % n] += tc_copy[i]
    r, _ = pearsonr(tc, plaid[1] - baseline)
    return r

def signal_to_noise(tuning_curve):
    """
    Computes signal-to-noise ratio for a tuning curve.
    """
    eps = 1e-9
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

def main():
    # Root directory containing all site subfolders.
    root_directory = './raw_data/tc_data'
    
    # Output folder for text files only.
    output_folder = os.path.join('metric_data', 'all_roi')
    os.makedirs(output_folder, exist_ok=True)
    
    # Gather all site subfolders (e.g., site002, site003, etc.) in alphabetical order.
    site_folders = sorted(os.listdir(root_directory))
    
    for site_name in site_folders:
        site_path = os.path.join(root_directory, site_name)
        if not os.path.isdir(site_path):
            continue  # Skip non-directory items
        
        # Get ROI files in alphabetical order.
        roi_files = sorted(os.listdir(site_path))
        
        roi_metric_values = []  # Will store tuples of (M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X)
        
        for idx, filename in enumerate(roi_files):
            file_path = os.path.join(site_path, filename)
            grat, plaid, base = xplot_read_xori(file_path)
            
            M_S = metric_suppression(plaid, grat, base)              # Raw metric suppression
            M_C = metric_correlation(plaid, grat, base)              # Metric correlation
            SNR_g = signal_to_noise(grat)                            # SNR for grat
            SNR_p = signal_to_noise(plaid)                           # SNR for plaid
            M_S_norm = metric_suppression_normalized(plaid, grat, base)  # Normalized M_S
            M_S_ratio = metric_suppression_ratio(plaid, grat, base)  # M_S ratio
            M_X = metric_x(plaid, grat, base)                        # Peak suppression metric
            
            roi_metric_values.append((M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X))
            
            print(f"Processed {site_name}/{filename}: "
                  f"M_S = {M_S:.6f}, "
                  f"M_C = {M_C:.6f}, "
                  f"SNR_g = {SNR_g:.4f}, "
                  f"SNR_p = {SNR_p:.4f}, "
                  f"M_S_norm = {M_S_norm:.3f}%, "
                  f"M_S_ratio = {M_S_ratio:.6f}, "
                  f"M_X = {M_X:.6f}")
        
        # Write the results to a text file in a cleaner table format.
        output_textfile = os.path.join(output_folder, f"metrics_{site_name}.txt")
        with open(output_textfile, 'w') as fout:
            # Define column widths for neat alignment (wider spacing)
            roi_col_w = 8
            ms_col_w = 12
            mc_col_w = 12
            snrg_col_w = 12
            snrp_col_w = 12
            ms_norm_col_w = 12
            ms_ratio_col_w = 12
            mx_col_w = 12
            
            # Create format strings
            header_fmt = (f"{{:<{roi_col_w}}}"       # ROI left-aligned
                         f"{{:>{ms_col_w}}}"         # M_S right-aligned
                         f"{{:>{mc_col_w}}}"         # M_C right-aligned
                         f"{{:>{snrg_col_w}}}"       # SNR_g right-aligned
                         f"{{:>{snrp_col_w}}}"       # SNR_p right-aligned
                         f"{{:>{ms_norm_col_w}}}"    # M_S_norm right-aligned
                         f"{{:>{ms_ratio_col_w}}}"   # M_S_ratio right-aligned
                         f"{{:>{mx_col_w}}}\n")      # M_X right-aligned
            
            # Format with 5 significant figures using g format
            data_fmt = (f"{{:<{roi_col_w}}}"
                       f"{{:>{ms_col_w}.5g}}"
                       f"{{:>{mc_col_w}.5g}}"
                       f"{{:>{snrg_col_w}.5g}}"
                       f"{{:>{snrp_col_w}.5g}}"
                       f"{{:>{ms_norm_col_w}.5g}}"
                       f"{{:>{ms_ratio_col_w}.5g}}"
                       f"{{:>{mx_col_w}.5g}}\n")
            
            # Write header line
            fout.write(header_fmt.format("ROI", "M_S", "M_C", "SNR_g", "SNR_p", "M_S_norm", "M_S_ratio", "M_X"))
            
            # Write separator
            sep_len = roi_col_w + ms_col_w + mc_col_w + snrg_col_w + snrp_col_w + ms_norm_col_w + ms_ratio_col_w + mx_col_w
            fout.write("-" * sep_len + "\n")
            
            # Write data rows
            for i, (M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X) in enumerate(roi_metric_values):
                fout.write(data_fmt.format(i, M_S, M_C, SNR_g, SNR_p, M_S_norm, M_S_ratio, M_X))
        
        print(f"Metrics for {site_name} written to {output_textfile}\n")

if __name__ == '__main__':
    main()