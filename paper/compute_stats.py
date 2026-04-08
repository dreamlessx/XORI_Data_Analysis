#!/usr/bin/env python3
"""
Compute all statistics for the XORI manuscript using log2(P) as the primary metric.

Outputs exact numbers for manuscript tables:
  - Partial correlations (Table 2)
  - Mediation analysis (Table 1)
  - Bootstrap CIs
  - Summary stats

P = M_S_ratio (observed/predicted ratio)
All correlations computed on log2(P) at site level.
"""

import numpy as np
from numpy.linalg import lstsq
from scipy.stats import pearsonr, linregress
from pathlib import Path
import json

ROOT = Path('.')
RAW_DATA = ROOT / 'raw_data' / 'bm_data'
METRIC_DATA = ROOT / 'metric_data' / 'all_roi'


def load_site_depths():
    depths = {}
    with open(RAW_DATA / 'site_depth.txt') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                depths[parts[0]] = float(parts[1])
    return depths


def load_metrics(site_name):
    path = METRIC_DATA / f'metrics_{site_name}.txt'
    data = {'ROI': [], 'M_S': [], 'M_C': [], 'SNR_g': [], 'SNR_p': [],
            'M_S_norm': [], 'M_S_ratio': [], 'M_X': []}
    with open(path) as f:
        for line in f.readlines()[2:]:
            p = line.strip().split()
            if len(p) >= 8:
                try:
                    data['ROI'].append(int(p[0]))
                    data['M_S'].append(float(p[1]))
                    data['M_C'].append(float(p[2]))
                    data['SNR_g'].append(float(p[3]))
                    data['SNR_p'].append(float(p[4]))
                    data['M_S_norm'].append(float(p[5]))
                    data['M_S_ratio'].append(float(p[6]))
                    data['M_X'].append(float(p[7]))
                except ValueError:
                    continue
    return {k: np.array(v) for k, v in data.items()}


def load_covariate(filename, col_idx=3):
    data = {}
    with open(RAW_DATA / filename) as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= col_idx + 1:
                try:
                    site = int(float(parts[0]))
                    roi = int(float(parts[1]))
                    val = float(parts[col_idx])
                    if site not in data:
                        data[site] = {}
                    data[site][roi] = val
                except (ValueError, IndexError):
                    continue
    return data


def load_roi_stat():
    data = {}
    with open(RAW_DATA / 'roi_stat.txt') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 14:
                try:
                    site = int(parts[0])
                    roi = int(parts[1])
                    if site not in data:
                        data[site] = {}
                    data[site][roi] = {
                        'npix': int(parts[2]),
                        'radius': float(parts[6]),
                    }
                except (ValueError, IndexError):
                    continue
    return data


def partial_corr(x, y, covariates):
    """Partial correlation between x and y controlling for covariates."""
    cov = np.column_stack([np.ones(len(covariates)), covariates])
    beta_x, _, _, _ = lstsq(cov, x, rcond=None)
    x_resid = x - cov @ beta_x
    beta_y, _, _, _ = lstsq(cov, y, rcond=None)
    y_resid = y - cov @ beta_y
    return pearsonr(x_resid, y_resid)


def build_site_table(sd):
    """Build site-level table with all covariates using log2(P)."""
    sf_data = load_covariate('roi_sf.txt')
    osi_data = load_covariate('roi_osi.txt', col_idx=5)
    lhi_data = load_covariate('roi_lhi.txt')
    hw_data = load_covariate('roi_hw_orth.txt')
    stat_data = load_roi_stat()

    rows = []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        try:
            m = load_metrics(site)
        except FileNotFoundError:
            continue

        site_num = int(site.replace('site', ''))

        # log2(P) for positive ratios
        pos_mask = m['M_S_ratio'] > 0
        if pos_mask.sum() == 0:
            continue

        log2_p = np.log2(m['M_S_ratio'][pos_mask])
        mc_vals = m['M_C'][pos_mask]
        snr_vals = m['SNR_g'][pos_mask]
        rois = m['ROI'][pos_mask]

        # Covariates (matched to positive-ratio ROIs)
        def get_cov(cov_dict, rois_arr):
            vals = []
            for roi in rois_arr:
                if site_num in cov_dict and int(roi) in cov_dict[site_num]:
                    vals.append(cov_dict[site_num][int(roi)])
                else:
                    vals.append(np.nan)
            return np.array(vals)

        sf_vals = get_cov(sf_data, rois)
        osi_vals = get_cov(osi_data, rois)
        lhi_vals = get_cov(lhi_data, rois)
        hw_vals = get_cov(hw_data, rois)

        # Radius
        radius_vals = []
        for roi in rois:
            if site_num in stat_data and int(roi) in stat_data[site_num]:
                radius_vals.append(stat_data[site_num][int(roi)]['radius'])
            else:
                radius_vals.append(np.nan)
        radius_vals = np.array(radius_vals)

        rows.append({
            'site': site,
            'depth': depth,
            'log2_p': np.nanmean(log2_p),
            'mc': np.nanmean(mc_vals),
            'snr': np.nanmean(snr_vals),
            'sf': np.nanmean(sf_vals),
            'osi': np.nanmean(osi_vals),
            'lhi': np.nanmean(lhi_vals),
            'hw': np.nanmean(hw_vals),
            'radius': np.nanmean(radius_vals),
            'n_roi': int(pos_mask.sum()),
        })

    return rows


def main():
    print('=' * 60)
    print('XORI Manuscript Statistics — log2(P)')
    print('=' * 60)

    sd = load_site_depths()
    rows = build_site_table(sd)

    # Extract arrays
    depth = np.array([r['depth'] for r in rows])
    log2_p = np.array([r['log2_p'] for r in rows])
    mc = np.array([r['mc'] for r in rows])
    snr = np.array([r['snr'] for r in rows])
    sf = np.array([r['sf'] for r in rows])
    osi = np.array([r['osi'] for r in rows])
    lhi = np.array([r['lhi'] for r in rows])
    hw = np.array([r['hw'] for r in rows])
    radius = np.array([r['radius'] for r in rows])
    n = len(depth)

    results = {'n_sites': n}

    # ========================================================
    # 1. Zero-order correlations
    # ========================================================
    print(f'\nn = {n} sites\n')
    print('--- Zero-order correlations ---')

    r_p, p_p = pearsonr(depth, log2_p)
    r_c, p_c = pearsonr(depth, mc)
    print(f'log2(P) ~ depth: r = {r_p:+.3f}, p = {p_p:.2e}')
    print(f'C       ~ depth: r = {r_c:+.3f}, p = {p_c:.2e}')

    results['zero_order'] = {
        'P_depth': {'r': round(r_p, 3), 'p': f'{p_p:.2e}'},
        'C_depth': {'r': round(r_c, 3), 'p': f'{p_c:.2e}'},
    }

    # ========================================================
    # 2. Partial correlations (Table 2)
    # ========================================================
    print('\n--- Partial correlations: P-depth ---')
    controls = [
        ('ROI radius', radius.reshape(-1, 1)),
        ('SNR', snr.reshape(-1, 1)),
        ('OSI', osi.reshape(-1, 1)),
        ('SF', sf.reshape(-1, 1)),
        ('Radius + SNR', np.column_stack([radius, snr])),
        ('Radius + SNR + OSI + SF', np.column_stack([radius, snr, osi, sf])),
    ]

    partial_results_p = []
    partial_results_c = []
    for label, cov in controls:
        # Drop NaN rows
        mask = ~np.any(np.isnan(cov), axis=1) if cov.ndim > 1 else ~np.isnan(cov.ravel())
        d_clean = depth[mask]
        p_clean = log2_p[mask]
        c_clean = mc[mask]
        cov_clean = cov[mask]

        rp, pp = partial_corr(d_clean, p_clean, cov_clean)
        rc, pc = partial_corr(d_clean, c_clean, cov_clean)
        print(f'  P ~ depth | {label:30s}: r = {rp:+.3f}, p = {pp:.3e} (n={len(d_clean)})')
        partial_results_p.append({
            'control': label, 'r': round(rp, 3), 'p': f'{pp:.3e}', 'n': int(len(d_clean))
        })
        partial_results_c.append({
            'control': label, 'r': round(rc, 3), 'p': f'{pc:.3e}', 'n': int(len(d_clean))
        })

    print('\n--- Partial correlations: C-depth ---')
    for res in partial_results_c:
        print(f'  C ~ depth | {res["control"]:30s}: r = {res["r"]:+.3f}, p = {res["p"]}')

    results['partial_P'] = partial_results_p
    results['partial_C'] = partial_results_c

    # ========================================================
    # 3. Mediation analysis (Table 1)
    # ========================================================
    print('\n--- Mediation: individual covariates on P-depth ---')
    mediators = [
        ('Spatial frequency', sf.reshape(-1, 1)),
        ('Half-width', hw.reshape(-1, 1)),
        ('LHI', lhi.reshape(-1, 1)),
        ('SNR', snr.reshape(-1, 1)),
        ('OSI', osi.reshape(-1, 1)),
    ]

    mediation_results = []
    for label, cov in mediators:
        mask = ~np.isnan(cov.ravel())
        d_clean = depth[mask]
        p_clean = log2_p[mask]
        cov_clean = cov[mask]

        rp, pp = partial_corr(d_clean, p_clean, cov_clean)
        base_r = abs(r_p)
        drop = abs(r_p) - abs(rp)
        pct = (drop / base_r) * 100 if base_r != 0 else 0
        print(f'  Control {label:20s}: partial r = {rp:+.3f}, drop = {drop:.3f}, '
              f'% explained = {pct:.1f}%')
        mediation_results.append({
            'covariate': label, 'partial_r': round(rp, 3),
            'p': f'{pp:.3e}', 'pct_explained': round(pct, 1),
        })

    results['mediation'] = mediation_results

    # ========================================================
    # 4. Bootstrap CIs (10,000 iterations)
    # ========================================================
    print('\n--- Bootstrap CIs (10,000 iterations) ---')
    np.random.seed(42)
    n_boot = 10000
    boot_r_p, boot_r_c = [], []
    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        rp, _ = pearsonr(depth[idx], log2_p[idx])
        rc, _ = pearsonr(depth[idx], mc[idx])
        boot_r_p.append(rp)
        boot_r_c.append(rc)

    boot_r_p = np.array(boot_r_p)
    boot_r_c = np.array(boot_r_c)

    ci_p = (np.percentile(boot_r_p, 2.5), np.percentile(boot_r_p, 97.5))
    ci_c = (np.percentile(boot_r_c, 2.5), np.percentile(boot_r_c, 97.5))

    print(f'  P ~ depth: r = {r_p:+.3f}  95% CI [{ci_p[0]:+.3f}, {ci_p[1]:+.3f}]')
    print(f'  C ~ depth: r = {r_c:+.3f}  95% CI [{ci_c[0]:+.3f}, {ci_c[1]:+.3f}]')

    results['bootstrap'] = {
        'P_depth': {'r': round(r_p, 3), 'ci_low': round(ci_p[0], 3), 'ci_high': round(ci_p[1], 3)},
        'C_depth': {'r': round(r_c, 3), 'ci_low': round(ci_c[0], 3), 'ci_high': round(ci_c[1], 3)},
    }

    # ========================================================
    # 5. P-C anticorrelation
    # ========================================================
    r_pc, p_pc = pearsonr(log2_p, mc)
    print(f'\nP-C anticorrelation: r = {r_pc:+.3f}, p = {p_pc:.2e}')
    results['pc_anticorr'] = {'r': round(r_pc, 3), 'p': f'{p_pc:.2e}'}

    # ========================================================
    # 6. Covariate-depth correlations
    # ========================================================
    print('\n--- Covariate-depth correlations ---')
    for label, vals in [('SF', sf), ('Half-width', hw), ('LHI', lhi),
                        ('OSI', osi), ('Radius', radius), ('SNR', snr)]:
        mask = ~np.isnan(vals)
        r_val, p_val = pearsonr(depth[mask], vals[mask])
        print(f'  {label:15s} ~ depth: r = {r_val:+.3f}, p = {p_val:.2e}')

    # ========================================================
    # Save to JSON
    # ========================================================
    out_path = ROOT / 'paper' / 'manuscript_stats.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nResults saved to {out_path}')

    # ========================================================
    # Print LaTeX-ready table rows
    # ========================================================
    print('\n' + '=' * 60)
    print('LaTeX Table 1 (Mediation):')
    print('=' * 60)
    print(f'    None (zero-order) & ${r_p:+.3f}$ & --- \\\\')
    for m in mediation_results:
        print(f'    + {m["covariate"]} & ${m["partial_r"]:+.3f}$ & {m["pct_explained"]:.0f}\\% \\\\')

    print('\n' + '=' * 60)
    print('LaTeX Table 2 (Partial correlations):')
    print('=' * 60)
    for res in partial_results_p:
        print(f'    $P$ & {res["control"]} & ${res["r"]:+.3f}$ & ${res["p"]}$ \\\\')
    for res in partial_results_c:
        print(f'    $C$ & {res["control"]} & ${res["r"]:+.3f}$ & ${res["p"]}$ \\\\')


if __name__ == '__main__':
    main()
