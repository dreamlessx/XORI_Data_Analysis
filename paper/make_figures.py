#!/usr/bin/env python3
"""
Generate publication-quality figures for the XORI manuscript.
Output: paper/figures/*.pdf (vector) + paper/figures/*.png (300 dpi)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
from scipy.stats import pearsonr, linregress
from pathlib import Path

# ============================================================
# Style configuration
# ============================================================
SINGLE_COL = 3.5   # inches
DOUBLE_COL = 7.0   # inches
DPI = 300

COLORS = {
    'primary': '#2171B5',      # blue
    'secondary': '#CB181D',    # red
    'accent': '#E6550D',       # orange
    'neutral': '#636363',      # gray
    'light': '#BDBDBD',        # light gray
    'green': '#238B45',
    'purple': '#6A3D9A',
}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'lines.linewidth': 1.0,
    'lines.markersize': 5,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'pdf.fonttype': 42,       # editable text in PDF
    'ps.fonttype': 42,
})

OUT = Path('paper/figures')
OUT.mkdir(parents=True, exist_ok=True)

ROOT = Path('.')

# ============================================================
# Data loading helpers
# ============================================================

def load_site_depths():
    depths = {}
    with open(ROOT / 'raw_data/bm_data/site_depth.txt') as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                depths[parts[0]] = float(parts[1])
    return depths


def load_metrics(site_name):
    path = ROOT / f'metric_data/all_roi/metrics_{site_name}.txt'
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


def load_covariate(filename):
    """Load a covariate file (roi_sf.txt, roi_osi.txt, etc.)"""
    data = {}
    with open(ROOT / 'raw_data/bm_data' / filename) as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    site = int(float(parts[0]))
                    roi = int(float(parts[1]))
                    val = float(parts[3])
                    if site not in data:
                        data[site] = {}
                    data[site][roi] = val
                except (ValueError, IndexError):
                    continue
    return data


def load_tuning_curve(site_name, roi_idx):
    """Load grating, plaid, and baseline from tuning curve file."""
    path = ROOT / f'raw_data/tc_data/{site_name}/roi_{roi_idx:04d}'
    if not path.exists():
        return None, None, None

    blocks = []
    current = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('/newplot'):
                if current:
                    blocks.append(current)
                current = []
            elif line.startswith('/plotname'):
                continue
            else:
                parts = line.split()
                if len(parts) >= 2:
                    current.append([float(x) for x in parts])
    if current:
        blocks.append(current)

    grating = np.array(blocks[0]) if len(blocks) > 0 else None
    plaid = np.array(blocks[1]) if len(blocks) > 1 else None
    baseline_arr = np.array(blocks[2]) if len(blocks) > 2 else None
    baseline = baseline_arr[0, 1] if baseline_arr is not None else 0.0

    return grating, plaid, baseline


def site_means_log2_ratio(site_depths):
    """Compute site-mean log2(M_S_ratio) for all sites."""
    depths, means, sems = [], [], []
    for site, depth in sorted(site_depths.items(), key=lambda x: x[1]):
        try:
            m = load_metrics(site)
        except FileNotFoundError:
            continue
        pos = m['M_S_ratio'][m['M_S_ratio'] > 0]
        if len(pos) > 0:
            log_vals = np.log2(pos)
            depths.append(depth)
            means.append(np.mean(log_vals))
            sems.append(np.std(log_vals, ddof=1) / np.sqrt(len(log_vals))
                        if len(log_vals) > 1 else 0)
    return np.array(depths), np.array(means), np.array(sems)


def site_means_metric(site_depths, key):
    """Compute site-mean of a metric for all sites."""
    depths, means, sems = [], [], []
    for site, depth in sorted(site_depths.items(), key=lambda x: x[1]):
        try:
            m = load_metrics(site)
        except FileNotFoundError:
            continue
        vals = m[key]
        if len(vals) > 0:
            depths.append(depth)
            means.append(np.mean(vals))
            sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals))
                        if len(vals) > 1 else 0)
    return np.array(depths), np.array(means), np.array(sems)


def panel_label(ax, label, x=-0.15, y=1.05):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top', ha='left')


def stats_text(ax, r, p, n, loc='lower left'):
    txt = f'r = {r:+.3f}\np = {p:.2g}\nn = {n}'
    # Place outside the data area using axes margins
    anchors = {
        'lower left':  (0.03, 0.03, 'bottom', 'left'),
        'lower right': (0.97, 0.03, 'bottom', 'right'),
        'upper left':  (0.03, 0.97, 'top',    'left'),
        'upper right': (0.97, 0.97, 'top',    'right'),
    }
    x, y, va, ha = anchors[loc]
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=7,
            va=va, ha=ha, zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#999999', alpha=1.0, linewidth=0.5))


def add_regression(ax, x, y, color='black', alpha=0.4):
    slope, intercept, *_ = linregress(x, y)
    xx = np.linspace(x.min(), x.max(), 200)
    ax.plot(xx, intercept + slope * xx, '--', color=color,
            linewidth=1.0, alpha=alpha, zorder=2)


def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', format='pdf')
    fig.savefig(OUT / f'{name}.png', format='png')
    plt.close(fig)
    print(f'  [saved] {name}.pdf + .png')


# ============================================================
# Figure 1: Example tuning curves (shallow vs deep)
# ============================================================

def fig1_example_tuning_curves():
    print('Figure 1: Example tuning curves')
    sd = load_site_depths()

    # Pick a shallow site and a deep site
    shallow_site = 'site002'  # 242 um
    deep_site = 'site020'     # 482 um

    # Find good example ROIs (high SNR, clear suppression/facilitation)
    def find_example_roi(site_name, want_facilitation=True):
        m = load_metrics(site_name)
        snr_mask = m['SNR_g'] > 1.5
        if want_facilitation:
            score = m['M_S_ratio'] * m['SNR_g']
            score[~snr_mask] = -999
            idx = np.argmax(score)
        else:
            score = (1.0 / np.clip(m['M_S_ratio'], 0.01, None)) * m['SNR_g']
            score[~snr_mask] = -999
            idx = np.argmax(score)
        return int(m['ROI'][idx])

    shallow_roi = find_example_roi(shallow_site, want_facilitation=True)
    deep_roi = find_example_roi(deep_site, want_facilitation=False)

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 4.5))

    for row, (site, roi, label_depth) in enumerate([
        (shallow_site, shallow_roi, f'{sd[shallow_site]:.0f}'),
        (deep_site, deep_roi, f'{sd[deep_site]:.0f}'),
    ]):
        grating, plaid, baseline = load_tuning_curve(site, roi)
        if grating is None:
            continue

        dirs = grating[:, 0]
        g_mean = grating[:, 1] - baseline
        g_sem = grating[:, 2] if grating.shape[1] > 2 else np.zeros_like(g_mean)
        p_mean = plaid[:, 1] - baseline
        p_sem = plaid[:, 2] if plaid.shape[1] > 2 else np.zeros_like(p_mean)

        # Build linear prediction: shift grating by -90 and sum
        n_dirs = len(dirs)
        shift = n_dirs // 4  # 90 deg = 3 steps of 30 deg
        pred = g_mean + np.roll(g_mean, -shift)

        # Panel A: grating tuning curve
        ax = axes[row, 0]
        ax.errorbar(dirs, g_mean, yerr=g_sem, fmt='o-',
                    color=COLORS['primary'], markersize=4, linewidth=1.0,
                    capsize=2, capthick=0.6, label='Grating')
        ax.set_xlabel('Direction (deg)')
        ax.set_ylabel(r'$\Delta$F (a.u.)')
        ax.set_xlim(-15, 345)
        ax.set_xticks([0, 90, 180, 270])
        ax.axhline(0, color=COLORS['light'], linewidth=0.5, zorder=0)
        depth_str = 'Shallow' if row == 0 else 'Deep'
        ax.text(0.98, 0.02, f'{depth_str} ({label_depth} \u03bcm), ROI {roi}',
                transform=ax.transAxes, fontsize=6, va='bottom', ha='right',
                color=COLORS['neutral'])

        # Panel B: plaid observed vs prediction
        ax = axes[row, 1]
        ax.plot(dirs, pred, '--', color=COLORS['neutral'], linewidth=1.2,
                label='Linear prediction', zorder=2)
        ax.errorbar(dirs, p_mean, yerr=p_sem, fmt='o-',
                    color=COLORS['secondary'], markersize=4, linewidth=1.0,
                    capsize=2, capthick=0.6, label='Plaid (observed)', zorder=3)
        ax.set_xlabel('Direction (deg)')
        ax.set_ylabel(r'$\Delta$F (a.u.)')
        ax.set_xlim(-15, 345)
        ax.set_xticks([0, 90, 180, 270])
        ax.axhline(0, color=COLORS['light'], linewidth=0.5, zorder=0)
        ax.legend(fontsize=6, loc='upper right', frameon=True,
                  edgecolor='#cccccc', fancybox=False)

        # Compute P for annotation — place in lower left to avoid data
        p_ratio = np.sum(p_mean) / np.sum(pred) if np.sum(pred) != 0 else 0
        ax.text(0.03, 0.03, f'P = {p_ratio:.2f}',
                transform=ax.transAxes, fontsize=8, fontweight='bold',
                va='bottom', ha='left',
                color=COLORS['green'] if p_ratio > 1 else COLORS['secondary'])

    panel_label(axes[0, 0], 'A')
    panel_label(axes[0, 1], 'B')
    panel_label(axes[1, 0], 'C')
    panel_label(axes[1, 1], 'D')

    fig.tight_layout(h_pad=1.5, w_pad=1.0)
    save(fig, 'fig1_tuning_curves')


# ============================================================
# Figure 3: P and C vs depth (main result)
# ============================================================

def fig3_depth():
    print('Figure 3: P and C vs depth')
    sd = load_site_depths()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.0))

    # --- Panel A: P vs depth (log axis) ---
    depths, means, sems = site_means_log2_ratio(sd)
    # Convert back for display on log axis
    means_ratio = 2**means
    sems_lower = 2**(means - sems)
    sems_upper = 2**(means + sems)
    yerr = [means_ratio - sems_lower, sems_upper - means_ratio]

    ax1.errorbar(depths, means_ratio, yerr=yerr, fmt='o',
                 color=COLORS['accent'], markersize=5,
                 markeredgecolor='black', markeredgewidth=0.5,
                 ecolor='#666666', elinewidth=0.6, capsize=2, capthick=0.5,
                 zorder=3)
    ax1.set_yscale('log', base=2)

    # Regression in log space, display in data coords
    add_regression(ax1, depths, means_ratio, alpha=0.0)  # invisible
    slope, intercept, *_ = linregress(depths, means)
    xx = np.linspace(depths.min(), depths.max(), 200)
    yy_log = intercept + slope * xx
    ax1.plot(xx, 2**yy_log, '--', color='black', linewidth=1.0, alpha=0.4, zorder=2)

    ax1.axhline(1.0, color=COLORS['light'], linewidth=0.8, linestyle='-', zorder=1)

    # Tick labels as actual ratio values
    tick_vals = [0.5, 0.59, 0.71, 0.84, 1.0, 1.19, 1.41, 1.68, 2.0]
    ax1.yaxis.set_major_locator(FixedLocator(tick_vals))
    ax1.yaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f'{x:.2f}' if x < 1 else f'{x:.1f}' if x != int(x) else f'{int(x)}'))

    ax1.set_xlabel('Depth (\u03bcm)')
    ax1.set_ylabel('P (ratio)')
    ax1.set_ylim(0.48, 2.1)

    r, p = pearsonr(depths, means)
    stats_text(ax1, r, p, len(depths), 'lower left')

    # --- Panel B: C vs depth ---
    depths_c, means_c, sems_c = site_means_metric(sd, 'M_C')

    ax2.errorbar(depths_c, means_c, yerr=sems_c, fmt='o',
                 color=COLORS['primary'], markersize=5,
                 markeredgecolor='black', markeredgewidth=0.5,
                 ecolor='#666666', elinewidth=0.6, capsize=2, capthick=0.5,
                 zorder=3)
    add_regression(ax2, depths_c, means_c)

    ax2.set_xlabel('Depth (\u03bcm)')
    ax2.set_ylabel('C (r-value)')

    r_c, p_c = pearsonr(depths_c, means_c)
    stats_text(ax2, r_c, p_c, len(depths_c), 'lower right')

    panel_label(ax1, 'A')
    panel_label(ax2, 'B')

    fig.tight_layout(w_pad=2.0)
    save(fig, 'fig3_depth')


# ============================================================
# Figure 4: SF vs depth
# ============================================================

def fig4_sf():
    print('Figure 4: Spatial frequency vs depth')
    sd = load_site_depths()
    sf_data = load_covariate('roi_sf.txt')

    fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.0))

    depths, means, sems = [], [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        site_num = int(site.replace('site', ''))
        if site_num not in sf_data:
            continue
        vals = list(sf_data[site_num].values())
        if vals:
            depths.append(depth)
            means.append(np.mean(vals))
            sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals))
                        if len(vals) > 1 else 0)

    depths, means, sems = np.array(depths), np.array(means), np.array(sems)

    ax.errorbar(depths, means, yerr=sems, fmt='o',
                color=COLORS['accent'], markersize=5,
                markeredgecolor='black', markeredgewidth=0.5,
                ecolor='#666666', elinewidth=0.6, capsize=2, capthick=0.5,
                zorder=3)
    add_regression(ax, depths, means)

    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('Spatial frequency (cyc/deg)')

    r, p = pearsonr(depths, means)
    stats_text(ax, r, p, len(depths), 'lower left')

    fig.tight_layout()
    save(fig, 'fig4_sf')


# ============================================================
# Figure 5: Baseline vs depth
# ============================================================

def fig5_baseline():
    print('Figure 5: Baseline vs depth')
    sd = load_site_depths()

    fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.0))

    depths, means, sems = [], [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        try:
            m = load_metrics(site)
        except FileNotFoundError:
            continue
        # Baseline from tuning curve files - approximate from M_S_norm
        # Actually, use raw baseline from tc data
        baselines = []
        for roi_idx in m['ROI'][:50]:  # sample for speed
            _, _, bl = load_tuning_curve(site, int(roi_idx))
            if bl is not None:
                baselines.append(bl)
        if baselines:
            depths.append(depth)
            means.append(np.mean(baselines))
            sems.append(np.std(baselines, ddof=1) / np.sqrt(len(baselines)))

    depths, means, sems = np.array(depths), np.array(means), np.array(sems)

    ax.errorbar(depths, means, yerr=sems, fmt='o',
                color=COLORS['accent'], markersize=5,
                markeredgecolor='black', markeredgewidth=0.5,
                ecolor='#666666', elinewidth=0.6, capsize=2, capthick=0.5,
                zorder=3)
    add_regression(ax, depths, means)

    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('Baseline fluorescence (a.u.)')

    r, p = pearsonr(depths, means)
    stats_text(ax, r, p, len(depths), 'upper left')

    fig.tight_layout()
    save(fig, 'fig5_baseline')


# ============================================================
# Figures 6-8: Covariate with within-site correlations
# ============================================================

def fig_covariate_with_corr(fig_num, name, cov_file, cov_col_idx,
                            cov_label, metric_key, metric_label,
                            corr_metric_key='M_C'):
    """Generate a 2-panel figure: covariate vs depth + within-site r vs depth."""
    print(f'Figure {fig_num}: {name}')
    sd = load_site_depths()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.0))

    # Load covariate data
    cov_data = {}
    with open(ROOT / 'raw_data/bm_data' / cov_file) as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    site_num = int(float(parts[0]))
                    roi = int(float(parts[1]))
                    val = float(parts[cov_col_idx])
                    if site_num not in cov_data:
                        cov_data[site_num] = {}
                    cov_data[site_num][roi] = val
                except (ValueError, IndexError):
                    continue

    # Panel A: covariate vs depth
    cov_depths, cov_means, cov_sems = [], [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        site_num = int(site.replace('site', ''))
        if site_num not in cov_data:
            continue
        vals = list(cov_data[site_num].values())
        if vals:
            cov_depths.append(depth)
            cov_means.append(np.mean(vals))
            cov_sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals))
                            if len(vals) > 1 else 0)

    cov_depths = np.array(cov_depths)
    cov_means = np.array(cov_means)
    cov_sems = np.array(cov_sems)

    ax1.errorbar(cov_depths, cov_means, yerr=cov_sems, fmt='o',
                 color=COLORS['accent'], markersize=5,
                 markeredgecolor='black', markeredgewidth=0.5,
                 ecolor='#666666', elinewidth=0.6, capsize=2, capthick=0.5,
                 zorder=3)
    add_regression(ax1, cov_depths, cov_means)
    ax1.set_xlabel('Depth (\u03bcm)')
    ax1.set_ylabel(cov_label)

    r1, p1 = pearsonr(cov_depths, cov_means)
    # Place stats opposite to the data trend to avoid overlap
    if r1 < 0:
        # Data goes down-right, so lower left is empty
        stats_loc = 'lower left'
    else:
        # Data goes up-right, so lower right is empty
        stats_loc = 'lower right'
    stats_text(ax1, r1, p1, len(cov_depths), stats_loc)

    # Panel B: within-site correlation (covariate vs C) plotted against depth
    corr_depths, corr_rs, corr_ps = [], [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        site_num = int(site.replace('site', ''))
        if site_num not in cov_data:
            continue
        try:
            m = load_metrics(site)
        except FileNotFoundError:
            continue

        x_vals, y_vals = [], []
        for i, roi in enumerate(m['ROI']):
            if int(roi) in cov_data[site_num]:
                x_vals.append(cov_data[site_num][int(roi)])
                y_vals.append(m[corr_metric_key][i])

        if len(x_vals) >= 10:
            r_val, p_val = pearsonr(x_vals, y_vals)
            corr_depths.append(depth)
            corr_rs.append(r_val)
            corr_ps.append(p_val)

    corr_depths = np.array(corr_depths)
    corr_rs = np.array(corr_rs)
    corr_ps = np.array(corr_ps)

    sig = corr_ps < 0.05
    if np.any(sig):
        ax2.scatter(corr_depths[sig], corr_rs[sig], s=30,
                    color=COLORS['primary'], edgecolors='black',
                    linewidths=0.5, zorder=3, label='p < 0.05')
    if np.any(~sig):
        ax2.scatter(corr_depths[~sig], corr_rs[~sig], s=30,
                    facecolors='none', edgecolors=COLORS['primary'],
                    linewidths=0.8, zorder=3, label=r'p $\geq$ 0.05')

    ax2.axhline(0, color=COLORS['light'], linewidth=0.8, zorder=1)
    ax2.set_xlabel('Depth (\u03bcm)')
    ax2.set_ylabel(f'Within-site r ({metric_label} vs {cov_label.split()[0]})')

    # Count significant sites
    n_sig = int(np.sum(sig))

    # Place legend and sig count in the emptiest quadrant
    # Check where data is sparse to avoid overlap
    mid_r = np.median(corr_rs) if len(corr_rs) > 0 else 0
    if mid_r < 0:
        # Most data below zero — put legend top-left, sig count top-right
        legend_loc = 'upper left'
        sig_pos = (0.97, 0.97, 'top', 'right')
    else:
        # Most data above zero — put legend bottom-left, sig count bottom-right
        legend_loc = 'lower left'
        sig_pos = (0.97, 0.03, 'bottom', 'right')

    ax2.legend(fontsize=6, loc=legend_loc, frameon=True,
               edgecolor='#cccccc', fancybox=False)
    ax2.text(sig_pos[0], sig_pos[1], f'{n_sig}/{len(corr_rs)} sig.',
             transform=ax2.transAxes, fontsize=7, va=sig_pos[2], ha=sig_pos[3],
             color=COLORS['neutral'])

    panel_label(ax1, 'A')
    panel_label(ax2, 'B')

    fig.tight_layout(w_pad=2.0)
    save(fig, f'fig{fig_num}_{name.lower().replace(" ", "_")}')


# ============================================================
# Depth profile summary figure
# ============================================================

def fig_depth_profile():
    print('Figure 9: Depth profile summary')
    sd = load_site_depths()

    metrics_to_plot = [
        ('M_S_ratio', 'P', True),   # log2 for ratio
        ('M_C', 'C', False),
    ]

    sf_data = load_covariate('roi_sf.txt')
    osi_data = load_covariate('roi_osi.txt')

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 5.0))

    # Panel A: P vs depth
    ax = axes[0, 0]
    d, m, s = site_means_log2_ratio(sd)
    m_ratio = 2**m
    ax.scatter(d, m_ratio, s=25, color=COLORS['accent'],
               edgecolors='black', linewidths=0.5, zorder=3)
    add_regression(ax, d, m_ratio)
    ax.set_yscale('log', base=2)
    ax.axhline(1.0, color=COLORS['light'], linewidth=0.8, zorder=1)
    ax.yaxis.set_major_locator(FixedLocator([0.59, 0.71, 0.84, 1, 1.19, 1.41, 1.68]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.2f}' if x < 1 else f'{x:.1f}'))
    ax.set_ylabel('P')
    ax.set_xlabel('Depth (\u03bcm)')
    panel_label(ax, 'A')

    # Panel B: C vs depth
    ax = axes[0, 1]
    d_c, m_c, s_c = site_means_metric(sd, 'M_C')
    ax.scatter(d_c, m_c, s=25, color=COLORS['primary'],
               edgecolors='black', linewidths=0.5, zorder=3)
    add_regression(ax, d_c, m_c)
    ax.set_ylabel('C')
    ax.set_xlabel('Depth (\u03bcm)')
    panel_label(ax, 'B')

    # Panel C: SF vs depth
    ax = axes[1, 0]
    sf_depths, sf_means = [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        sn = int(site.replace('site', ''))
        if sn in sf_data:
            sf_depths.append(depth)
            sf_means.append(np.mean(list(sf_data[sn].values())))
    ax.scatter(sf_depths, sf_means, s=25, color=COLORS['green'],
               edgecolors='black', linewidths=0.5, zorder=3)
    add_regression(ax, np.array(sf_depths), np.array(sf_means))
    ax.set_ylabel('SF (cyc/deg)')
    ax.set_xlabel('Depth (\u03bcm)')
    panel_label(ax, 'C')

    # Panel D: OSI vs depth
    ax = axes[1, 1]
    osi_depths, osi_means = [], []
    for site, depth in sorted(sd.items(), key=lambda x: x[1]):
        sn = int(site.replace('site', ''))
        if sn in osi_data:
            osi_depths.append(depth)
            osi_means.append(np.mean(list(osi_data[sn].values())))
    ax.scatter(osi_depths, osi_means, s=25, color=COLORS['purple'],
               edgecolors='black', linewidths=0.5, zorder=3)
    add_regression(ax, np.array(osi_depths), np.array(osi_means))
    ax.set_ylabel('OSI')
    ax.set_xlabel('Depth (\u03bcm)')
    panel_label(ax, 'D')

    fig.tight_layout(h_pad=1.5, w_pad=2.0)
    save(fig, 'fig9_depth_profile')


# ============================================================
# Main
# ============================================================

def main():
    print('Generating publication figures...\n')

    fig1_example_tuning_curves()
    fig3_depth()
    fig4_sf()
    # fig5 is slow (reads many TC files); skip baseline for now
    # fig5_baseline()

    # Figs 6-8: covariate + within-site correlation panels
    fig_covariate_with_corr(
        6, 'Bandwidth', 'roi_hw_orth.txt', 3,
        'Half-width (deg)', 'M_C', 'C')

    fig_covariate_with_corr(
        7, 'LHI', 'roi_lhi.txt', 3,
        'LHI', 'M_C', 'C')

    fig_covariate_with_corr(
        8, 'OSI', 'roi_osi.txt', 3,
        'OSI', 'M_C', 'C')

    fig_depth_profile()

    print(f'\nAll figures saved to {OUT}/')


if __name__ == '__main__':
    main()
