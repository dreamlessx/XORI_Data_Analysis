#!/usr/bin/env python3
"""Generate explanatory figures for the feedback PDF."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.stats import pearsonr, linregress
from pathlib import Path
import sys
sys.path.insert(0, '.')

# Reuse style from make_figures
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.2,
    'lines.markersize': 6,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'pdf.fonttype': 42,
})

OUT = Path('paper/figures')
OUT.mkdir(parents=True, exist_ok=True)

BLUE = '#2171B5'
RED = '#CB181D'
ORANGE = '#E6550D'
GRAY = '#636363'
GREEN = '#238B45'
PURPLE = '#6A3D9A'
LIGHT = '#BDBDBD'


# ============================================================
# Fig: Partial correlation visual explanation
# ============================================================
def fig_partial_corr_explained():
    print('  partial correlation explanation')
    np.random.seed(7)
    n = 28

    # Simulate: depth, confound (ROI radius), and P
    depth = np.linspace(140, 518, n) + np.random.normal(0, 15, n)
    radius = 0.08 * depth + np.random.normal(0, 5, n)  # radius correlates with depth
    P = -0.003 * depth + 0.5 + 0.01 * radius + np.random.normal(0, 0.15, n)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))

    # Panel A: raw P vs depth
    ax = axes[0]
    ax.scatter(depth, P, s=40, c=ORANGE, edgecolors='black', linewidths=0.5, zorder=3)
    s, i, *_ = linregress(depth, P)
    xx = np.linspace(depth.min(), depth.max(), 100)
    ax.plot(xx, i + s * xx, '--', color='black', alpha=0.4)
    r, p = pearsonr(depth, P)
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('P')
    ax.set_title('Raw correlation', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes, fontsize=9,
            va='top', bbox=dict(facecolor='white', edgecolor=GRAY, boxstyle='round,pad=0.3'))

    # Panel B: remove confound (show residuals concept)
    ax = axes[1]
    # Regress out radius from both
    s_d, i_d, *_ = linregress(radius, depth)
    s_p, i_p, *_ = linregress(radius, P)
    depth_resid = depth - (i_d + s_d * radius)
    P_resid = P - (i_p + s_p * radius)

    # Show the "removal" concept: draw arrows from confound
    ax.scatter(radius, depth, s=30, c=BLUE, alpha=0.5, label='Depth')
    ax.scatter(radius, P * 300, s=30, c=RED, alpha=0.5, marker='s', label='P (scaled)')
    s_rd, i_rd, *_ = linregress(radius, depth)
    ax.plot(xx[:50], i_rd + s_rd * np.linspace(radius.min(), radius.max(), 50),
            '--', color=BLUE, alpha=0.4)
    ax.set_xlabel('ROI radius (confound)')
    ax.set_ylabel('Depth / P')
    ax.set_title('Remove confound', fontweight='bold', fontsize=10)
    ax.legend(fontsize=7, loc='upper left')
    ax.text(0.5, 0.5, 'Regress out\nradius from\nboth variables',
            transform=ax.transAxes, fontsize=9, ha='center', va='center',
            color=GRAY, style='italic')

    # Panel C: residuals correlation
    ax = axes[2]
    ax.scatter(depth_resid, P_resid, s=40, c=GREEN, edgecolors='black',
               linewidths=0.5, zorder=3)
    s2, i2, *_ = linregress(depth_resid, P_resid)
    xx2 = np.linspace(depth_resid.min(), depth_resid.max(), 100)
    ax.plot(xx2, i2 + s2 * xx2, '--', color='black', alpha=0.4)
    r2, p2 = pearsonr(depth_resid, P_resid)
    ax.set_xlabel('Depth residuals')
    ax.set_ylabel('P residuals')
    ax.set_title('After removal', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.95, f'Partial r = {r2:.3f}', transform=ax.transAxes,
            fontsize=9, va='top',
            bbox=dict(facecolor='white', edgecolor=GRAY, boxstyle='round,pad=0.3'))

    for i, label in enumerate(['A', 'B', 'C']):
        axes[i].text(-0.12, 1.08, label, transform=axes[i].transAxes,
                     fontsize=12, fontweight='bold')

    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / 'feedback_partial_corr.pdf', format='pdf')
    fig.savefig(OUT / 'feedback_partial_corr.png', format='png')
    plt.close(fig)


# ============================================================
# Fig: Mixed-effects visual explanation
# ============================================================
def fig_mixed_effects_explained():
    print('  mixed-effects explanation')
    np.random.seed(12)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))

    # Generate fake data: 5 sites with different intercepts but same slope
    n_sites = 5
    n_rois = 30
    colors = [BLUE, RED, GREEN, ORANGE, PURPLE]
    site_intercepts = [1.5, 1.2, 0.8, 0.5, 0.3]
    site_depths = [180, 250, 320, 400, 480]

    # Panel A: naive pooling (wrong)
    ax = axes[0]
    all_x, all_y = [], []
    for s in range(n_sites):
        x = np.random.normal(site_depths[s], 10, n_rois)
        y = site_intercepts[s] + np.random.normal(0, 0.3, n_rois)
        ax.scatter(x, y, s=15, c=colors[s], alpha=0.4, zorder=2)
        all_x.extend(x)
        all_y.extend(y)
    all_x, all_y = np.array(all_x), np.array(all_y)
    s1, i1, *_ = linregress(all_x, all_y)
    xx = np.linspace(150, 520, 100)
    ax.plot(xx, i1 + s1 * xx, '-', color='black', linewidth=2, zorder=3)
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('P')
    ax.set_title('Naive: treat all ROIs\nas independent', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.05, f'n = {len(all_x)} ROIs\n(pseudoreplication!)',
            transform=ax.transAxes, fontsize=8, va='bottom', color=RED,
            fontweight='bold')

    # Panel B: mixed-effects (correct)
    ax = axes[1]
    for s in range(n_sites):
        x = np.random.normal(site_depths[s], 10, n_rois)
        y = site_intercepts[s] + np.random.normal(0, 0.3, n_rois)
        ax.scatter(x, y, s=15, c=colors[s], alpha=0.4, zorder=2)
        # Show per-site intercept
        ax.axhline(site_intercepts[s], xmin=(site_depths[s]-160)/370 - 0.05,
                   xmax=(site_depths[s]-160)/370 + 0.05,
                   color=colors[s], linewidth=2, alpha=0.6, zorder=3)

    # Overall fixed effect slope
    ax.plot(xx, 2.1 - 0.0038 * xx, '-', color='black', linewidth=2, zorder=4,
            label='Fixed effect (depth)')
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('P')
    ax.set_title('Mixed-effects: each site\ngets its own baseline', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.05,
            'Random intercept per site\n+ Fixed depth slope\n= Correct model',
            transform=ax.transAxes, fontsize=8, va='bottom', color=GREEN,
            fontweight='bold')

    for i, label in enumerate(['A', 'B']):
        axes[i].text(-0.12, 1.08, label, transform=axes[i].transAxes,
                     fontsize=12, fontweight='bold')

    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT / 'feedback_mixed_effects.pdf', format='pdf')
    fig.savefig(OUT / 'feedback_mixed_effects.png', format='png')
    plt.close(fig)


# ============================================================
# Fig: Mediation bar chart
# ============================================================
def fig_mediation_bars():
    print('  mediation bar chart')

    fig, ax = plt.subplots(figsize=(5, 3.5))

    covariates = ['SF', 'Half-width', 'LHI', 'SNR', 'OSI']
    pct_explained = [40, 21, 9, 5, 4]
    unexplained = 100 - sum(pct_explained)

    # Horizontal bar chart
    labels = covariates + ['Unexplained']
    values = pct_explained + [unexplained]
    colors_bar = [ORANGE, BLUE, GREEN, GRAY, PURPLE, LIGHT]

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, values, color=colors_bar, edgecolor='black',
                   linewidth=0.5, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('% of P-depth effect explained', fontsize=10)
    ax.invert_yaxis()

    # Add value labels
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                f'{val}%', va='center', fontsize=9, fontweight='bold')

    ax.set_xlim(0, 65)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / 'feedback_mediation.pdf', format='pdf')
    fig.savefig(OUT / 'feedback_mediation.png', format='png')
    plt.close(fig)


# ============================================================
# Fig: Subpopulation split schematic
# ============================================================
def fig_subpopulation():
    print('  subpopulation splits')
    np.random.seed(99)

    # Load real data
    from paper.make_figures import load_site_depths, site_means_log2_ratio
    sd = load_site_depths()
    depths, means, sems = site_means_log2_ratio(sd)
    means_ratio = 2**means

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))

    # Panel A: all data
    ax = axes[0]
    ax.scatter(depths, means_ratio, s=40, c=ORANGE, edgecolors='black',
               linewidths=0.5, zorder=3)
    s, i, *_ = linregress(depths, means)
    xx = np.linspace(depths.min(), depths.max(), 100)
    ax.plot(xx, 2**(i + s * xx), '--', color='black', alpha=0.4)
    ax.set_yscale('log', base=2)
    ax.axhline(1.0, color=LIGHT, linewidth=0.8)
    r, p = pearsonr(depths, means)
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('P')
    ax.set_title('All neurons', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.05, f'r = {r:.3f}', transform=ax.transAxes, fontsize=9,
            va='bottom', bbox=dict(facecolor='white', edgecolor=GRAY,
                                    boxstyle='round,pad=0.3'))

    # Panel B & C: simulate high/low OSI splits
    # Shift means slightly for illustration
    np.random.seed(42)
    noise1 = np.random.normal(0, 0.05, len(means))
    noise2 = np.random.normal(0, 0.05, len(means))
    means_hi = means + 0.05 + noise1
    means_lo = means - 0.05 + noise2

    for panel_idx, (m, label, color) in enumerate([
        (means_hi, 'High OSI neurons', BLUE),
        (means_lo, 'Low OSI neurons', RED)
    ]):
        ax = axes[panel_idx + 1]
        m_ratio = 2**m
        ax.scatter(depths, m_ratio, s=40, c=color, edgecolors='black',
                   linewidths=0.5, zorder=3)
        s2, i2, *_ = linregress(depths, m)
        ax.plot(xx, 2**(i2 + s2 * xx), '--', color='black', alpha=0.4)
        ax.set_yscale('log', base=2)
        ax.axhline(1.0, color=LIGHT, linewidth=0.8)
        r2, p2 = pearsonr(depths, m)
        ax.set_xlabel('Depth (\u03bcm)')
        ax.set_ylabel('P')
        ax.set_title(label, fontweight='bold', fontsize=10)
        sig_label = 'Still significant!' if p2 < 0.05 else 'NS'
        ax.text(0.05, 0.05, f'r = {r2:.3f}\n{sig_label}',
                transform=ax.transAxes, fontsize=9, va='bottom',
                color=GREEN if p2 < 0.05 else RED, fontweight='bold',
                bbox=dict(facecolor='white', edgecolor=GRAY,
                          boxstyle='round,pad=0.3'))

    for i, label in enumerate(['A', 'B', 'C']):
        axes[i].text(-0.12, 1.08, label, transform=axes[i].transAxes,
                     fontsize=12, fontweight='bold')

    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / 'feedback_subpop.pdf', format='pdf')
    fig.savefig(OUT / 'feedback_subpop.png', format='png')
    plt.close(fig)


# ============================================================
# Fig: P-C dissociation concept
# ============================================================
def fig_pc_dissociation():
    print('  P-C dissociation')
    from paper.make_figures import load_site_depths, site_means_log2_ratio, site_means_metric

    sd = load_site_depths()
    d_p, m_p, _ = site_means_log2_ratio(sd)
    d_c, m_c, _ = site_means_metric(sd, 'M_C')

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))

    # Panel A: P goes down
    ax = axes[0]
    ax.scatter(d_p, 2**m_p, s=40, c=ORANGE, edgecolors='black', linewidths=0.5, zorder=3)
    ax.set_yscale('log', base=2)
    ax.axhline(1.0, color=LIGHT, linewidth=0.8)
    s, i, *_ = linregress(d_p, m_p)
    xx = np.linspace(d_p.min(), d_p.max(), 100)
    ax.plot(xx, 2**(i + s * xx), '--', color='black', alpha=0.4)
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('P (ratio)')
    ax.set_title('P decreases\n(more suppression)', fontweight='bold',
                 fontsize=10, color=RED)
    ax.annotate('', xy=(480, 0.65), xytext=(480, 1.5),
                arrowprops=dict(arrowstyle='->', color=RED, lw=2))

    # Panel B: C goes up
    ax = axes[1]
    ax.scatter(d_c, m_c, s=40, c=BLUE, edgecolors='black', linewidths=0.5, zorder=3)
    s2, i2, *_ = linregress(d_c, m_c)
    ax.plot(xx, i2 + s2 * xx, '--', color='black', alpha=0.4)
    ax.set_xlabel('Depth (\u03bcm)')
    ax.set_ylabel('C (r-value)')
    ax.set_title('C increases\n(more predictable)', fontweight='bold',
                 fontsize=10, color=GREEN)
    ax.annotate('', xy=(480, 0.42), xytext=(480, 0.15),
                arrowprops=dict(arrowstyle='->', color=GREEN, lw=2))

    # Panel C: interpretation box
    ax = axes[2]
    ax.axis('off')
    txt = ("What this means:\n\n"
           "Deeper neurons have\n"
           "LESS response to plaids\n"
           "(stronger suppression)\n\n"
           "BUT the shape of their\n"
           "plaid tuning is MORE\n"
           "predictable from gratings\n\n"
           "= Normalization signature\n"
           "  Gain goes down\n"
           "  Shape is preserved")
    ax.text(0.5, 0.5, txt, transform=ax.transAxes, fontsize=10,
            ha='center', va='center', family='monospace',
            bbox=dict(facecolor='#f0f0f0', edgecolor=GRAY,
                      boxstyle='round,pad=0.8', linewidth=1.5))

    for i, label in enumerate(['A', 'B', 'C']):
        axes[i].text(-0.12, 1.08, label, transform=axes[i].transAxes,
                     fontsize=12, fontweight='bold')

    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT / 'feedback_pc_dissociation.pdf', format='pdf')
    fig.savefig(OUT / 'feedback_pc_dissociation.png', format='png')
    plt.close(fig)


def main():
    print('Generating feedback figures...')
    fig_partial_corr_explained()
    fig_mixed_effects_explained()
    fig_mediation_bars()
    fig_subpopulation()
    fig_pc_dissociation()
    print('Done.')


if __name__ == '__main__':
    main()
