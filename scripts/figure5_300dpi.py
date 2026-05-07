"""Regenerate Figure 5 at 300 DPI.

Reads the v7.0 canonical-fits CSV from data/ and writes Figure 5
(ΔBIC distribution histogram) to figures/.

Run from the package root:
    cd halo_shells_v7.0/
    python3 scripts/figure5_300dpi.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Package-relative paths
HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
INPUT_CSV = os.path.join(PACKAGE_ROOT, 'data', 'sparc_T2-T9_canonical_fits.csv')
OUTPUT_DIR = os.path.join(PACKAGE_ROOT, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(INPUT_CSV)
dbic_burk = df['fw_best_bic'] - df['burk_bic']
dbic_nfw = df['fw_best_bic'] - df['nfw_bic']

fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

def plot_hist(ax, dbic, title):
    dbic_capped = np.clip(dbic, -100, 30)
    bins = np.unique(np.concatenate([
        np.linspace(-100, -10, 7),
        np.linspace(-10, -2, 5),
        np.linspace(-2, 2, 3),
        np.linspace(2, 30, 5)
    ]))
    counts, edges = np.histogram(dbic_capped, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)
    
    colors = []
    for c in centers:
        if c < -10:    colors.append('#1a5490')
        elif c < -6:   colors.append('#3a7cba')
        elif c < -2:   colors.append('#7eaed3')
        elif c < 2:    colors.append('#d3d3d3')
        elif c < 6:    colors.append('#e8a87a')
        elif c < 10:   colors.append('#d97a3e')
        else:          colors.append('#a83c0c')
    
    ax.bar(centers, counts, width=widths, color=colors,
           edgecolor='black', linewidth=0.5, align='center')
    
    for x in [-10, -6, -2, 2, 6, 10]:
        ax.axvline(x, color='black', linestyle=':', linewidth=0.7, alpha=0.5)
    
    txt = [f"FW very strong: {(dbic < -10).sum()}",
           f"FW strong: {((dbic >= -10) & (dbic < -6)).sum()}",
           f"FW positive: {((dbic >= -6) & (dbic < -2)).sum()}",
           f"Inconclusive: {((dbic >= -2) & (dbic <= 2)).sum()}",
           f"other positive: {((dbic > 2) & (dbic < 6)).sum()}",
           f"other strong: {((dbic >= 6) & (dbic < 10)).sum()}",
           f"other very strong: {(dbic >= 10).sum()}"]
    
    ax.text(0.97, 0.97, '\n'.join(txt),
            transform=ax.transAxes, va='top', ha='right',
            fontsize=9, family='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      edgecolor='gray', alpha=0.95))
    
    ax.set_xlabel('Delta BIC = BIC(Framework) - BIC(comparison)', fontsize=11)
    ax.set_ylabel('Number of galaxies', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    ax.set_xlim(-105, 32)
    
    if (dbic >= 10).sum() == 0:
        ax.text(20, ax.get_ylim()[1]*0.45, '0 galaxies',
                fontsize=10, color='darkred', fontweight='bold', ha='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='mistyrose',
                          edgecolor='darkred'))

plot_hist(axes[0], dbic_burk, 'Framework vs Burkert')
plot_hist(axes[1], dbic_nfw, 'Framework vs NFW')

legend_elements = [
    Patch(facecolor='#1a5490', label='FW very strong (Delta BIC < -10)'),
    Patch(facecolor='#3a7cba', label='FW strong'),
    Patch(facecolor='#7eaed3', label='FW positive'),
    Patch(facecolor='#d3d3d3', label='Inconclusive'),
    Patch(facecolor='#e8a87a', label='other positive'),
    Patch(facecolor='#d97a3e', label='other strong'),
    Patch(facecolor='#a83c0c', label='other very strong'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4,
           fontsize=9, framealpha=0.95, bbox_to_anchor=(0.5, -0.06))

fig.suptitle('Delta BIC distribution across SPARC T=2-9 sample (n=102) - '
             '0 galaxies show very strong evidence against framework',
             fontsize=12, fontweight='bold', y=1.05)

plt.savefig(os.path.join(OUTPUT_DIR, 'figure5_dbic_histogram.png'),
            dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'figure5_dbic_histogram.pdf'),
            bbox_inches='tight')
print(f"Saved Figure 5 at 300 DPI (PNG + PDF) to {OUTPUT_DIR}")
plt.close()
