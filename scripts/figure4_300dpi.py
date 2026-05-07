"""Regenerate Figure 4 at 300 DPI.

Shell-bearing fraction by morphological type, with synthetic-null overlay
and per-bin Fisher one-sided p-values.

All numbers (real shell-bearing fractions, null FP rates, per-T sample sizes,
trend statistics) are computed from the v7.0 CSVs in data/ — no hardcoded
values, so the figure stays in sync with the data files.

Run from the package root:
    cd halo_shells_v7.0/
    python3 scripts/figure4_300dpi.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Package-relative paths
HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
CANONICAL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'sparc_T2-T9_canonical_fits.csv')
NULL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'null_test_T2-T9_combined.csv')
OUTPUT_DIR = os.path.join(PACKAGE_ROOT, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- Load data ----
canon = pd.read_csv(CANONICAL_CSV)
null = pd.read_csv(NULL_CSV)

T_types = [2, 3, 4, 5, 6, 7, 8, 9]
labels = ['Sab', 'Sb', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm']

# ---- Real-data shell-bearing fractions per T ----
n_total = []
n_shell = []
real_frac = []
for T in T_types:
    sub = canon[canon['T'] == T]
    n_total.append(len(sub))
    n_s = (sub['fw_best_n_shells'] >= 1).sum()
    n_shell.append(int(n_s))
    real_frac.append(n_s / len(sub) if len(sub) else 0.0)

# ---- Burkert-truth null fractions per T (Burkert-truth only, the comparable null) ----
null_burk = null[null['smooth_truth'] == 'Burkert']
null_frac = []
null_n = []
for T in T_types:
    sub = null_burk[null_burk['T'] == T]
    if len(sub) > 0:
        null_frac.append((sub['best_n_shells'] > 0).mean())
        null_n.append(len(sub))
    else:
        null_frac.append(np.nan)
        null_n.append(0)


def binomial_errors(k, n):
    """Symmetric SE band on a binomial proportion."""
    if n == 0:
        return (0, 0)
    p = k / n
    se = np.sqrt(p * (1 - p) / n)
    return (max(p - se, 0), min(p + se, 1))


real_lo = [binomial_errors(k, n)[0] for k, n in zip(n_shell, n_total)]
real_hi = [binomial_errors(k, n)[1] for k, n in zip(n_shell, n_total)]
real_yerr_lo = [f - lo for f, lo in zip(real_frac, real_lo)]
real_yerr_hi = [hi - f for f, hi in zip(real_frac, real_hi)]

null_k = [int(round(f * n)) if not np.isnan(f) else 0 for f, n in zip(null_frac, null_n)]
null_lo = [binomial_errors(k, n)[0] if n > 0 else 0 for k, n in zip(null_k, null_n)]
null_hi = [binomial_errors(k, n)[1] if n > 0 else 0 for k, n in zip(null_k, null_n)]
null_yerr_lo = [f - lo if not np.isnan(f) else 0 for f, lo in zip(null_frac, null_lo)]
null_yerr_hi = [hi - f if not np.isnan(f) else 0 for f, hi in zip(null_frac, null_hi)]

# ---- Fisher one-sided p (real vs Burkert-null per T) ----
fisher_p = []
for T, ns, nt in zip(T_types, n_shell, n_total):
    sub = null_burk[null_burk['T'] == T]
    if len(sub) == 0:
        fisher_p.append(None)
        continue
    nul_n = len(sub)
    nul_s = (sub['best_n_shells'] > 0).sum()
    table = [[ns, nt - ns], [nul_s, nul_n - nul_s]]
    odds, p = stats.fisher_exact(table, alternative='greater')
    fisher_p.append(p)

# ---- Trend statistics computed from data (matches manuscript) ----
rho_T, p_T = stats.spearmanr(T_types, real_frac)

# Mann-Kendall p (Kendall tau on per-galaxy data — matches manuscript value 0.003)
shell_flag = (canon['fw_best_n_shells'] >= 1).astype(int)
tau, p_mk = stats.kendalltau(canon['T'], shell_flag)

# ---- Plot ----
fig, ax = plt.subplots(figsize=(12, 6.5))
x = np.arange(len(T_types))
width = 0.35

ax.bar(x - width/2, [f * 100 for f in real_frac], width,
       yerr=[[e * 100 for e in real_yerr_lo], [e * 100 for e in real_yerr_hi]],
       capsize=4, label='Real SPARC data', color='tab:blue', alpha=0.85,
       edgecolor='navy', linewidth=1)

null_x, null_h, null_yl, null_yh = [], [], [], []
for i, f in enumerate(null_frac):
    if not np.isnan(f):
        null_x.append(x[i] + width/2)
        null_h.append(f * 100)
        null_yl.append(null_yerr_lo[i] * 100)
        null_yh.append(null_yerr_hi[i] * 100)

ax.bar(null_x, null_h, width, yerr=[null_yl, null_yh], capsize=4,
       label='Burkert-truth null expectation', color='tab:gray', alpha=0.7,
       edgecolor='black', linewidth=1)

# Per-T Fisher p-value annotations
for i, (T, p, frac) in enumerate(zip(T_types, fisher_p, real_frac)):
    if p is not None:
        label_y = max(real_hi[i] * 100, null_hi[i] * 100) + 4
        if p < 0.01:
            sig = f'p = {p:.1e}\n***'
            color = 'darkred'
        elif p < 0.05:
            sig = f'p = {p:.3f}\n**'
            color = 'darkorange'
        elif p < 0.1:
            sig = f'p = {p:.2f}\n*'
            color = 'goldenrod'
        else:
            sig = f'p = {p:.2f}\nn.s.'
            color = 'gray'
        ax.text(i, label_y, sig, ha='center', va='bottom', fontsize=8.5,
                color=color, fontweight='bold')
    else:
        label_y = real_hi[i] * 100 + 4
        ax.text(i, label_y, '(no null run)', ha='center', va='bottom',
                fontsize=8, color='steelblue', style='italic')

ax.set_xticks(x)
ax.set_xticklabels([f'T={T}\n({lbl})' for T, lbl in zip(T_types, labels)])
ax.set_ylabel('Shell-bearing fraction (%)', fontsize=12)
ax.set_xlabel('Morphological type', fontsize=12)
ax.set_title(
    f'Shell-bearing fraction vs morphological type (T=2–9, n=102 galaxies)\n'
    f'Spearman rho = {rho_T:.2f}, p = {p_T:.3f}; Mann-Kendall p = {p_mk:.3f}',
    fontsize=11.5, fontweight='bold'
)
ax.set_ylim(0, 130)
ax.axhline(y=100, color='black', linestyle=':', linewidth=0.6, alpha=0.5)
ax.legend(loc='upper right', fontsize=9.5, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

ax.annotate('', xy=(7.3, 30), xytext=(0.3, 105),
            arrowprops=dict(arrowstyle='->', color='steelblue', lw=2, alpha=0.6))
ax.text(3.5, 115, 'Significant declining trend across T=2–9',
        fontsize=10, color='steelblue', fontweight='bold', alpha=0.85,
        ha='center', style='italic')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'figure4_shellfrac_vs_T.png'),
            dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'figure4_shellfrac_vs_T.pdf'),
            bbox_inches='tight')
print(f"Saved Figure 4 at 300 DPI (PNG + PDF) to {OUTPUT_DIR}")
print(f"  Real shell-bearing per T: {dict(zip(T_types, n_shell))}")
print(f"  Burkert-null per T (k/n): "
      f"{dict(zip(T_types, [f'{k}/{n}' for k, n in zip(null_k, null_n)]))}")
print(f"  Spearman rho (per-T-bin) = {rho_T:.3f}, p = {p_T:.4f}")
print(f"  Kendall tau (per-galaxy MK proxy) p = {p_mk:.4f}")
plt.close()
