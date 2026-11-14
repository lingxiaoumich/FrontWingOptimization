"""Compare optimization chains that differ in their FFD parameterization.

For the November 1 runs from the NACA 0012 baseline to the corner of the box with 5 and
with 10 upper/lower FFD pairs. Each chain is the --out folder of plot_run_chain.py (its
chain_summary.json) and the surface of its best design. Writes under <out>/plots:
  pairs_comparison  downforce over run time; drag polar; the best sections, true scale
  python scripts/plot_pairs_comparison.py --out optimizations/ffd10_summary_20261101 \\
      --chain '5 pairs' optimizations/box_summary_20261101 optimizations/box_remesh2_xmax1.0_gap2mm_20261101/best_feasible/surface_geometry.csv \\
      --chain '10 pairs' optimizations/ffd10_summary_20261101 optimizations/ffd10_remesh3b_xmax1.0_gap2mm_20261101/best_feasible/surface_geometry.csv
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
p.add_argument('--chain', nargs=3, action='append', required=True, metavar=('LABEL', 'SUMMARY_DIR', 'BEST_SURFACE_CSV'))
a = p.parse_args()
out = a.out.resolve(); (out/'plots').mkdir(parents=True, exist_ok=True)

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c = '#2b2b2b', '#6b6b6b', '#e4e4e4'
COLORS = ['#2a78d6', '#eb6834']
box = dict(xmin=-.125, xmax=.625, ymin=.05, ymax=.35)
LIMITS = ['leading edge to the box', 'trailing edge to the box', 'clearance above the box floor',
          'top of the section to the box top', 'gap above the 2 mm floor']

chains = []
for label, folder, surface in a.chain:
    s = json.loads((Path(folder)/'chain_summary.json').read_text())
    its = s['iterates']
    # The corner: every box limit and the gap floor within 1 mm and 0.1 mm.
    corner = next((q for q in its if all(q[k] <= 1. for k in LIMITS[:4]) and q[LIMITS[4]] <= .1), None)
    best = min(its, key=lambda q: q['cl_c'])
    chains.append(dict(label=label, its=its, corner=corner, best=best,
                       surf=np.loadtxt(surface, delimiter=',', skiprows=1)))

fig, axs = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 1.25]))
ax = axs[0]
for c, col in zip(chains, COLORS):
    t = [q['t_min'] for q in c['its']]; d = [-q['cl_c'] for q in c['its']]
    ax.plot(t, d, c=col, lw=1.6, marker='o', ms=2.8, mec='white', mew=.3, label=c['label'])
    if c['corner']:
        q = c['corner']
        ax.plot(q['t_min'], -q['cl_c'], 'o', ms=9, mfc='none', mec=col, mew=1.4)
        ax.annotate(f"corner reached,\n{q['t_min']:.0f} min", (q['t_min'], -q['cl_c']), xytext=(9, -26), textcoords='offset points',
                    ha='left', fontsize=8, color=ink, arrowprops=dict(arrowstyle='-', color=muted, lw=.6))
    b = c['best']
    ax.annotate(f"{c['label']}: {-b['cl_c']:.3f}", (b['t_min'], -b['cl_c']), xytext=(4, 4 if col == COLORS[1] else -12),
                textcoords='offset points', fontsize=8.5, color=ink)
ax.set(xlabel='Run time, cumulative (min; hand-overs between runs excluded)', ylabel='Downforce $-C_Lc$',
       title='(a) Progress: PSQP iterates')
ax.legend(fontsize=8.5, frameon=False, loc='lower right')

ax = axs[1]
for c, col in zip(chains, COLORS):
    ax.plot([-q['cl_c'] for q in c['its']], [q['cd_c'] for q in c['its']], c=col, lw=1.6, marker='o', ms=2.8,
            mec='white', mew=.3, label=c['label'])
    b = c['best']
    ax.annotate(f"{b['cd_c']:.3f}", (-b['cl_c'], b['cd_c']), xytext=(6, 0 if col == COLORS[1] else -10),
                textcoords='offset points', ha='left', va='center', fontsize=8.5, color=ink)
ax.set(xlabel='Downforce $-C_Lc$', ylabel='Drag $C_Dc$', title='(b) Drag polar (cap $C_Dc$ = 5, never near)')
ax.legend(fontsize=8.5, frameon=False, loc='upper left')

ax = axs[2]
ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=1.1))
ax.text(box['xmin'], box['ymax']+.004, 'Allowed surface box', ha='left', va='bottom', color=muted, fontsize=9)
ax.axhline(0, c=ink, lw=2); ax.text(box['xmax'], .004, 'Slip ground', ha='right', va='bottom', color=muted, fontsize=9)
for c, col, ls in zip(chains, COLORS, ['-', (0, (5, 2))]):
    s = c['surf']; b = c['best']
    ax.plot(np.r_[s[:, 0], s[0, 0]], np.r_[s[:, 1], s[0, 1]], c=col, lw=1.5, ls=ls,
            label=f"{c['label']}: $-C_Lc$ {-b['cl_c']:.3f}, $C_Dc$ {b['cd_c']:.3f}, incidence {b['incidence_deg']:.1f}°")
ax.set_aspect('equal'); ax.set_xlim(-.16, .66); ax.set_ylim(-.01, .40)
ax.set(xlabel='x (m)', ylabel='y (m)', title='(c) Best sections, true scale')
ax.legend(fontsize=8, frameon=False, loc='upper left', bbox_to_anchor=(.07, .84))  # inside the box, above the sections
for ax in axs:
    ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True); ax.set_title(ax.get_title(), fontsize=10.5)
fig.tight_layout()
fig.savefig(out/'plots/pairs_comparison.png', dpi=170, bbox_inches='tight')
fig.savefig(out/'plots/pairs_comparison.svg', bbox_inches='tight'); plt.close(fig)
print(json.dumps({c['label']: {'corner': c['corner'] and {k: c['corner'][k] for k in ['k', 't_min', 'cl_c', 'cd_c']},
                               'best': {k: c['best'][k] for k in ['k', 'run', 'evaluation', 't_min', 'cl_c', 'cd_c', 'incidence_deg']}}
                  for c in chains}, indent=1))
