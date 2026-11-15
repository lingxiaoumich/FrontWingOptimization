"""Compare finished optimizations: progress, the drag and ground-gap sensitivities, the start shape.

Reads each run's summary/chain_summary.json (plot_run_chain.py) and the surface of its best
iterate. Writes <out>/campaign.png|svg and <out>/campaign.json (one row per run).
  python scripts/plot_campaign.py --out optimizations/comparisons/sensitivity_20261115
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

O = Path(__file__).resolve().parents[1]/'optimizations'
RUNS = {  # run: (label, colour)
    'run_002_reducedBox': ('run_002: 10 pairs', '#9a9a9a'),
    'run_003_ffd20Pairs': ('run_003: 20 pairs (base)', '#2b2b2b'),
    'run_005_dragCap005': ('run_005: C_D c ≤ 0.05', '#2a78d6'),
    'run_008_dragCap003': ('run_008: C_D c ≤ 0.03', '#86b6ef'),
    'run_009_dragCap008': ('run_009: C_D c ≤ 0.08', '#0d366b'),
    'run_007_floor30mm': ('run_007: floor 30 mm', '#eb6834'),
    'run_011_naca0010LevelFront': ('run_011: NACA 0010 at 0°', '#1baf7a'),
    'run_006_levelNoTilt': ('run_006: level, vertical-only, tilt locked', '#8e44ad'),
}
p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)

rows, data = [], {}
for run, (label, col) in RUNS.items():
    f = O/run/'summary/chain_summary.json'
    if not f.exists(): continue
    s = json.loads(f.read_text()); its = s['iterates']
    best = min(its, key=lambda q: q['cl_c'])
    names = [r['run'] for r in s['runs']]
    it = O/run/names[best['run']]/'iterates'/f"eval_{best['evaluation']:03d}"
    settings = json.loads((O/run/names[0]/'settings.json').read_text())
    data[run] = dict(label=label, col=col, its=its, best=best, xy=np.loadtxt(it/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2],
                     box=settings['box_m'], cap=settings['drag_cap']*settings['reference_chord_m'], ranks=settings['mpi_ranks'])
    rows.append(dict(run=run, best_cl_c=best['cl_c'], best_cd_c=best['cd_c'], at_min=best['t_min'], grids=len(s['runs']),
                     drag_cap_cd_c=data[run]['cap'], floor_m=settings['box_m']['ymin'], ranks=settings['mpi_ranks'],
                     incidence_deg=best['incidence_deg'], best_in=f"{names[best['run']]} evaluation {best['evaluation']}"))
(out/'campaign.json').write_text(json.dumps(rows, indent=1)+'\n')

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c = '#2b2b2b', '#6b6b6b', '#e4e4e4'
fig, axs = plt.subplots(3, 2, figsize=(15, 14), layout='constrained')

ax = axs[0, 0]
for run, d in data.items():
    t = [q['t_min'] for q in d['its']]; c = [-q['cl_c'] for q in d['its']]
    ax.plot(t, c, c=d['col'], lw=1.6, label=f"{d['label']} ({d['ranks']} ranks)")
    ax.plot(d['best']['t_min'], -d['best']['cl_c'], 'o', ms=6, mfc='white', mec=d['col'], mew=1.5)
ax.set(xlabel='Wall time, cumulative (min); 24-rank runs take about twice as long per solve', ylabel='Downforce $-C_Lc$')
ax.set_title('(a) Progress of the PSQP iterates (best marked)', loc='left'); ax.legend(fontsize=8, frameon=False, loc='lower right')

ax = axs[0, 1]
for run, d in data.items():
    ax.plot(d['best']['cd_c'], -d['best']['cl_c'], 'o', ms=9, c=d['col'], label=d['label'])
    if d['cap'] < 1: ax.axvline(d['cap'], c=d['col'], lw=.8, ls=(0, (3, 2)))
caps = sorted((d['cap'], -d['best']['cl_c']) for d in data.values() if d['cap'] < 1 and d['box']['ymin'] == .05)
if len(caps) > 1: ax.plot(*zip(*caps), c='#2a78d6', lw=1, alpha=.5)
ax.set(xlabel='$C_Dc$ of the best design (dashed: its drag limit)', ylabel='Downforce $-C_Lc$', xlim=(0, None))
ax.set_title('(b) Downforce against drag: the drag-limit sensitivity', loc='left'); ax.legend(fontsize=8, frameon=False, loc='lower right')

def sections(ax, runs, title):
    b = data[runs[0]]['box']
    for run in runs:
        d = data[run]; bx = d['box']
        ax.add_patch(Rectangle((bx['xmin'], bx['ymin']), bx['xmax']-bx['xmin'], bx['ymax']-bx['ymin'], fill=False, ec=d['col'], lw=.8, ls=(0, (4, 2))))
        ax.plot(*np.r_[d['xy'], d['xy'][:1]].T, c=d['col'], lw=1.6,
                label=f"{d['label']}: $-C_Lc$ {-d['best']['cl_c']:.3f}, $C_Dc$ {d['best']['cd_c']:.3f}, tilt {-d['best']['incidence_deg']:.1f}°")
    ax.axhline(0, c=ink, lw=2); ax.set_aspect('equal'); ax.set(xlim=(-.03, .53), ylim=(-.005, .2), xlabel='x (m)', ylabel='height (m)')
    ax.legend(fontsize=8, frameon=False, loc='upper left'); ax.set_title(title, loc='left')
sections(axs[1, 0], [r for r in ['run_003_ffd20Pairs', 'run_007_floor30mm'] if r in data], '(c) Ground gap: best sections, floor 50 mm and 30 mm (dashed boxes)')
sections(axs[1, 1], [r for r in ['run_003_ffd20Pairs', 'run_011_naca0010LevelFront'] if r in data], '(d) Start shape: best sections from a NACA 0012 at 3° and a NACA 0010 at 0°')
sections(axs[2, 0], [r for r in ['run_003_ffd20Pairs', 'run_006_levelNoTilt'] if r in data], '(e) Control points moving in x and y vs only vertically (tilt locked)')
ax = axs[2, 1]; ax.axis('off'); ax.set_title('(f) Best design of each run', loc='left')
ax.table(cellText=[[d['label'], f"{d['best']['cl_c']:.3f}", f"{d['best']['cd_c']:.3f}", f"{d['ranks']}"] for d in data.values()],
         colLabels=['run', '$C_Lc$', '$C_Dc$', 'ranks'], loc='upper left', cellLoc='left', colWidths=[.62, .13, .13, .1]).scale(1, 1.45)
for ax in axs.flat[:5]: ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
fig.savefig(out/'campaign.png', dpi=150); fig.savefig(out/'campaign.svg'); plt.close(fig)
for r in rows: print(f"{r['run']:28s} C_L c {r['best_cl_c']:7.3f}  C_D c {r['best_cd_c']:.3f}  at {r['at_min']:5.0f} min on {r['ranks']} ranks, {r['grids']} grids, floor {r['floor_m']*1e3:.0f} mm, cap {r['drag_cap_cd_c']:g}")
