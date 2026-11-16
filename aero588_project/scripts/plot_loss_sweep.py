"""Outlet Cpt loss sweep on the run_006 baseline: downforce against the loss, against drag (next to the
drag-cap sweep on run_003), and the best sections.

Reads each run's best design inside every limit (<run>/r*/best_feasible, lowest C_L c over its grids);
run_006's loss comes from run_006_levelNoTilt/cptloss_eval. Writes <out>/loss_sweep.png|svg|json.
  python scripts/plot_loss_sweep.py --out optimizations/comparisons/cptloss_20261116
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
RUNS = {  # run: (label, colour); one blue, darker = looser cap
    'run_006_levelNoTilt': ('run_006 (no cap)', '#0d366b'),
    'run_012_cptLoss075': ('run_012: 75 %', '#1f5fae'),
    'run_013_cptLoss050': ('run_013: 50 %', '#2a78d6'),
    'run_014_cptLoss025': ('run_014: 25 %', '#5d9ae6'),
    'run_015_cptLoss0125': ('run_015: 12.5 %', '#93bdf0'),
}
DRAG_SWEEP = ['run_008_dragCap003', 'run_005_dragCap005', 'run_009_dragCap008', 'run_003_ffd20Pairs']
p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)


def best_feasible(run):
    best = None
    for r in sorted((O/run).glob('r[0-9]*')):
        f = r/'best_feasible/result.json'
        if f.exists():
            d = json.loads(f.read_text())
            if best is None or d['cl'] < best[1]['cl']: best = (r, d)
    return best


rows, data = [], {}
for run, (label, col) in RUNS.items():
    r, d = best_feasible(run)
    s = json.loads((O/run/'r00/settings.json').read_text())
    if 'cpt_loss_outlet' not in d:
        e = json.loads((O/run/'cptloss_eval/best_feasible/result.json').read_text())
        d = {**d, **{k: e[k] for k in ('cpt_loss_outlet', 'cpt_loss_outlet_flux_c_m')}}
    xy = np.loadtxt(r/'best_feasible/surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]
    data[run] = dict(label=label, col=col, cl_c=d['cl']/2, cd_c=d['cd']/2, loss=d['cpt_loss_outlet'],
                     flux=d['cpt_loss_outlet_flux_c_m'], cap=s.get('loss_cap'), box=s['box_m'], xy=xy)
    rows.append(dict(run=run, best_in=f"{r.name} evaluation {d['evaluation']}", cl_c=d['cl']/2, cd_c=d['cd']/2,
                     cpt_loss_outlet=d['cpt_loss_outlet'], loss_cap=s.get('loss_cap'),
                     loss_flux_over_cd_c=d['cpt_loss_outlet_flux_c_m']/(d['cd']/2)))
camp = {q['run']: q for q in json.loads((O/'comparisons/sensitivity_20261115/campaign.json').read_text())}
(out/'loss_sweep.json').write_text(json.dumps(rows, indent=1)+'\n')

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
muted, grid_c, orange = '#6b6b6b', '#e4e4e4', '#eb6834'
fig = plt.figure(figsize=(14, 9), layout='constrained')
gs = fig.add_gridspec(2, 2, height_ratios=[1, .8])
base = data['run_006_levelNoTilt']

ax = fig.add_subplot(gs[0, 0])
L = [d['loss'] for d in data.values()]; C = [-d['cl_c'] for d in data.values()]
ax.plot(L, C, c=muted, lw=1, zorder=1)
for run, d in data.items():
    ax.scatter(d['loss'], -d['cl_c'], s=70, c=d['col'], edgecolors='white', linewidths=2, zorder=3)
    if d['cap']: ax.axvline(d['cap'], c=d['col'], lw=.8, ls=(0, (3, 2)), zorder=0)
    ax.annotate(f"{d['label']}\n$-C_Lc$ {-d['cl_c']:.3f}", (d['loss'], -d['cl_c']), xytext=(8, -4),
                textcoords='offset points', fontsize=8.5, va='top')
ax.set(xlabel='Outlet $C_{pt}$ loss of the best design (dashed: caps)', ylabel='Downforce $-C_Lc$', xlim=(0, .026), ylim=(1, 2.45))
ax.set_title('(a) Downforce against the outlet loss', loc='left'); ax.grid(c=grid_c, lw=.6)

ax = fig.add_subplot(gs[0, 1])
b3 = camp['run_003_ffd20Pairs']
dx = [camp[r]['best_cd_c'] for r in DRAG_SWEEP]; dy = [100*camp[r]['best_cl_c']/b3['best_cl_c'] for r in DRAG_SWEEP]
ax.plot(dx, dy, 'o-', c=orange, lw=1.4, ms=7, mec='white', mew=1.5, label='drag cap, base run_003 (x + y, tilt free)')
lx = [d['cd_c'] for d in data.values()][::-1]; ly = [100*d['cl_c']/base['cl_c'] for d in data.values()][::-1]
ax.plot(lx, ly, 'o-', c='#2a78d6', lw=1.4, ms=7, mec='white', mew=1.5, label='outlet loss cap, base run_006 (y only, tilt locked)')
ax.set(xlabel='$C_Dc$ of the best design', ylabel='Downforce, % of its own base', xlim=(0, .18), ylim=(45, 103))
ax.set_title('(b) Same trade per unit drag: loss cap vs drag cap', loc='left'); ax.grid(c=grid_c, lw=.6)
ax.legend(fontsize=8.5, frameon=False, loc='lower right')

ax = fig.add_subplot(gs[1, :])
bx = base['box']
ax.add_patch(Rectangle((bx['xmin'], bx['ymin']), bx['xmax']-bx['xmin'], bx['ymax']-bx['ymin'], fill=False, ec=muted, lw=.8, ls=(0, (4, 3))))
for run, d in data.items():
    ax.fill(*d['xy'].T, fc=d['col'], alpha=.12, ec='none'); ax.plot(*d['xy'].T, c=d['col'], lw=1.4,
            label=f"{d['label']}: $C_Lc$ {d['cl_c']:.3f}, $C_Dc$ {d['cd_c']:.3f}, loss {d['loss']:.4f}")
ax.set(xlabel='x (m)', ylabel='height (m)', xlim=(-.03, .56), ylim=(.035, .2)); ax.set_aspect('equal')
ax.set_title('(c) Best sections (dashed: the design box; its floor is the 50 mm ground clearance)', loc='left'); ax.legend(fontsize=8.5, frameon=False, loc='upper left', ncol=2)

fig.suptitle('Outlet $C_{pt}$ loss sweep on run_006 (caps as a share of run_006\'s 0.0221)', x=.01, ha='left', fontsize=12)
for ext in ('png', 'svg'): fig.savefig(out/f'loss_sweep.{ext}', dpi=110)
for r in rows: print(f"{r['run']:22s} C_L c {r['cl_c']:7.3f}  C_D c {r['cd_c']:.4f}  loss {r['cpt_loss_outlet']:.5f}  cap {r['loss_cap']}  flux/C_D c {r['loss_flux_over_cd_c']:.2f}  ({r['best_in']})")
