"""Animate a chain of runs, one frame per PSQP iterate, in four panels:
  top left      total-pressure coefficient Cpt (as plot_iterate_flows.py), u < 0 outlined
  top right     pressure coefficient Cp: maize at +1, white at 0, blue at -10
  bottom left   the FFD control points and the section; their trails fade with age
  bottom right  C_L c against C_D c
The iterates are those of <out>/chain_summary.json (plot_run_chain.py); each needs its flow
file. Writes <out>/plots/chain_iterates.gif. Run with the MACH-Aero env:
  python scripts/plot_chain_gif.py --out optimizations/run_NNN_<name>/summary optimizations/run_NNN_<name>/r00 ...
"""
import argparse
import json
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm, to_rgba
from matplotlib.patches import Rectangle
from PIL import Image

from flow_fields import add_cpt, read_zones, tri

p = argparse.ArgumentParser()
p.add_argument('runs', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True)
p.add_argument('--ms', type=int, default=250, help='ms per frame; the last one stays 3 s')
p.add_argument('--jobs', type=int, default=8, help='frames drawn at once (keep cores free for a running solver)')
a = p.parse_args()
out = a.out.resolve(); runs = [r.resolve() for r in a.runs]
settings = json.loads((runs[0]/'settings.json').read_text())
mach, box, c0 = settings['flow_derived']['mach'], settings['box_m'], settings['reference_chord_m']

def lattice(path):
    lines = path.read_text().split('\n'); dims = [int(v) for v in lines[1].split()]
    vals = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
    return np.stack([v.reshape(dims, order='F') for v in vals], -1)  # (station, lower/upper, span, xyz)
# The control points at a design, as remesh_design.py: rotate the run's lattice about the
# incidence pivot (the original one in every run), then add the shape moves.
xy0 = np.loadtxt(runs[0]/'source/airfoil_coordinates.csv', delimiter=',', skiprows=1)
pivot = np.array([.25*c0, xy0[136, 1]+.25*c0*np.sin(np.deg2rad(3.))])
def control_points(L0, design):
    t = np.deg2rad(-(design['incidence'][0]+3.)); R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    P = (L0[:, :, 0, :2]-pivot)@R.T+pivot
    dx = design.get('shape_x', np.zeros(P[..., 0].size))  # --shape-directions y: only 'shape', the y moves
    dy = design['shape_y'] if 'shape_y' in design else design['shape']
    return P+np.stack([np.reshape(dx, P.shape[:2]), np.reshape(dy, P.shape[:2])], -1)

its = json.loads((out/'chain_summary.json').read_text())['iterates']
L0 = [lattice(r/'source/ffd.xyz') for r in runs]
for q in its:
    q['dir'] = runs[q['run']]/'iterates'/f"eval_{q['evaluation']:03d}"
    q['P'] = control_points(L0[q['run']], json.loads((q['dir']/'design.json').read_text()))
    q['xy'] = np.loadtxt(q['dir']/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]
N = len(its)
P = np.array([q['P'] for q in its]); t = np.array([q['t_min'] for q in its])
cl, cd, g = (np.array([q[k] for q in its]) for k in ('cl_c', 'cd_c', 'run'))

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c, orange = '#2b2b2b', '#6b6b6b', '#e4e4e4', '#eb6834'
COLORS = ['#2a78d6', '#1baf7a', '#0d366b', '#86b6ef']  # as plot_run_chain.py: original grid, then each remesh
grid_name = lambda k: 'original grid' if k == 0 else f'remesh {k}'
CPT = dict(cmap=LinearSegmentedColormap.from_list('cpt', ['#0d366b', '#2a78d6', '#86b6ef', '#f4f8fd']),
           norm=Normalize(-1.2, 1.02), levels=np.linspace(-1.2, 1.02, 75))
# Michigan blue at Cp = -10 (through a mid blue, not grey), white at 0, maize at +1
CP = dict(cmap=LinearSegmentedColormap.from_list('cp', [(0, '#00274c'), (.3, '#2f6fbf'), (.5, 'white'), (1, '#ffcb05')]),
          norm=TwoSlopeNorm(0., -10., 1.), levels=np.r_[np.linspace(-10, 0, 51), np.linspace(0, 1, 11)[1:]],
          ticks=[-10, -8, -6, -4, -2, 0, .5, 1])
CPT['ticks'] = None

def scene(ax, xy):
    ax.fill(*xy.T, color='#9a9a9a', lw=.6, ec=ink, zorder=3)
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=ink, lw=.8, zorder=4))
    ax.axhline(0, c=ink, lw=2, zorder=4); ax.set_aspect('equal')

def field(ax, zones, key, style, xy, label):
    for tag in ['background', 'near']:
        f = zones[tag]
        ax.tricontourf(tri(f, tag == 'near'), np.clip(f[key], -50, 2).ravel(), extend='both', zorder=2, rasterized=True,
                       cmap=style['cmap'], norm=style['norm'], levels=style['levels'])
    scene(ax, xy); ax.set(xlim=(-.2, 1.05), ylim=(0, .42), xlabel='x (m)', ylabel='height (m)')
    plt.colorbar(ScalarMappable(style['norm'], style['cmap']), ax=ax, shrink=.85, pad=.02, label=label, extend='both', ticks=style['ticks'])

def frame(n):
    q = its[n]
    zones = read_zones(q['dir']/'flow_000_vol.cgns'); add_cpt(zones, mach)
    fig, axs = plt.subplots(2, 2, figsize=(14, 8.4), dpi=90, layout='constrained', gridspec_kw=dict(height_ratios=[1, 1.45]))
    ax = axs[0, 0]; field(ax, zones, 'Cpt', CPT, q['xy'], '$C_{pt}$'); ax.set_title('(a) Total pressure; u < 0 outlined', loc='left')
    for tag in ['background', 'near']:
        f = zones[tag]; ax.tricontour(tri(f, tag == 'near'), f['VelocityX'].ravel(), levels=[0], colors=orange, linewidths=1.1, zorder=3)
    ax = axs[0, 1]; field(ax, zones, 'CoefPressure', CP, q['xy'], '$C_p$'); ax.set_title('(b) Pressure', loc='left')

    ax = axs[1, 0]
    ax.plot(*np.r_[its[0]['xy'], its[0]['xy'][:1]].T, c=muted, lw=.8, ls=(0, (3, 2)), zorder=2, label='starting section')
    scene(ax, q['xy'])
    seg = np.stack([P[:n], P[1:n+1]], -2).transpose(1, 2, 0, 3, 4).reshape(-1, 2, 2) if n else np.empty((0, 2, 2))  # (point, step, 2, xy)
    age = np.tile(np.arange(n)[::-1], P.shape[1]*P.shape[2])
    rgba = np.array([to_rgba(COLORS[g[i+1] % 4]) for i in np.tile(np.arange(n), P.shape[1]*P.shape[2])]) if n else np.empty((0, 4))
    if n: rgba[:, 3] = np.maximum(.05, .9*.95**age)  # older steps fade
    ax.add_collection(LineCollection(seg, colors=rgba, linewidths=1.3, zorder=5))
    C = q['P']
    for j in range(2): ax.plot(*C[:, j].T, c=ink, lw=.7, alpha=.5, zorder=6)
    ax.plot(C[:, :, 0].T, C[:, :, 1].T, c=ink, lw=.7, alpha=.5, zorder=6)
    ax.plot(*C.reshape(-1, 2).T, 'o', ms=5, mfc='white', mec=ink, mew=1.2, zorder=7, label=f'control points ({C.shape[0]} upper/lower pairs)')
    ax.set(xlim=(-.1, .63), ylim=(-.03, .45), xlabel='x (m)', ylabel='height (m)')
    ax.legend(loc='upper left', fontsize=8.5, frameon=False); ax.set_title('(c) FFD control points and their paths (fading with age)', loc='left')

    ax = axs[1, 1]
    for k in range(len(runs)):
        s = g == k; done = s & (np.arange(N) <= n)
        ax.plot(cd[s], cl[s], c=COLORS[k % 4], lw=1.2, alpha=.25)
        ax.plot(cd[done], cl[done], c=COLORS[k % 4], lw=1.8, marker='o', ms=2.5, label=grid_name(k))
    ax.plot(cd[n], cl[n], 'o', ms=9, mfc='white', mec=ink, mew=1.5, zorder=5)
    ax.set(xlabel='$C_Dc$', ylabel='$C_Lc$', xlim=(0, cd.max()*1.08), ylim=(cl.min()-.08*np.ptp(cl), cl.max()+.08*np.ptp(cl)))
    ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True); ax.legend(loc='upper right', fontsize=8.5, frameon=False)
    ax.set_title('(d) $C_Lc$ against $C_Dc$ (more downforce is lower)', loc='left')
    fig.suptitle(f"{out.parent.name} · PSQP iteration {n} of {N-1} · {runs[q['run']].name} ({grid_name(q['run'])}) · {t[n]:.0f} min · "
                 f"$C_Lc$ {cl[n]:.3f}, $C_Dc$ {cd[n]:.3f}", fontsize=11.5)
    path = Path(tmp)/f'{n:04d}.png'; fig.savefig(path); plt.close(fig)
    return path

with tempfile.TemporaryDirectory() as tmp, Pool(a.jobs) as pool:
    frames = [Image.open(f).convert('RGB') for f in pool.map(frame, range(N))]
(out/'plots').mkdir(exist_ok=True)
frames[0].save(out/'plots/chain_iterates.gif', save_all=True, append_images=frames[1:],
               duration=[a.ms]*(N-1)+[3000], loop=0, optimize=True)
print(out/'plots/chain_iterates.gif', N, 'frames', frames[0].size)
