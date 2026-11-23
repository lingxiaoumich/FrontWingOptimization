"""Figures for the campaign report (reports/report_20261122.html).

Per-run figures go to optimizations/<run>/summary/plots/:
  best_section.png      the best design inside every limit, on a fixed background grid
  setup_geometry.png    run_016 only: start section, FFD control points, box, ground
  mesh.png              run_016 only: the overset grid, whole domain and near the wing
  result_overview.png   run_016 only: start and best sections, surface Cp, trailing edge, progress
Cross-run figures go to optimizations/comparisons/report_20261122/:
  drag_sweep, loss_sweep, floor_sweep, cpmin_sweep, misc_sections (.png), and contact_<group>.png
Numbers come from <scratch>/runs.json (gather.py). Run with the MACH-Aero env from scripts/:
  python plot_report.py <runs.json>
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

from flow_fields import read_zones, tri
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

O = Path(__file__).resolve().parents[1]/'optimizations'
CMP = O/'comparisons/report_20261122'; CMP.mkdir(parents=True, exist_ok=True)
RUNS = json.loads(Path(sys.argv[1]).read_text())
C0 = 0.5
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': 'white',
                     'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans', 'axes.titlesize': 10.5, 'axes.titlelocation': 'left'})
ink, muted, grid_c, ground_c = '#2b2b2b', '#6b6b6b', '#e4e4e4', '#d4d4d4'
BLUES = ['#93bff2', '#4387dd', '#1f62bd', '#154a8f', '#0d3366']
ORANGE = '#eb6834'
BOX = RUNS['run_016_area050']['box']

def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig); print(path)

def best_dir(R):
    b = RUNS[R]['best']; return O/R/b['round']/'best_feasible'

def surface(d):
    return np.loadtxt(d/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]

def slices(d):
    """Surface x, y, Cp, y+ of a slice file, put in the loop order of the design's surface vertices (nearest vertex)."""
    lines = (d/'flow_000_slices.dat').read_text().split('\n')
    n = int(lines[3].split()[2])
    P = np.array([l.split() for l in lines[5:5+n]], float)
    xy = surface(d)
    near = np.argmin((P[:, 0, None]-xy[None, :, 0])**2+(P[:, 1, None]-xy[None, :, 1])**2, axis=1)
    P = P[np.argsort(near, kind='stable')]; near = np.sort(near)
    lo, up = halves(xy)
    side = np.isin(near, lo)
    return dict(x=P[:, 0], y=P[:, 1], cp=P[:, 6], yplus=P[:, 8], suction=side)

def halves(xy):
    """Index arrays of the suction (lower) and pressure (upper) sides of a closed loop that starts at the trailing edge."""
    le = int(np.argmin(xy[:, 0])); a, b = np.arange(0, le+1), np.arange(le, len(xy))
    return (a, b) if xy[a, 1].mean() < xy[b, 1].mean() else (b, a)

def scene(ax, box=BOX, free_air=False, xlim=(-.03, .56)):
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=.8, zorder=2))
    if free_air: ax.axhline(0, c=muted, lw=.8, ls=':')
    else:
        ax.add_patch(Rectangle((xlim[0]-1, -1), xlim[1]-xlim[0]+2, 1, color=ground_c, lw=0, zorder=1)); ax.axhline(0, c=ink, lw=1.4)
    ax.set_aspect('equal')

def shape_grid(ax):
    """The same background grid in every run, for comparing shapes: 50 mm in x, 25 mm in y."""
    ax.xaxis.set_major_locator(MultipleLocator(.1)); ax.xaxis.set_minor_locator(MultipleLocator(.05))
    ax.yaxis.set_major_locator(MultipleLocator(.05)); ax.yaxis.set_minor_locator(MultipleLocator(.025))
    ax.grid(which='both', color=grid_c, lw=.6); ax.set_axisbelow(True)

def draw_section(ax, R, title=None, xlim=(-.03, .56), ylim=(-.012, .20), numbers=True):
    d = best_dir(R); xy = surface(d); b = RUNS[R]['best']
    free = (RUNS[R]['ground'] or '').startswith('none')
    scene(ax, RUNS[R]['box'], free, xlim)
    ax.fill(*xy.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4)
    ax.set(xlim=xlim, ylim=ylim); shape_grid(ax); ax.tick_params(labelsize=8)
    nums = f"$-C_Lc$ {-b['clc']:.2f}, $C_Dc$ {b['cdc']:.3f}, lowest point {b['clearance_mm']:.0f} mm"
    if title: ax.set_title(f"{title}\n{nums}" if numbers else title, fontsize=9)

def contact(group, items, ncol=3, ylim=(-.012, .20)):
    n = len(items); nrow = int(np.ceil(n/ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.0*ncol, 2.4*nrow+.3), squeeze=False)
    for ax, (R, title) in zip(axs.ravel(), items): draw_section(ax, R, title, ylim=ylim)
    for ax in axs.ravel()[n:]: ax.axis('off')
    fig.text(.5, -.01, 'Best feasible design of each run. Dashed: regulation box. Grid: 50 mm in x, 25 mm in y.', ha='center', fontsize=8.5, color=muted)
    fig.tight_layout(); save(fig, CMP/f'contact_{group}.png')

def per_run_section(R):
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    b = RUNS[R]['best']
    draw_section(ax, R, f"{R}: best inside every limit ({b['round']}, evaluation {b['evaluation']})", ylim=(-.015, .22))
    ax.set(xlabel='x (m)', ylabel='y (m)')
    save(fig, O/R/'summary/plots/best_section.png')

def polar(ax, items, caps=(), label_dx=.003):
    """Best points: items = [(R, label, colour)]; caps = [(C_D c, colour)] dashed."""
    for c, col in caps: ax.axvline(c, c=col, lw=.9, ls='--', zorder=1)
    for R, lab, col, *off in items:
        b = RUNS[R]['best']; dx, dy, ha = off[0] if off else (6, -3, 'left')
        ax.plot(b['cdc'], -b['clc'], 'o', ms=8, mfc=col, mec=ink, mew=.8, zorder=3)
        ax.annotate(lab, (b['cdc'], -b['clc']), xytext=(dx, dy), textcoords='offset points', fontsize=8.5, va='top', ha=ha)
    ax.set(xlabel='Drag $C_Dc$', ylabel='Downforce $-C_Lc$'); ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)

def front(ax, Rs, col=muted, label=None):
    P = sorted((RUNS[R]['best']['cdc'], -RUNS[R]['best']['clc']) for R in Rs)
    ax.plot(*np.array(P).T, c=col, lw=1, zorder=2, label=label)

# ------------------------------------------------------------------ run_016: setup, mesh, result
R16 = 'run_016_area050'
start_d = O/R16/'r00/baseline'
xy0 = surface(start_d)
lat = (O/R16/'r00/source/ffd.xyz').read_text().split('\n'); dims = [int(v) for v in lat[1].split()]
vals = np.array(' '.join(lat[2:]).split(), float).reshape(3, -1)
LX, LY = (v.reshape(dims, order='F')[:, :, 0] for v in vals[:2])
fig, ax = plt.subplots(figsize=(9, 3.6))
scene(ax, xlim=(-.06, .58))
ax.fill(*xy0.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4)
for j in range(2): ax.plot(LX[:, j], LY[:, j], c=BLUES[2], lw=.8, zorder=5)
ax.plot(LX.T, LY.T, c=BLUES[2], lw=.6, alpha=.6, zorder=5)
ax.plot(LX.ravel(), LY.ravel(), 'o', ms=4.5, mfc='white', mec=BLUES[3], mew=1.1, zorder=6, label='FFD control points (20 pairs, move up/down)')
ax.plot([], [], ls='--', c=muted, lw=.8, label='regulation box: x −12.3…512 mm, y 50…188.1 mm')
ax.annotate('', (BOX['xmin'], .03), (BOX['xmax'], .03), arrowprops=dict(arrowstyle='<->', color=muted, lw=.8))
ax.text(.25, .033, 'chord 524.3 mm (box length)', ha='center', va='bottom', fontsize=8.5, color=muted)
ax.annotate('', (.535, 0), (.535, BOX['ymin']), arrowprops=dict(arrowstyle='<->', color=muted, lw=.8))
ax.text(.54, .025, '50 mm', fontsize=8.5, color=muted, va='center')
ax.text(.45, -.008, 'ground (slip wall)', fontsize=8.5, color=muted, va='top', ha='center')
ax.set(xlim=(-.06, .6), ylim=(-.02, .2), xlabel='x (m)', ylabel='y (m)')
ax.legend(loc='upper center', bbox_to_anchor=(.5, -.2), fontsize=8.5, frameon=False, ncol=2)
save(fig, O/R16/'summary/plots/setup_geometry.png')

z = read_zones(start_d/'flow_000_vol.cgns')
fig, axs = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw=dict(width_ratios=[1.35, 1]))
bx, by = z['background']['CoordinateX'], z['background']['CoordinateY']
nx, ny = z['near']['CoordinateX'], z['near']['CoordinateY']
ax = axs[0]
ax.plot(bx[::3, :].T, by[::3, :].T, c='#9a9a9a', lw=.3); ax.plot(bx[:, ::3], by[:, ::3], c='#9a9a9a', lw=.3)
ax.fill(nx[:, -1], ny[:, -1], c=BLUES[2], lw=0)
ax.axvline(7.3, c=ORANGE, lw=1.2); ax.text(7.25, 4.6, 'loss plane\nx = 7.3 m', ha='right', va='top', fontsize=8.5, color=ORANGE)
ax.set(aspect='equal', xlim=(-5.2, 7.7), ylim=(-.2, 5.2), xlabel='x (m)', ylabel='y (m)', title='(a) Whole domain (every 3rd background line)')
ax.text(-5.1, 2.5, 'far field', rotation=90, va='center', fontsize=8.5, color=muted); ax.text(1, 5.05, 'far field', fontsize=8.5, color=muted)
ax.text(1, .08, 'ground: slip wall', fontsize=8.5, color=muted)
ax = axs[1]
ax.plot(bx.T, by.T, c='#b4b4b4', lw=.3); ax.plot(bx, by, c='#b4b4b4', lw=.3)
ax.plot(nx[:, ::2], ny[:, ::2], c=BLUES[2], lw=.3); ax.plot(nx[::2, :].T, ny[::2, :].T, c=BLUES[2], lw=.3)
ax.fill(*xy0.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4); ax.axhline(0, c=ink, lw=1.4)
ax.set(aspect='equal', xlim=(-.06, .58), ylim=(-.005, .17), xlabel='x (m)', ylabel='y (m)', title='(b) Near the wing: O-grid (blue) over the background (grey)')
fig.tight_layout(); save(fig, O/R16/'summary/plots/mesh.png')

b16 = RUNS[R16]['best']; d16 = best_dir(R16); xyb = surface(d16)
s0, s1 = slices(start_d), slices(d16)
fig = plt.figure(figsize=(13, 7.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
ax = fig.add_subplot(gs[0, :])
scene(ax, xlim=(-.03, .56))
ax.plot(*np.r_[xy0, xy0[:1]].T, c=muted, lw=1, ls=(0, (3, 2)), label=f"start: level NACA 0012, $-C_Lc$ {-RUNS[R16]['start']['clc']:.2f}, $C_Dc$ {RUNS[R16]['start']['cdc']:.3f}")
ax.fill(*xyb.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4, label=f"best: $-C_Lc$ {-b16['clc']:.2f}, $C_Dc$ {b16['cdc']:.3f}, area {b16['area_pct']:.0f} % of start")
ax.set(xlim=(-.03, .56), ylim=(-.012, .2), xlabel='x (m)', ylabel='y (m)', title='(a) Start and best section')
ax.legend(loc='lower left', bbox_to_anchor=(.0, .07), fontsize=8.5, frameon=False)
ax = fig.add_subplot(gs[1, 0])
for s, lab, ls in ((s0, 'start', (0, (3, 2))), (s1, 'best', '-')):
    lo, up = s['suction'], ~s['suction']
    ax.plot(s['x'][lo]/C0, s['cp'][lo], c=BLUES[3] if lab == 'best' else muted, ls=ls, lw=1.4, label=f'{lab}, suction (lower) side')
    ax.plot(s['x'][up]/C0, s['cp'][up], c=ORANGE if lab == 'best' else '#b4b4b4', ls=ls, lw=1.2, label=f'{lab}, pressure (upper) side')
ax.invert_yaxis(); ax.set(xlabel='x/c', ylabel='$C_p$', title='(b) Surface pressure'); ax.grid(c=grid_c, lw=.6)
ax.legend(fontsize=8, frameon=False, loc='upper right')
ax = fig.add_subplot(gs[1, 1])
H = []
for r in RUNS[R16]['rounds']:
    for l in open(O/R16/r/'history.jsonl'): H.append(json.loads(l) | {'round': r})
cs = json.loads((O/R16/'summary/chain_summary.json').read_text())['iterates']
t = np.array([q['t_min'] for q in cs]); dn = -np.array([q['cl_c'] for q in cs]); g = np.array([q['run'] for q in cs])
for k, r in enumerate(RUNS[R16]['rounds']):
    m = g == k
    ax.plot(t[m], dn[m], 'o-', ms=2.5, lw=1.4, c=BLUES[min(k, 4)], label=f'mesh {k}')
ax.set(xlabel='wall time (min, 24 cores)', ylabel='Downforce $-C_Lc$', title='(c) Progress (a new mesh starts each colour)')
ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='lower right')
fig.tight_layout(); save(fig, O/R16/'summary/plots/result_overview.png')

# trailing edge close-up, start and best
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(*np.r_[xy0, xy0[:1]].T, c=muted, lw=1, ls=(0, (3, 2)), label='start')
ax.fill(*xyb.T, color='#9a9a9a', ec=ink, lw=.7, label='best')
up = halves(xyb)[1]; P = xyb[up]; P = P[np.argsort(P[:, 0])]
ang = np.degrees(np.arctan2(np.gradient(P[:, 1]), np.gradient(P[:, 0])))
i45 = len(P)-1-np.argmax(np.abs(ang[::-1]) < 45)  # last point flatter than 45 deg; steeper from there to the edge
x45, y45, yte = P[i45, 0], P[i45, 1], P[:, 1].max()
h = yte-y45
ax.plot([x45], [y45], 'o', c=ORANGE, ms=4)
ax.annotate('', (.522, y45), (.522, yte), arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1))
ax.text(.526, (y45+yte)/2, f'{1e3*h:.0f} mm\n= {100*h/C0:.0f} % c\nsteeper\nthan 45°', color=ORANGE, va='center', fontsize=9)
ax.set(aspect='equal', xlim=(.36, .56), ylim=(.04, .2), xlabel='x (m)', ylabel='y (m)', title='Trailing edge of the best run_016 design')
ax.add_patch(Rectangle((BOX['xmin'], BOX['ymin']), BOX['xmax']-BOX['xmin'], BOX['ymax']-BOX['ymin'], fill=False, ls='--', ec=muted, lw=.8))
ax.legend(fontsize=8.5, frameon=False, loc='upper left')
save(fig, O/R16/'summary/plots/trailing_edge.png')
print('turn-up (steeper than 45 deg) mm', 1e3*h, 'from x', x45, 'y+ max start/best', s0['yplus'].max(), s1['yplus'].max(), 'Cp min start/best', s0['cp'].min(), s1['cp'].min())

# ------------------------------------------------------------------ per-run best sections
for R, d in RUNS.items():
    if d.get('best') and (O/R/'summary/plots').exists(): per_run_section(R)

# ------------------------------------------------------------------ drag sweep
DRAG = [(R16, 'no limit'), ('run_018_area050DragCap008', '≤ 0.08'), ('run_028_area050DragCap005', '≤ 0.05'),
        ('run_017_area050DragCap003', '≤ 0.03'), ('run_029_area050DragCap002', '≤ 0.02')]
fig, ax = plt.subplots(figsize=(6.6, 4.4))
cols = [ink]+BLUES[1:]
polar(ax, [(R, f"{R[4:7]}  $C_Dc$ {lab}", c) for (R, lab), c in zip(DRAG, cols)],
      caps=[(RUNS[R]['drag_cap_cdc'], c) for (R, _), c in zip(DRAG[1:], cols[1:])])
front(ax, [R for R, _ in DRAG])
ax.set(xlim=(0, .18), ylim=(1.2, 2.4), title='Drag sweep: best design of each run (dashed: drag limit)')
save(fig, CMP/'drag_sweep.png')
contact('drag', [(R, f"{R[4:7]}: drag {lab}") for R, lab in DRAG])

# ------------------------------------------------------------------ loss sweep (+ min-loss runs)
LOSS = [(R16, 'no limit'), ('run_019_area050CptLoss075', '≤ 75 %'), ('run_020_area050CptLoss050', '≤ 50 %'),
        ('run_021_area050CptLoss025', '≤ 25 %'), ('run_030_area050CptLoss015', '≤ 15 %')]
MINLOSS = [('run_031_area050MinLossClc18CapStart', '031'), ('run_032_area050MinLossClc16CapStart', '032'),
           ('run_036_area050MinLossClc10CpMin4CapStart', '036 (Cp ≥ −4)'), ('run_037_area050MinLossClc08CpMin4CapStart', '037 (Cp ≥ −4)')]
fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.4))
ax = axs[0]
front(ax, [R for R, _ in DRAG], col='#c4c4c4', label='drag sweep (previous section)')
polar(ax, [(R, f"{R[4:7]}  loss {lab}", c) for (R, lab), c in zip(LOSS, cols)])
front(ax, [R for R, _ in LOSS]); ax.set(xlim=(0, .18), ylim=(1.2, 2.4), title='(a) Downforce vs drag')
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax = axs[1]
for (R, lab), c in zip(LOSS, cols):
    b = RUNS[R]['best']; ax.plot(100*b['loss'], -b['clc'], 'o', ms=8, mfc=c, mec=ink, mew=.8, zorder=3)
    ax.annotate(f"{R[4:7]}  {lab}", (100*b['loss'], -b['clc']), xytext=(6, -3), textcoords='offset points', fontsize=8.5, va='top')
P = sorted((100*RUNS[R]['best']['loss'], -RUNS[R]['best']['clc']) for R, _ in LOSS); ax.plot(*np.array(P).T, c=muted, lw=1)
for R, lab in MINLOSS:
    b = RUNS[R].get('best')
    if not b: continue
    ax.plot(100*b['loss'], -b['clc'], 's', ms=7, mfc='white', mec=BLUES[3], mew=1.3, zorder=4)
    ax.annotate(lab, (100*b['loss'], -b['clc']), xytext=(-6, 4), textcoords='offset points', fontsize=8, ha='right', color=BLUES[3])
ax.plot([], [], 's', ms=7, mfc='white', mec=BLUES[3], mew=1.3, label='loss minimized, downforce limit (CapStart runs)')
ax.set(xlabel='Outlet $C_{pt}$ loss (% of $q_\\infty$)', ylabel='Downforce $-C_Lc$', title='(b) Downforce vs outlet loss', xlim=(0, 2.3), ylim=(.6, 2.4))
ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True); ax.legend(fontsize=8, frameon=False, loc='lower right')
fig.tight_layout(); save(fig, CMP/'loss_sweep.png')
contact('loss', [(R, f"{R[4:7]}: loss {lab}") for R, lab in LOSS])

# ------------------------------------------------------------------ box floor (ground clearance) and free air
FLOOR = [(R16, 50), ('run_022_area050Floor37p5mm', 37.5), ('run_023_area050Floor25mm', 25), ('run_024_area050Floor12p5mm', 12.5)]
FREE = 'run_035_area050FreeAir'
fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.4))
ax = axs[0]
fc = [ink, BLUES[1], BLUES[2], BLUES[4]]
OFF_A = {R16: (-8, -6, 'right'), 'run_022_area050Floor37p5mm': (0, -10, 'center'), 'run_023_area050Floor25mm': (0, 18, 'center'),
         'run_024_area050Floor12p5mm': (6, -6, 'left')}
polar(ax, [(R, f"{R[4:7]}  floor {h:g} mm", c, OFF_A[R]) for (R, h), c in zip(FLOOR, fc)] + [(FREE, '035  free air', 'white')])
ax.set(xlim=(0, .28), ylim=(1.0, 2.5), title='(a) Downforce vs drag')
ax = axs[1]
for (R, h), c in zip(FLOOR, fc):
    b = RUNS[R]['best']; hc = b['clearance_mm']/1e3/C0
    ax.plot(hc, -b['clc'], 'o', ms=8, mfc=c, mec=ink, mew=.8, zorder=3)
    dx, dy, ha = {R16: (0, 10, 'center'), 'run_022_area050Floor37p5mm': (0, -14, 'center'), 'run_023_area050Floor25mm': (0, 10, 'center'),
                  'run_024_area050Floor12p5mm': (-8, -14, 'center')}[R]
    ax.annotate(f"{R[4:7]}  floor {h:g} mm" + (f"\n(sits at {b['clearance_mm']:.0f} mm)" if b['clearance_mm'] > h+1 else ''), (hc, -b['clc']),
                xytext=(dx, dy), textcoords='offset points', fontsize=8.5, ha=ha, va='bottom' if dy > 0 else 'top')
bf = RUNS[FREE]['best']; ax.axhline(-bf['clc'], c=muted, lw=1, ls='--'); ax.text(.105, -bf['clc']+.02, 'free air (035)', fontsize=8.5, color=muted, ha='right')
ax.set(xlabel='lowest point above ground, h/c (c = 0.5 m)', ylabel='Downforce $-C_Lc$', title='(b) Downforce vs height of the best design',
       xlim=(0, .11), ylim=(1.0, 2.5)); ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
fig.tight_layout(); save(fig, CMP/'floor_sweep.png')
contact('floor', [(R, f"{R[4:7]}: floor {h:g} mm") for R, h in FLOOR] + [(FREE, '035: free air, box floor 50 mm')], ncol=3)

# ------------------------------------------------------------------ Cp min
CPR = [(R16, None), ('run_025_area050CpMin8', -8), ('run_026_area050CpMin6', -6), ('run_027_area050CpMin4', -4)]
fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
ax = axs[0]
polar(ax, [(R, f"{R[4:7]}  " + ('no limit' if L is None else f'$C_p$ ≥ {L:g}'), c) for (R, L), c in zip(CPR, cols)])
front(ax, [R for R, _ in DRAG], col='#c4c4c4', label='drag sweep'); ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.set(xlim=(0, .18), ylim=(1.2, 2.4), title='(a) Downforce vs drag')
ax = axs[1]
for (R, L), c in zip(CPR, cols):
    s = slices(best_dir(R)); lo = s['suction']
    ax.plot(s['x'][lo]/C0, s['cp'][lo], c=c, lw=1.3, label=f"{R[4:7]} " + ('no limit' if L is None else f'≥ {L:g}') + f" (min {s['cp'].min():.1f})")
ax.invert_yaxis(); ax.set(xlabel='x/c', ylabel='$C_p$, suction side', title='(b) Suction-side pressure of the best designs', xlim=(-.05, 1.05))
ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='upper right')
ax = axs[2]
for (R, L), c in zip(CPR, cols):
    b = RUNS[R]['best']; ax.plot(b['area_pct'], -b['clc'], 'o', ms=8, mfc=c, mec=ink, mew=.8)
    ax.annotate(f"{R[4:7]}", (b['area_pct'], -b['clc']), xytext=(6, -3), textcoords='offset points', fontsize=8.5, va='top')
ax.axvline(50, c=muted, lw=.9, ls='--'); ax.text(51, 1.25, 'area limit 50 %', fontsize=8.5, color=muted)
ax.set(xlabel='section area (% of start)', ylabel='Downforce $-C_Lc$', title='(c) Downforce vs section area', xlim=(40, 100), ylim=(1.2, 2.4))
ax.grid(c=grid_c, lw=.6)
fig.tight_layout(); save(fig, CMP/'cpmin_sweep.png')
contact('cpmin', [(R, f"{R[4:7]}: " + ('no Cp limit' if L is None else f'Cp ≥ {L:g}')) for R, L in CPR], ncol=2)

# ------------------------------------------------------------------ misc
MISC = [('run_003_ffd20Pairs', '003: x + y moves, tilt free, start NACA 0012 at 3°'),
        ('run_011_naca0010LevelFront', '011: as 003, start NACA 0010 at 0°'),
        ('run_006_levelNoTilt', '006: y moves only, tilt locked, level start'),
        (R16, '016: as 006, area ≥ 50 %'),
        ('run_031_area050MinLossClc18CapStart', '031: min loss, downforce ≥ 1.8, loss-cap start'),
        ('run_007_floor30mm', '007: as 003, floor 30 mm')]
contact('misc', MISC, ncol=2, ylim=(-.012, .20))

# ------------------------------------------------------------------ ground effect: Cp around the best designs, with and without ground
cp_at = lambda v: (v+10)/11
CPMAP = LinearSegmentedColormap.from_list('cp', [(0, '#00274c'), (cp_at(-4), '#2f6fbf'), (cp_at(0), 'white'), (1, '#ffcb05')])
CPN = Normalize(-10., 1.); LEV = np.r_[np.linspace(-10, 0, 51), np.linspace(0, 1, 11)[1:]]
fig, axs = plt.subplots(1, 2, figsize=(13, 3.9))
for ax, R, title in ((axs[0], R16, '(a) In ground effect: run_016'), (axs[1], 'run_035_area050FreeAir', '(b) Free air: run_035')):
    d = best_dir(R); z = read_zones(d/'flow_000_vol.cgns'); xy = surface(d); b = RUNS[R]['best']
    for tag in ['background', 'near']:
        f = z[tag]; ax.tricontourf(tri(f, tag == 'near'), np.clip(f['CoefPressure'], -50, 2).ravel(), levels=LEV, cmap=CPMAP, norm=CPN, extend='both', zorder=1)
    free = R != R16
    ax.fill(*xy.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4)
    if free: ax.axhline(0, c=muted, lw=.8, ls=':')
    else: ax.add_patch(Rectangle((-1, -1), 3, 1, color=ground_c, lw=0, zorder=3)); ax.axhline(0, c=ink, lw=1.4, zorder=4)
    sm = slices(d)
    ax.set(aspect='equal', xlim=(-.15, .75), ylim=(-.08, .3), xlabel='x (m)', ylabel='y (m)',
           title=f"{title}\nbest design: $-C_Lc$ {-b['clc']:.2f}, lowest $C_p$ {sm['cp'].min():.1f}")
fig.colorbar(ScalarMappable(CPN, CPMAP), ax=axs, label='$C_p$', ticks=[-10, -8, -6, -4, -2, 0, 1], shrink=.75, pad=.02)
save(fig, CMP/'cp_ground_vs_free.png')
