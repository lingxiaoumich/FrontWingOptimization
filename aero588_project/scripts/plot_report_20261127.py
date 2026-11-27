"""Figures new in reports/report_20261127.html (sections 8 and 9); the figures of sections 1-7 are those of
report_20261122 (plot_report.py), unchanged. Output: optimizations/comparisons/report_20261127/.
  kkt_sections.png   the KKT designs of runs 036-039
  kkt_results.png    loss vs downforce limit, downforce vs drag against the drag sweep, suction-side Cp
  kkt_history.png    loss and PSQP's G over the KKT push (25-27 Nov), with the fixes marked
  replays.png        line-search replays: far field moving vs held (a), solver tolerance 1e-12 vs 1e-13 (b)
  rotsuction.png     rotated start (040 vs 016) and minimum suction (042)
Numbers come from runs.json (report_gather.py). MACH-Aero env, from scripts/:
  python plot_report_20261127.py ../optimizations/comparisons/report_20261127/runs.json
"""
import datetime as dt
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

O = Path(__file__).resolve().parents[1]/'optimizations'
CMP = O/'comparisons/report_20261127'; CMP.mkdir(parents=True, exist_ok=True)
RUNS = json.loads(Path(sys.argv[1]).read_text())
C0, LOSS_SCALE = 0.5, 0.005
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': 'white',
                     'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans', 'axes.titlesize': 10.5, 'axes.titlelocation': 'left'})
ink, muted, grid_c, ground_c = '#2b2b2b', '#6b6b6b', '#e4e4e4', '#d4d4d4'
BLUES = ['#93bff2', '#4387dd', '#1f62bd', '#154a8f', '#0d3366']
ORANGE = '#eb6834'
BOX = RUNS['run_016_area050']['box']
KKT = [('run_037_area050MinLossClc08CpMin4CapStart', 0.8), ('run_036_area050MinLossClc10CpMin4CapStart', 1.0),
       ('run_039_area050MinLossClc12CpMin4', 1.2), ('run_038_area050MinLossClc14CpMin4CapStart', 1.4)]
KC = dict(zip([r for r, _ in KKT], BLUES[1:]))  # ordered by the downforce limit: the report's sequential blues

def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches='tight'); plt.close(fig); print(path)

# ------------------------------------------------------------------ helpers (as plot_report.py)
def best_dir(R, b=None):
    b = b or RUNS[R]['best']; return O/R/b['round']/'best_feasible'

def surface(d):
    return np.loadtxt(d/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]

def halves(xy):
    le = int(np.argmin(xy[:, 0])); a, b = np.arange(0, le+1), np.arange(le, len(xy))
    return (a, b) if xy[a, 1].mean() < xy[b, 1].mean() else (b, a)

def slices(d):
    lines = (d/'flow_000_slices.dat').read_text().split('\n')
    n = int(lines[3].split()[2])
    P = np.array([l.split() for l in lines[5:5+n]], float)
    xy = surface(d)
    near = np.argmin((P[:, 0, None]-xy[None, :, 0])**2+(P[:, 1, None]-xy[None, :, 1])**2, axis=1)
    P = P[np.argsort(near, kind='stable')]; near = np.sort(near)
    lo, up = halves(xy)
    return dict(x=P[:, 0], y=P[:, 1], cp=P[:, 6], suction=np.isin(near, lo))

def scene(ax, box=BOX, xlim=(-.03, .56)):
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=.8, zorder=2))
    ax.add_patch(Rectangle((xlim[0]-1, -1), xlim[1]-xlim[0]+2, 1, color=ground_c, lw=0, zorder=1)); ax.axhline(0, c=ink, lw=1.4)
    ax.set_aspect('equal')

def shape_grid(ax):
    ax.xaxis.set_major_locator(MultipleLocator(.1)); ax.xaxis.set_minor_locator(MultipleLocator(.05))
    ax.yaxis.set_major_locator(MultipleLocator(.05)); ax.yaxis.set_minor_locator(MultipleLocator(.025))
    ax.grid(which='both', color=grid_c, lw=.6); ax.set_axisbelow(True)

def section_axes(ax, xlim=(-.03, .56), ylim=(-.012, .20)):
    scene(ax, xlim=xlim); ax.set(xlim=xlim, ylim=ylim); shape_grid(ax); ax.tick_params(labelsize=8)

# ------------------------------------------------------------------ 1 KKT designs
fig, axs = plt.subplots(2, 2, figsize=(10.4, 5.4))
for ax, (R, lim) in zip(axs.ravel(), KKT):
    b = RUNS[R]['best']; xy = surface(best_dir(R))
    section_axes(ax); ax.fill(*xy.T, color='#9a9a9a', ec=ink, lw=.6, zorder=4)
    ax.set_title(f"{R[4:7]}: downforce ≥ {lim:g}, $C_p$ ≥ −4\n$-C_Lc$ {-b['clc']:.2f}, $C_Dc$ {b['cdc']:.4f}, loss {100*b['loss']:.3f} %, "
                 f"lowest $C_p$ {b['cp_min_ks']:.2f}", fontsize=9)
fig.text(.5, -.01, 'KKT design of each run. Dashed: regulation box. Grid: 50 mm in x, 25 mm in y.', ha='center', fontsize=8.5, color=muted)
fig.tight_layout(); save(fig, CMP/'kkt_sections.png')

# ------------------------------------------------------------------ 2 KKT results
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
ax = axs[0]
P = np.array([(-RUNS[R]['best']['clc'], 100*RUNS[R]['best']['loss']) for R, _ in KKT])
ax.plot(*P.T, c=muted, lw=1, zorder=2)
for (R, lim), (x, y) in zip(KKT, P):
    ax.plot(x, y, 'o', ms=8, mfc=KC[R], mec=ink, mew=.8, zorder=3)
    ax.annotate(R[4:7], (x, y), xytext=(8, -4) if R.startswith('run_038') else (-8, 4), textcoords='offset points', fontsize=8.5,
                ha='left' if R.startswith('run_038') else 'right', va='top' if R.startswith('run_038') else 'baseline')
for a in RUNS['run_038_area050MinLossClc14CpMin4CapStart'].get('kkt_alternatives', []):
    ax.plot(-a['clc'], 100*a['loss'], 'o', ms=8, mfc='white', mec=KC['run_038_area050MinLossClc14CpMin4CapStart'], mew=1.3, zorder=3)
    ax.annotate('038, second KKT point (hollow)', (-a['clc'], 100*a['loss']), xytext=(-10, 8), textcoords='offset points', fontsize=8, ha='right', color=muted)
ax.set(xlabel='Downforce $-C_Lc$ (= its limit)', ylabel='Outlet $C_{pt}$ loss (% of $q_\\infty$), corrected', title='(a) Loss at the KKT point',
       xlim=(.7, 1.5), ylim=(0, .5)); ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
ax = axs[1]
DRAG = ['run_016_area050', 'run_018_area050DragCap008', 'run_028_area050DragCap005', 'run_017_area050DragCap003', 'run_029_area050DragCap002']
F = sorted((RUNS[R]['best']['cdc'], -RUNS[R]['best']['clc']) for R in DRAG)
ax.plot(*np.array(F).T, c='#c4c4c4', lw=1, zorder=1, label='drag sweep 016–029: most downforce, drag limit')
for R in DRAG: ax.plot(RUNS[R]['best']['cdc'], -RUNS[R]['best']['clc'], 'o', ms=6, mfc='#c4c4c4', mec='#9a9a9a', zorder=2)
b27 = RUNS['run_027_area050CpMin4']['best']
ax.plot(b27['cdc'], -b27['clc'], 's', ms=7, mfc='white', mec=muted, mew=1.2, zorder=3, label='027: most downforce, $C_p$ ≥ −4')
for R, lim in KKT:
    b = RUNS[R]['best']; ax.plot(b['cdc'], -b['clc'], 'o', ms=8, mfc=KC[R], mec=ink, mew=.8, zorder=4)
    ax.annotate(R[4:7], (b['cdc'], -b['clc']), xytext=(6, -3), textcoords='offset points', fontsize=8.5, va='top')
ax.plot([], [], 'o', ms=8, mfc=BLUES[2], mec=ink, mew=.8, label='036–039: least loss, $C_p$ ≥ −4')
ax.set(xlabel='Drag $C_Dc$', ylabel='Downforce $-C_Lc$', title='(b) Downforce vs drag', xlim=(0, .1), ylim=(.6, 2.1))
ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True); ax.legend(fontsize=7.8, frameon=False, loc='lower right')
ax = axs[2]
for R, lim in KKT:
    s = slices(best_dir(R)); lo = s['suction']
    ax.plot(s['x'][lo]/C0, s['cp'][lo], c=KC[R], lw=1.4, label=f"{R[4:7]} (≥ {lim:g})")
ax.axhline(-4, c=muted, lw=.9, ls='--'); ax.text(1.03, -4.08, '$C_p$ limit −4 (KS aggregate)', fontsize=8, color=muted, ha='right', va='bottom')
ax.set(xlabel='x/c (c = 0.5 m)', ylabel='$C_p$, suction side', title='(c) Suction-side pressure of the KKT designs', xlim=(-.05, 1.05), ylim=(.7, -4.6))
ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='center right', title='downforce limit', title_fontsize=8)
fig.tight_layout(); save(fig, CMP/'kkt_results.png')

# ------------------------------------------------------------------ 3 history of the KKT push
LINE = re.compile(r'NIT=\s*(\d+)\s+NFV=\s*(\d+)\s+NFG=\s*\d+\s+F=\s*(\S+)\s+C=\s*(\S+)\s+G=\s*(\S+)')

def series(root):
    """[(time, F, C, G)] per PSQP iteration of a round. PSQP's first function call is the last evaluation whose
    objective matches its first F (the derivative check also solves the start design); iteration k is then NFV_k -
    NFV_0 evaluations later, as in plot_run_chain.py. (The history's objective_value does not always equal PSQP's F
    later in a round, so the evaluations are counted, not matched.)"""
    rows = {}
    for m in LINE.finditer((root/'PSQP.out').read_text()): rows[int(m.group(1))] = (int(m.group(2)),)+tuple(map(float, m.group(3, 4, 5)))
    H = [json.loads(l) for l in open(root/'history.jsonl') if l.strip()]
    if not rows or not H: return []
    t0 = (root/'history.jsonl').stat().st_mtime-H[-1]['elapsed_wall_time_s']
    L = [rows[k] for k in sorted(rows)]; nfv0, F0 = L[0][0], L[0][1]; span = L[-1][0]-nfv0
    cands = [i for i, h in enumerate(H) if abs(h.get('objective_value', 1e9)-F0) < 1e-5*max(1., abs(F0)) and i+span < len(H)+1]
    e0 = max(cands) if cands else 0
    return [(dt.datetime.fromtimestamp(t0+H[min(e0+nfv-nfv0, len(H)-1)]['elapsed_wall_time_s']), F, C, G) for nfv, F, C, G in L]

def push_rounds(link):
    rs = []
    for c in link.rglob('command.txt'):
        txt = c.read_text()
        if ('KKT continuation' in txt or 'same-grid continuation' in txt) and (c.parent/'PSQP.out').exists() and (c.parent/'history.jsonl').exists():
            rs.append(c.parent)
    return rs

FIXES = [(dt.datetime(2026, 11, 25, 22, 17), 'loss against an inlet plane'), (dt.datetime(2026, 11, 25, 23, 2), 'solves to 10$^{-12}$'),
         (dt.datetime(2026, 11, 26, 1, 19), 'loss-plane points off the cell centres'), (dt.datetime(2026, 11, 26, 14, 48), 'far field held fixed'),
         (dt.datetime(2026, 11, 26, 19, 36), '038: steps ≤ 0.5 mm'), (dt.datetime(2026, 11, 27, 0, 11), '038: 2 mm branch'),
         (dt.datetime(2026, 11, 27, 15, 13), '038: solves to 10$^{-13}$, same grid')]
KKT_ROUNDS = {(R, r) for R, _ in KKT for r in [RUNS[R]['best']['round']]+[a['round'] for a in RUNS[R].get('kkt_alternatives', [])]}
fig, axs = plt.subplots(2, 1, figsize=(14, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.15]))
for R, lim in KKT:
    link = O/R
    for root in push_rounds(link):
        s = series(root)
        if not s: continue
        t = [p[0] for p in s]; F = np.array([p[1] for p in s]); G = np.array([p[3] for p in s])
        axs[0].plot(t, 100*LOSS_SCALE*F, c=KC[R], lw=1.1)
        ok = G < 1e10; axs[1].plot(np.array(t)[ok], G[ok], c=KC[R], lw=.9)
        rel = str(root.relative_to(link))
        st = json.loads((root/'run_complete.json').read_text())['inform'] if (root/'run_complete.json').exists() else {}
        if isinstance(st, dict) and st.get('value') == 4:
            final = (R, rel) in KKT_ROUNDS
            for ax, v in ((axs[0], 100*LOSS_SCALE*F[-1]), (axs[1], G[ok][-1])):
                ax.plot(t[-1], v, '*', ms=13, mfc=KC[R] if final else 'white', mec=ink, mew=.8, zorder=5)
    axs[1].plot([], [], c=KC[R], lw=1.4, label=f"{R[4:7]}: downforce ≥ {lim:g}")
axs[1].plot([], [], '*', ms=12, mfc=BLUES[2], mec=ink, mew=.8, ls='none', label='KKT met (PSQP exit 4)')
axs[1].plot([], [], '*', ms=12, mfc='white', mec=ink, mew=.8, ls='none', label='KKT met before the far-field fix (036 r06), superseded')
axs[1].axhline(1e-4, c=ORANGE, lw=1, ls='--'); axs[1].text(dt.datetime(2026, 11, 25, 21, 30), 1.15e-4, 'KKT tolerance 10$^{-4}$', color=ORANGE, fontsize=8.5, va='bottom')
for k, (t, lab) in enumerate(FIXES, 1):
    for ax in axs: ax.axvline(t, c=muted, lw=.8, ls=':')
    axs[0].text(t, .95 if k % 2 else .86, str(k), transform=axs[0].get_xaxis_transform(), ha='center', va='center', fontsize=9, color=ink,
                bbox=dict(boxstyle='circle,pad=.18', fc='white', ec=muted, lw=.7), zorder=6)
axs[0].set(ylabel='Outlet $C_{pt}$ loss (% of $q_\\infty$)', title='(a) Loss of each PSQP iterate (its definition changes at 1 and 3)', ylim=(0, 1.0))
axs[1].set(yscale='log', ylabel='PSQP G: largest component of\nthe Lagrangian gradient', ylim=(1e-5, 3), title='(b) Distance from KKT')
for ax in axs: ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
axs[1].legend(fontsize=8, frameon=False, loc='upper right', ncol=2)
axs[1].xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 6, 12, 18])); axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%d Nov\n%H:%M'))
fig.text(.01, -.02, '   '.join(f'{k} {lab}' for k, (_, lab) in enumerate(FIXES, 1)), fontsize=8.5, color=muted, wrap=True)
fig.tight_layout(); save(fig, CMP/'kkt_history.png')

# ------------------------------------------------------------------ 4 replays
RP = CMP/'replays'
rd = lambda t: json.loads((RP/t/'linesearch.json').read_text())
fig, axs = plt.subplots(1, 2, figsize=(13, 4.3))
ax = axs[0]; T0 = 3e-5
for tag, col, lab in (('G12', ORANGE, 'far field moving with the wing (before 26 Nov)'), ('G12fft', BLUES[2], 'far field held, top tilted (after)')):
    d = rd(tag); base = d[0]['loss']
    t = np.array([max(p['t'], T0) for p in d]); y = 1e7*(np.array([p['loss'] for p in d])-base)
    o = np.argsort(t, kind='stable')
    ax.plot(t, y, 'o', ms=6, mfc=col, mec=ink, mew=.6, label=lab, zorder=3)
    ax.plot(t[o][t[o] > T0], y[o][t[o] > T0], c=col, lw=.9, alpha=.7)
ax.set(xscale='log', xlim=(2e-5, 2), xlabel='step fraction t along the line search (t = 0 drawn at the left edge)',
       ylabel='loss − loss at t = 0  (10$^{-7}$)', title='(a) Replay of a failed line search (037 r07, evaluation 172)')
ax.axhline(0, c=muted, lw=.7); ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='lower left')
ax.annotate('a jump of 5×10$^{-7}$ between t = 5×10$^{-4}$\nand 10$^{-3}$: every trial point reads\nworse than t = 0, so the step fails', (1e-3, 4.6), xytext=(1.5e-3, -9),
            fontsize=8, color=ORANGE, arrowprops=dict(arrowstyle='->', color=ORANGE, lw=.8))
ax = axs[1]
for tag, col, lab in (('N12', ORANGE, 'solves to 10$^{-12}$ of the free-stream residual'), ('N13', BLUES[2], 'solves to 10$^{-13}$')):
    d = rd(tag); base = np.mean([p['loss'] for p in d if p['t'] == 0])
    t = np.array([p['t'] for p in d]); y = 1e12*(np.array([p['loss'] for p in d])-base)
    ax.plot(t, y, 'o', ms=6, mfc=col, mec=ink, mew=.6, label=lab, zorder=3)
ax.set(xlabel='step fraction t (0.36 mm at t = 1)', ylabel='loss − mean at t = 0  (10$^{-12}$)',
       title='(b) Loss along 0.36 mm near 038\'s KKT point, solves in sequence')
ax.axhline(0, c=muted, lw=.7); ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='upper left')
fig.tight_layout(); save(fig, CMP/'replays.png')

# ------------------------------------------------------------------ 5 rotated start and minimum suction
R16, R40, R42, R27 = 'run_016_area050', 'run_040_area050RotStart', 'run_042_area050MinSuctionClc12CapStart', 'run_027_area050CpMin4'
fig, axs = plt.subplots(2, 2, figsize=(13, 6.6))
ax = axs[0, 0]; section_axes(ax, ylim=(-.012, .245))
for R, col, lab in ((R16, muted, 'level start (016–039)'), (R40, BLUES[2], 'rotated start (040, 041): 8.3° nose down, TE at 2/3 of the box height')):
    xy = surface(O/R/'r00/baseline'); ax.plot(*np.r_[xy, xy[:1]].T, c=col, lw=1.3, label=lab)
ax.legend(fontsize=8, frameon=False, loc='upper left'); ax.set_title('(a) Start sections', fontsize=9.5)
ax = axs[0, 1]; section_axes(ax, ylim=(-.012, .245))
for R, col, fill in ((R16, muted, False), (R40, BLUES[2], True)):
    b = RUNS[R]['best']; xy = surface(best_dir(R))
    lab = f"{R[4:7]}: $-C_Lc$ {-b['clc']:.3f}, $C_Dc$ {b['cdc']:.3f}, {len(RUNS[R]['rounds'])} meshes"
    if fill: ax.fill(*xy.T, color=col, alpha=.35, ec=col, lw=.8, zorder=4, label=lab)
    else: ax.plot(*np.r_[xy, xy[:1]].T, c=ink, lw=1, ls=(0, (3, 2)), zorder=5, label=lab)
ax.legend(fontsize=8, frameon=False, loc='upper left'); ax.set_title("(b) Best designs under run_016's limits", fontsize=9.5)
ax = axs[1, 0]; section_axes(ax, ylim=(-.012, .245))
b42 = RUNS[R42]['best']; s0 = O/R42/'r00/baseline'; r0 = json.loads((s0/'result.json').read_text())
xy0 = surface(s0); xyb = surface(best_dir(R42))
ax.plot(*np.r_[xy0, xy0[:1]].T, c=muted, lw=1.2, ls=(0, (3, 2)), label=f"start (027's best, new grid): $C_{{p,min}}$ {r0['cp_min_ks']:.2f}, $-C_Lc$ {-.5*r0['cl']:.2f}, $C_Dc$ {.5*r0['cd']:.3f}")
ax.fill(*xyb.T, color=BLUES[2], alpha=.35, ec=BLUES[2], lw=.8, zorder=4,
        label=f"042 best: $C_{{p,min}}$ {b42['cp_min_ks']:.2f}, $-C_Lc$ {-b42['clc']:.2f}, $C_Dc$ {b42['cdc']:.3f}")
ax.legend(fontsize=8, frameon=False, loc='upper left'); ax.set_title('(c) Least suction peak with downforce ≥ 1.2 (042)', fontsize=9.5)
ax = axs[1, 1]
for d, col, ls, lab in ((s0, muted, (0, (3, 2)), 'start'), (best_dir(R42), BLUES[2], '-', '042 best'),
                        (best_dir('run_039_area050MinLossClc12CpMin4'), BLUES[3], ':', '039 KKT (least loss, $C_p$ ≥ −4, same downforce)')):
    s = slices(d)
    for side, lw in ((s['suction'], 1.5), (~s['suction'], .9)):
        ax.plot(s['x'][side]/C0, s['cp'][side], c=col, ls=ls, lw=lw, label=lab if lw > 1 else None)
ax.invert_yaxis(); ax.grid(c=grid_c, lw=.6); ax.legend(fontsize=8, frameon=False, loc='upper right')
ax.set(xlabel='x/c (c = 0.5 m)', ylabel='$C_p$ (thick: suction side)', title='(d) Surface pressure', xlim=(-.05, 1.05))
fig.tight_layout(); save(fig, CMP/'rotsuction.png')
