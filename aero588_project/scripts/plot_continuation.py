"""Plot the PSQP continuation chains started on October 31 (run_continuation.py).

A chain's path is its PSQP major iterates, stage after stage: an iterate is an
evaluation at which the optimizer asked for a gradient (gradients.jsonl), so
line-search trials, failed solves and rejected geometries are left out. Iterations
are counted cumulatively from the sweep restart point (the October 29 final design).
Chains made with --min-local-thickness name their parent chain's evaluation in
start_source.txt and are drawn from that point on the parent's path.

Every iterate's surface is rebuilt with the run's own pyGeo setup, so the minimum
chord-normal local thickness (vertex pairs k / 272-k, k = 15..124, the quantity the
2 mm constraint bounds) is shown for all chains, including those run without it.
Run with the MACH-Aero environment active:
  python scripts/plot_continuation.py --out optimizations/continuation_summary_20261031 \
      optimizations/continuation_xmax2.0_20261031 optimizations/continuation_xmax5.0_20261031 \
      optimizations/continuation_tmin2mm_xmax2.0_20261031 optimizations/continuation_tmin2mm_xmax5.0_20261031
"""
import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

P = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('chains', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
out = a.out.resolve()/'plots'; out.mkdir(parents=True, exist_ok=True)
chains = [c.resolve() for c in a.chains]

# The optimizer's own geometry setup (everything before the CFD options), executed in a
# scratch root so its setup files are not written into a run folder.
scratch = Path(tempfile.mkdtemp())
shutil.copytree(P/'setups/slip_ground_xy5_cd005/source', scratch/'source')
shutil.copy(P/'scripts/run_slip_optimization.py', scratch/'source')
sys.path.insert(0, str(scratch/'source'))
argv, sys.argv = sys.argv, ['x', '--root', str(scratch), '--ffd-stations', '5', '--shape-directions', 'xy',
                            '--min-local-thickness', '0.002']
setup = {}
exec(compile((scratch/'source/run_slip_optimization.py').read_text().split('options={')[0], 'setup', 'exec'), setup)
sys.argv = argv
geo, metrics, local_thickness, c0 = setup['geo'], setup['metrics'], setup['local_thickness'], setup['c0']

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

def surface(design):
    geo.setDesignVars({k: np.array(v) for k, v in design.items()})
    return geo.update('surface_check')

def iterates(stage):
    h = [json.loads(l) for l in (stage/'history.jsonl').read_text().splitlines()]
    g = [json.loads(l)['evaluation_count']-1 for l in (stage/'gradients.jsonl').read_text().splitlines()] \
        if (stage/'gradients.jsonl').exists() else []
    # Before the INVALID_GEOMETRY fix, a rejected trial re-ran the adjoint at the last
    # valid evaluation, repeating its evaluation_count; count each evaluation once.
    g = list(dict.fromkeys(g))
    return [h[i] for i in g if i < len(h)], h

def path(chain):
    """Major iterates of every stage; a later stage's first iterate repeats its start design."""
    its, stages = [], []
    for k, stage in enumerate(sorted(chain.glob('stage_*'))):
        major, h = iterates(stage)
        s = json.loads((stage/'stage_summary.json').read_text()) if (stage/'stage_summary.json').exists() else {}
        cmd = (stage/'command.txt').read_text()
        cap = float(cmd.split('--drag-cap ')[1].split()[0])
        if k: major = major[1:]
        for r in major: its.append(dict(r, stage=k, cap=cap, stage_dir=stage))
        stages.append(dict(s, first=len(its)-len(major), cap=cap, rejected=(stage/'run.log').read_text().count('INVALID_GEOMETRY')))
    return its, stages

def enrich(its, remeshed):
    """Surface of every iterate. A remeshed chain's design variables are relative to its own
    grid, so its surfaces come from the iterate files it saved (--save-iterate-flows) and its
    thickness from the optimizer's record."""
    for r in its:
        if remeshed:
            f = r['stage_dir']/f"iterates/eval_{r['evaluation']:03d}/surface_geometry.csv"
            if r['evaluation'] == 0 or not f.exists(): f = r['stage_dir']/'baseline/surface_geometry.csv'
            pts = np.loadtxt(f, delimiter=',', skiprows=1)
        else:
            pts = surface(r['design'])
        r['pts'] = pts
        t = local_thickness(pts)
        r['tloc'] = r.get('minimum_local_thickness_m', float(t.min())) if remeshed else float(t.min())
        r['tloc_k'] = int(setup['local_pairs'][np.argmin(t)])
    return its

runs = {}
for c in chains:
    its, stages = path(c)
    if not its: print('skip', c.name, '(no iterates yet)'); continue
    remeshed = '--baseline' in (sorted(c.glob('stage_*'))[0]/'command.txt').read_text()
    kind = 'remesh' if remeshed else 'tmin' if 'tmin' in c.name else 'free'
    src = c/'start_source.txt'
    words = src.read_text().split() if src.exists() else []
    parent = (words[0] if len(words) == 4 else c.name.replace('tmin2mm_', '')) if words else None
    runs[c.name] = dict(its=enrich(its, remeshed), stages=stages, offset=None if parent else 0, kind=kind,
                        xmax=c.name.split('xmax')[1].split('_')[0], constrained=kind != 'free',
                        parent=parent, start=(int(words[-3].split('_')[1]), int(words[-1])) if words else None, dir=c)
# A chain starts from an iterate of its parent: place it at that iterate's position.
while any(r['offset'] is None for r in runs.values()):
    for r in runs.values():
        if r['offset'] is None and runs[r['parent']]['offset'] is not None:
            par = runs[r['parent']]; k, ev = r['start']
            idx = [i for i, q in enumerate(par['its']) if q['stage'] == k and q['evaluation'] <= ev]
            r['offset'] = par['offset']+idx[-1]
start = [r for r in runs.values() if r['kind'] == 'free'][0]['its'][0]
naca = surface(setup['xinitial'])

# Style follows plot_optimization_trajectory.py: blue light->dark for iteration order.
ramp = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
cmap = LinearSegmentedColormap.from_list('early_to_late', ramp)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid = '#2b2b2b', '#6b6b6b', '#e4e4e4'
hue = {'2.0': '#2a78d6', '5.0': '#eb6834', '1.0': '#1baf7a'}  # chain identity: XMAX (categorical slots 1-3)
style = {'free': dict(marker='o', ms=4.5, lw=1.1, ls=(0, (4, 2))), 'tmin': dict(marker='s', ms=4.5, lw=1.6, ls='-'),
         'remesh': dict(marker='^', ms=5.5, lw=2.2, ls='-')}
pinched = '#b5b5b5'
def save(fig, name):
    fig.savefig(out/(name+'.png'), dpi=170, bbox_inches='tight')
    fig.savefig(out/(name+'.svg'), bbox_inches='tight'); plt.close(fig)
def label(r):
    return f"XMAX {r['xmax']}, " + {'free': 'no local-thickness limit', 'tmin': 'local t ≥ 2 mm, warped grid',
                                   'remesh': 'local t ≥ 2 mm, remeshed grid'}[r['kind']]

# ---------------------------------------------------------------- performance history
panels = [('(a) Downforce', 'Downforce $-C_Lc$', lambda q: -q['cl']*c0),
          ('(b) Drag against its cap (cap doubles when hit)', '$C_Dc$', lambda q: q['cd']*c0),
          ('(c) Thinnest local section (pairs k = 15–124)', 'Minimum local thickness (mm)', lambda q: q['tloc']*1e3),
          ('(d) Distance to the box x limits (goal: 0)', 'x margin (mm)', None),
          ('(e) Drag polar of the iterates', '$C_Dc$', None),
          ('(f) Ground clearance', 'Minimum clearance (mm)', lambda q: q['minimum_clearance_m']*1e3)]
fig, axs = plt.subplots(2, 3, figsize=(16, 9.2)); axs = axs.ravel()
for name, r in runs.items():
    its = r['its']; n = np.arange(len(its))+r['offset']; col = hue[r['xmax']]
    ok = np.array([q['tloc'] >= .002-5e-6 for q in its])  # the driver's restart tolerance
    mk = style[r['kind']]
    for ax, (title, ylab, f) in zip(axs, panels):
        if f is None: continue
        y = np.array([f(q) for q in its])
        ax.plot(n, y, c=col, **{k: v for k, v in mk.items() if k != 'marker'}, zorder=2)
        ax.scatter(n[ok], y[ok], c=col, marker=mk['marker'], s=22, lw=0, zorder=3)
        ax.scatter(n[~ok], y[~ok], facecolor='white', edgecolor=pinched if not r['constrained'] else col,
                   marker=mk['marker'], s=22, lw=1, zorder=3)
    lead = np.array([(q['minimum_x_m']+.125)*1e3 for q in its]); trail = np.array([(.625-q['maximum_x_m'])*1e3 for q in its])
    axs[3].plot(n, lead, c=col, **{k: v for k, v in mk.items() if k != 'marker'})
    axs[3].plot(n, trail, c=col, alpha=.45, **{k: v for k, v in mk.items() if k != 'marker'})
    axs[4].plot([-q['cl']*c0 for q in its], [q['cd']*c0 for q in its], c=col, **mk, mfc=col, mec='white', mew=.6)
    axs[0].plot([], [], c=col, **{k: v for k, v in mk.items()}, label=label(r))
# caps: one step line per unconstrained chain path (both chains doubled at the same iteration)
for name, r in runs.items():
    if r['constrained']: continue
    n = np.arange(len(r['its'])); cap = np.array([q['cap']*c0 for q in r['its']])
    axs[1].step(n, cap, where='post', c=muted, lw=1, ls='--')
    for x, y in zip(n[1:][np.diff(cap) > 0], cap[1:][np.diff(cap) > 0]):
        axs[1].annotate(f'cap doubled → {y:.3f}', (x, y), xytext=(6, 4), textcoords='offset points', fontsize=8, color=muted)
    break
# A remeshed chain re-solves its parent's design on a fresh grid: the jump is the grid, not progress.
for r in [r for r in runs.values() if r['kind'] == 'remesh'][:1]:
    k, ev = r['start']
    par = [json.loads(l) for l in (runs[r['parent']]['dir']/f'stage_{k:02d}/history.jsonl').read_text().splitlines()][ev]
    fresh = json.loads((sorted(r['dir'].glob('stage_*'))[0]/'history.jsonl').read_text().splitlines()[0])
    axs[0].annotate(f"remeshed: the same design gives\n{-fresh['cl']*c0:.3f} m on a fresh grid,\n{-par['cl']*c0:.3f} m on the warped one "
                    f"({fresh['cl']/par['cl']-1:+.1%})", (r['offset'], -fresh['cl']*c0), xytext=(-150, 10),
                    textcoords='offset points', fontsize=8, color=ink, arrowprops=dict(arrowstyle='->', color=muted, lw=.8))
axs[2].axhline(2, c=ink, lw=1, ls='--'); axs[2].text(.5, 2, ' 2 mm floor (constrained chains)', va='bottom', fontsize=8, color=ink)
axs[2].set_yscale('log')
axs[3].axhline(0, c=ink, lw=1.2)
axs[3].text(.99, .97, 'dark: leading edge (x ≥ −0.125 m)\nlight: trailing edge (x ≤ 0.625 m)', transform=axs[3].transAxes,
            ha='right', va='top', fontsize=8, color=muted)
axs[5].axhline(50, c=ink, lw=1, ls='--'); axs[5].text(.5, 50, ' 50 mm limit', va='bottom', fontsize=8, color=ink)
for capc in [.025, .05, .1, .2]:
    axs[4].axhline(capc, c=muted, lw=1, ls='--'); axs[4].text(axs[4].get_xlim()[0], capc, f' cap {capc}', va='bottom', fontsize=8, color=muted)
for ax, (title, ylab, _) in zip(axs, panels):
    ax.set_title(title, fontsize=10.5); ax.set_ylabel(ylab); ax.grid(c=grid, lw=.7); ax.set_axisbelow(True)
    ax.set_xlabel('Downforce $-C_Lc$' if ax is axs[4] else 'PSQP major iteration from the sweep restart point')
axs[0].scatter([], [], facecolor='white', edgecolor=pinched, s=22, label='iterate thinner than 2 mm somewhere (hollow)')
axs[0].legend(fontsize=8.5, frameon=False, loc='lower right')
fig.suptitle('Continuation from the October 29 design: performance, constraints and progress toward the box', fontsize=12.5, y=.995)
fig.text(.5, -.01, '$C_Lc = L\'/q_\\infty$, $C_Dc = D\'/q_\\infty$ in metres, so the reference chord drops out. Each stage restarts PSQP '
         '(10 iterations) from the previous stage\'s best design; line-search trials, failed solves and rejected (self-intersecting) '
         'sections are not plotted.', ha='center', fontsize=8.5, color=muted, wrap=True)
fig.tight_layout(); save(fig, 'performance_history')

# ---------------------------------------------------------------- control points
lattice0 = np.load(P/'optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/control_points.npz')['baseline'][:, :, 0, :2]
xy0 = np.loadtxt(P/'setups/slip_ground_xy5_cd005/source/airfoil_coordinates.csv', delimiter=',', skiprows=1)
pivot = np.array([.125, xy0[136, 1]+.125*np.sin(np.deg2rad(3.))])
def read_lattice(path):
    lines = path.read_text().split('\n'); dims = [int(v) for v in lines[1].split()]
    vals = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
    return np.stack([v.reshape(dims, order='F') for v in vals], -1)[:, :, 0, :2]
def lattice(design, base=None):
    """Control points (station, lower/upper, xy): rotate the reference lattice about the pivot,
    then add the shape moves. A remeshed chain's reference is its deformed starting lattice."""
    L0, piv = (lattice0, pivot) if base is None else base
    t = np.deg2rad(-(design['incidence'][0]+3.)); R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    shape = np.stack([np.reshape(design['shape_x'], L0.shape[:2]), np.reshape(design['shape_y'], L0.shape[:2])], -1)
    return (L0-piv)@R.T+piv+shape
for r in runs.values():
    cmd = (sorted(r['dir'].glob('stage_*'))[0]/'command.txt').read_text().split()
    if '--baseline' in cmd:
        bj = Path(cmd[cmd.index('--baseline')+1])
        r['base'] = (read_lattice(bj.parent/'ffd.xyz'), np.array(json.loads(bj.read_text())['pivot_m']))
    else: r['base'] = None
def lineage(r):
    """(run, iterate, cumulative index) from the Oct 29 design to the end of run r."""
    own = [(r, q, r['offset']+i) for i, q in enumerate(r['its'])]
    if r['parent'] is None: return own
    up = [e for e in lineage(runs[r['parent']]) if e[2] <= r['offset']]
    return up+own[1:]
box = dict(xmin=-.125, xmax=.625, ymin=.05, ymax=.35); nst = lattice0.shape[0]
names = [[f'L{i+1}', f'U{i+1}'] for i in range(nst)]
phase = {'free': 'no thickness limit', 'tmin': '2 mm floor, warped grid', 'remesh': '2 mm floor, remeshed grid'}
L0 = lattice(setup['xinitial'])
terminal = [r for r in runs.values() if not any(o['parent'] == r['dir'].name for o in runs.values())]
for r in terminal:
    seq = lineage(r)
    Pts = np.array([lattice(q['design'], run['base']) for run, q, _ in seq]); ncum = [i for _, _, i in seq]
    norm = Normalize(0, max(ncum)); colors = cmap(norm(ncum))
    fig = plt.figure(figsize=(13, 10.2))
    gs = fig.add_gridspec(3, nst, height_ratios=[2.3, 1, 1], hspace=.42, wspace=.12)
    ax = fig.add_subplot(gs[0, :])
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=1.1))
    ax.text(box['xmax']-.005, box['ymax']-.008, 'Allowed surface box', ha='right', va='top', color=muted, fontsize=9)
    ax.axhline(0, c=ink, lw=2); ax.text(box['xmax']-.005, .004, 'Slip ground', ha='right', va='bottom', color=muted, fontsize=9)
    ax.fill(naca[:, 0], naca[:, 1], color='#e6e6e6', lw=0, label='NACA 0012 baseline')
    ax.plot(start['pts'][:, 0], start['pts'][:, 1], c=muted, lw=.9, ls=':', label=f'Sweep start, Oct 29 final ({-start["cl"]*c0:.3f} m)')
    seen = []
    for run, q, _ in seq:  # last iterate of each phase in this lineage
        if run is not (seen[-1][0] if seen else None): seen.append([run, q])
        else: seen[-1][1] = q
    for run, q in seen:
        ax.plot(q['pts'][:, 0], q['pts'][:, 1], c=hue[run['xmax']], lw=style[run['kind']]['lw'], ls=style[run['kind']]['ls'],
                label=f"End of XMAX {run['xmax']}, {phase[run['kind']]}: {-q['cl']*c0:.3f} m, t_min {q['tloc']*1e3:.2f} mm")
    for lat, st in [(L0, dict(c='#c8c8c8', lw=.7, ls=':')), (Pts[-1], dict(c='#555555', lw=.7))]:
        for s_ in range(2): ax.plot(lat[:, s_, 0], lat[:, s_, 1], **st)
        for i in range(nst): ax.plot(lat[i, :, 0], lat[i, :, 1], **st)
    for i in range(nst):
        for s_ in range(2): ax.scatter(Pts[:, i, s_, 0], Pts[:, i, s_, 1], c=colors, s=9, lw=0, zorder=5)
        ax.annotate(names[i][1], Pts[-1, i, 1], xytext=(0, 7), textcoords='offset points', ha='center', fontsize=8, color=ink)
        ax.annotate(names[i][0], Pts[-1, i, 0], xytext=(0, -12), textcoords='offset points', ha='center', fontsize=8, color=ink)
    ax.plot([], [], c='#c8c8c8', lw=.7, ls=':', label='NACA FFD lattice'); ax.plot([], [], c='#555555', lw=.7, label='Latest FFD lattice')
    ax.set(xlim=(box['xmin']-.02, box['xmax']+.02), ylim=(-.012, box['ymax']+.012), xlabel='x (m)', ylabel='Height above ground (m)')
    ax.set_aspect('equal'); ax.legend(loc='upper left', bbox_to_anchor=(.02, .97), fontsize=8, frameon=False)
    path_txt = ' → '.join(f"XMAX {run['xmax']} {run['kind']}" for run, _ in seen)
    ax.set_title(f'FFD control points over {max(ncum)} PSQP iterations from the Oct 29 design ({path_txt})', fontsize=10.5)
    d = (Pts-L0)*1000; lim = np.abs(d).max()*1.08
    for row, s_ in [(1, 1), (2, 0)]:
        for i in range(nst):
            z = fig.add_subplot(gs[row, i])
            z.axhline(0, c=grid, lw=.8, zorder=0); z.axvline(0, c=grid, lw=.8, zorder=0)
            z.plot(d[:, i, s_, 0], d[:, i, s_, 1], c='#c8c8c8', lw=.6, zorder=1)
            z.scatter(d[:, i, s_, 0], d[:, i, s_, 1], c=colors, s=12, lw=0, zorder=2)
            z.set(xlim=(-lim, lim), ylim=(-lim, lim)); z.set_aspect('equal'); z.tick_params(labelsize=7.5)
            z.set_title(f'{names[i][s_]}  ({L0[i, s_, 0]:.2f}, {L0[i, s_, 1]:.3f}) m', fontsize=8.5)
            if i: z.set_yticklabels([])
            else: z.set_ylabel(('Upper' if s_ else 'Lower')+' points\nΔy (mm)', fontsize=8.5)
            if row == 2: z.set_xlabel('Δx (mm)', fontsize=8.5)
            else: z.set_xticklabels([])
    cb = fig.colorbar(ScalarMappable(norm, cmap), ax=fig.axes[1:], shrink=.85, pad=.015, aspect=30)
    cb.set_label('PSQP iteration from the sweep restart (light = early)')
    fig.text(.125, .03, 'Lower panels: displacement of each control point from the NACA 0012 lattice, same scale in every panel, including the '
             'rigid incidence rotation about the pivot (0.125, 0.080) m.\nTrails start at the Oct 29 design (iteration 0). A remeshed chain restarts '
             'its shape variables at zero on its new grid; its control points are shown in absolute position.', fontsize=8, color=muted)
    save(fig, f"control_point_history_{r['dir'].name.replace('continuation_', '').replace('_20261031', '')}")

# ---------------------------------------------------------------- section shapes and the pinch
fig, (ax, zx) = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw=dict(width_ratios=[2.2, 1]))
shapes = [(naca, '#e6e6e6', 'NACA 0012 baseline', dict(fill=True)), (start['pts'], muted, 'Oct 29 final (sweep start)', dict(ls=':'))]
for name, r in runs.items():
    last = r['its'][-1]
    if r['kind'] == 'free':
        pre = max([q for q in r['its'] if q['tloc'] >= .002], key=lambda q: -q['cl'])
        shapes.append((pre['pts'], hue[r['xmax']], f'XMAX {r["xmax"]}, no limit: best with t ≥ 2 mm, {-pre["cl"]*c0:.3f} m', dict(lw=1)))
        shapes.append((last['pts'], hue[r['xmax']], f'XMAX {r["xmax"]}, no limit: last, {-last["cl"]*c0:.3f} m, '
                       f't_min {last["tloc"]*1e3:.2f} mm', dict(lw=1, ls=(0, (4, 2)))))
    else:
        shapes.append((last['pts'], hue[r['xmax']], f'XMAX {r["xmax"]}, {phase[r["kind"]]}: last, {-last["cl"]*c0:.3f} m',
                       dict(lw=style[r['kind']]['lw'])))
def thickness_along_chord(pts):
    xy = pts[:, :2]; ch = (xy[0]+xy[272])/2-xy[136]; L = np.sqrt(ch@ch)
    k = np.arange(1, 136); mid = (xy[k]+xy[272-k])/2-xy[136]
    return (mid@ch)/L**2, local_thickness_all(pts)
def local_thickness_all(pts):
    xy = pts[:, :2]; ch = (xy[0]+xy[272])/2-xy[136]; ch = ch/np.sqrt(ch@ch); k = np.arange(1, 136)
    return -(xy[272-k]-xy[k])@np.array([-ch[1], ch[0]])
for pts, col, lab, st in shapes:
    if st.get('fill'): ax.fill(pts[:, 0], pts[:, 1], color=col, lw=0, label=lab)
    else: ax.plot(pts[:, 0], pts[:, 1], c=col, lw=st.get('lw', 1), ls=st.get('ls', '-'), label=lab)
    xc, t = thickness_along_chord(pts)
    zx.plot(xc, t*1e3, c='#bdbdbd' if st.get('fill') else col, lw=st.get('lw', 1.2), ls=st.get('ls', '-'))
ax.axhline(0.05, c=muted, lw=.8, ls='--'); ax.text(-.12, .051, 'clearance limit 0.05 m', fontsize=8, color=muted, va='bottom')
for xl in (-.125, .625): ax.axvline(xl, c=ink, lw=1, ls='--')
ax.set(xlim=(-.135, .635), ylim=(.03, max(sh[0][:, 1].max() for sh in shapes)+.01), xlabel='x (m)', ylabel='Height above ground (m)')
ax.set_aspect('equal')
ax.set_title('Sections, true scale; dashed verticals are the box x limits', fontsize=10.5)
ax.legend(fontsize=7.5, frameon=False, loc='upper center', bbox_to_anchor=(.5, -.18), ncol=2)
xs = [thickness_along_chord(sh[0])[0] for sh in shapes]
zx.axvspan(min(x[124-1] for x in xs), max(x[15-1] for x in xs), color='#f2f2f2', zorder=0)
zx.axhline(2, c=ink, lw=1, ls='--'); zx.text(.02, 2.15, '2 mm floor', fontsize=8, color=ink)
zx.set(yscale='log', ylim=(.04, 80), xlim=(0, 1), xlabel='Chordwise position x/c (current chord)', ylabel='Local thickness (mm)')
zx.set_title('Local chord-normal thickness (shaded: constrained pairs k = 15–124)', fontsize=10.5)
for z in (ax, zx): z.grid(c=grid, lw=.6); z.set_axisbelow(True)
fig.tight_layout(); save(fig, 'section_comparison')

summary = {name: dict(kind=r['kind'], xmax=r['xmax'], offset=r['offset'], parent=r['parent'],
                      stages=[{k: v for k, v in s.items() if k in ('stage', 'drag_cap', 'best_clc_m', 'best_cdc_m', 'gain', 'drag_stop',
                               'x_margin_lead_m', 'x_margin_trail_m', 'chord_m', 'wall_min', 'stop', 'rejected', 'first')} for s in r['stages']],
                      iterates=[dict(i=i+r['offset'], stage=q['stage'], evaluation=q['evaluation'], clc_m=q['cl']*c0, cdc_m=q['cd']*c0,
                                     cap_cdc_m=q['cap']*c0, min_local_thickness_mm=q['tloc']*1e3, chord_m=q['actual_chord_m'],
                                     x_margin_lead_mm=(q['minimum_x_m']+.125)*1e3, x_margin_trail_mm=(.625-q['maximum_x_m'])*1e3,
                                     clearance_mm=q['minimum_clearance_m']*1e3, incidence_deg=q['chord_incidence_deg'])
                                for i, q in enumerate(r['its'])]) for name, r in runs.items()}
(out.parent/'continuation_summary.json').write_text(json.dumps(summary, indent=1)+'\n')
shutil.rmtree(scratch)
print('wrote', sorted(p.name for p in out.iterdir()))
