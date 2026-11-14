"""Plot a chain of optimization runs, each started from the previous one's last iterate.

Colour marks the grid (the original one, then each remesh). Wall time is cumulative over
the runs' evaluations; a later run's PSQP iteration 0 repeats its start design and is not
counted again. Writes under <out>/plots:
  chain_history   downforce over wall time; cost and flow solves of every evaluation;
                  distance to every limit; drag polar; recirculation area
  chain_sections  the section at every PSQP iterate, true scale, in the box
and <out>/chain_summary.json. Recirculation areas come from each run's
iterate_flows_*.json (plot_iterate_flows.py) where present.
  python scripts/plot_run_chain.py --out optimizations/run_NNN_<name>/summary optimizations/run_NNN_<name>/r00 ...
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

P = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('runs', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True)
p.add_argument('--compare', type=Path, default=P/'optimizations/continuation_summary_20261031/iterate_flows_remesh2_xmax1.0.json',
               help='iterate_flows json of an earlier chain, drawn behind the drag polar and recirculation panels')
a = p.parse_args()
out = a.out.resolve(); (out/'plots').mkdir(parents=True, exist_ok=True)

def iterates(run, H):
    """PSQP iterates: the evaluations behind the NIT lines of PSQP.out."""
    grads = list(dict.fromkeys(json.loads(l)['evaluation_count']-1 for l in (run/'gradients.jsonl').read_text().splitlines()))
    text = (run/'PSQP.out').read_text().split('EXIT FROM PSQP')[0]  # the exit summary repeats the last line
    L = [(int(v), int(g)) for v, g in re.findall(r'NIT=\s*\d+\s+NFV=\s*(\d+)\s+NFG=\s*(\d+)\s+F=', text)]
    F0 = float(re.search(r'NIT=\s*\d+\s+NFV=\s*\d+\s+NFG=\s*\d+\s+F=\s*(\S+)', text).group(1))
    # The start is the latest gradient evaluation with cl = F[0] (the derivative check also
    # solves the start design); the printed F differs from cl where a constraint is violated.
    starts = [k for k, e in enumerate(grads) if abs(H[e]['cl']-F0) < 1e-4][::-1]
    # Match by count: each PSQP function call is one evaluation, so an iterate is NFV-NFV0
    # evaluations after the start, and every iterate lands on a gradient evaluation.
    for k0 in starts:
        its = [grads[k0]+v-L[0][0] for v, g in L]
        if all(e in grads for e in its): return its
    # After failed solves function calls drift from the evaluations; then count gradients:
    # an iterate is the (NFG-NFG0)-th gradient evaluation after the start.
    for k0 in starts:
        if k0+L[-1][1]-L[0][1] < len(grads): return [grads[k0+g-L[0][1]] for v, g in L]
    raise ValueError(f'{run}: PSQP.out iterations do not match gradients.jsonl')

runs, t0, grids = [], 0., []
for r in a.runs:
    r = r.resolve()
    H = [json.loads(l) for l in (r/'history.jsonl').read_text().splitlines()]
    cmd = (r/'command.txt').read_text().split('\n')[0].split()
    grid = cmd[cmd.index('--baseline')+1] if '--baseline' in cmd else 'original'
    if grid not in grids: grids.append(grid)
    area = {}
    for f in sorted(r.glob('iterate_flows_*.json')):
        for q in json.loads(f.read_text()): area[int(Path(q['iterate']).name.split('_')[1])] = q['reversed_area_cm2']
    for h in H:
        h.setdefault('solves', 1); h.setdefault('stalled_solves', int(h['cfd_failed']))  # a direct solve, no walk
        h['t_min'] = (t0+h['elapsed_wall_time_s'])/60
    its = iterates(r, H)
    runs.append(dict(dir=r, H=H, its=its if not runs else its[1:], grid=grids.index(grid), area=area))
    t0 += H[-1]['elapsed_wall_time_s']
settings = json.loads((a.runs[0]/'settings.json').read_text())
c0, box = settings['reference_chord_m'], settings['box_m']
tmin, te_base, te_angle = settings['min_local_thickness_m'], settings['te_min_thickness_m'], settings['te_min_angle_deg']

ramp = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
cmap = LinearSegmentedColormap.from_list('early_to_late', ramp)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c = '#2b2b2b', '#6b6b6b', '#e4e4e4'
blue, orange, aqua = '#2a78d6', '#eb6834', '#1baf7a'
GRID_COLORS = [blue, aqua, '#0d366b', '#86b6ef']  # original grid, then each remesh
if len(grids) > len(GRID_COLORS):  # a long chain: one ordered ramp, dark to light, in run order
    GRID_COLORS = [plt.get_cmap('viridis')(.88*k/(len(grids)-1)) for k in range(len(grids))]
NCOL = 1 if len(grids) <= 5 else 2  # legend columns for the per-grid entries
grid_name = lambda g: 'original grid, warped' if g == 0 else f'remesh {g}, warped'
def save(fig, name):
    fig.savefig(out/'plots'/(name+'.png'), dpi=170, bbox_inches='tight')
    fig.savefig(out/'plots'/(name+'.svg'), bbox_inches='tight'); plt.close(fig)

# Chain-wide arrays: every evaluation, and the iterates
E = [dict(h, run=i, grid=r['grid'], iterate=h['evaluation'] in r['its']) for i, r in enumerate(runs) for h in r['H']]
for j, e in enumerate(E): e['n'] = j
I = [e for e in E if e['iterate']]
for k, e in enumerate(I): e['k'] = k
col = lambda e: GRID_COLORS[e['grid'] % len(GRID_COLORS)]
v = lambda key, S=E: np.array([s[key] for s in S])
down = lambda S: -c0*v('cl', S)

fig, axs = plt.subplots(2, 3, figsize=(16.5, 9.6)); axs = axs.ravel()
ax = axs[0]
for g in range(len(grids)):
    S = [e for e in I if e['grid'] == g]
    ax.plot(v('t_min', S), down(S), c=GRID_COLORS[g], lw=1.4, marker='o', ms=3.2, mec='white', mew=.4, label=grid_name(g))
bad = [e for e in E if e['cfd_failed']]
ax.scatter(v('t_min', bad), down(bad), c=orange, marker='X', s=26, lw=0, zorder=4, label='evaluation whose solve failed')
for i, r in enumerate(runs[1:], 1):
    ax.axvline(r['H'][0]['t_min'], c=muted, lw=.8, ls=':')
    if NCOL == 1: ax.text(r['H'][0]['t_min'], .02, f' run {i+1}', transform=ax.get_xaxis_transform(), fontsize=7.5, color=muted)
for i, r in enumerate(runs[1:], 1):
    if r['grid'] == runs[i-1]['grid'] or NCOL > 1: continue  # a long chain has no room for these labels
    # the same design, re-solved on the fresh grid
    prev = next(h for h in runs[i-1]['H'] if h['evaluation'] == runs[i-1]['its'][-1]); new = r['H'][0]
    ax.annotate(f"remesh: {-c0*prev['cl']:.3f} → {-c0*new['cl']:.3f}\n({new['cl']/prev['cl']-1:+.1%}), same design", (new['t_min'], -c0*new['cl']),
                xytext=(8, -26) if new['t_min'] < .6*E[-1]['t_min'] else (-14, 30), ha='left' if new['t_min'] < .6*E[-1]['t_min'] else 'right',
                textcoords='offset points', fontsize=7.5, color=ink, arrowprops=dict(arrowstyle='->', color=muted, lw=.7))
ax.annotate(f"{down([I[-1]])[0]:.3f}", (I[-1]['t_min'], down([I[-1]])[0]), xytext=(-4, 6), textcoords='offset points', ha='right', fontsize=8.5)
ax.set(xlabel='Wall time, cumulative (min)', ylabel='Downforce $-C_Lc$', title='(a) Progress: PSQP iterates')
ax.legend(fontsize=8, frameon=False, loc='upper left', ncol=NCOL)

ax = axs[1]
kinds = [('converged in one solve', blue, lambda e: not e['cfd_failed'] and not e['stalled_solves']),
         ('converged after a stalled solve', aqua, lambda e: not e['cfd_failed'] and e['stalled_solves']),
         ('failed', orange, lambda e: e['cfd_failed'])]
for lab, c, f in kinds:
    S = [e for e in E if f(e)]
    ax.bar(v('n', S), v('evaluation_wall_time_s', S), width=.8, color=c, lw=0, label=f'{lab} ({len(S)})')
ax.set(xlabel='Evaluation, cumulative', ylabel='Flow-solve wall time (s)', title='(b) Cost of each evaluation')
ax.legend(fontsize=8, frameon=False, loc='upper right')

ax = axs[2]
conv = v('solves')-v('stalled_solves')
ax.bar(v('n'), conv, width=.8, color=blue, lw=0, label='converged solves')
ax.bar(v('n'), v('stalled_solves'), bottom=conv, width=.8, color=orange, lw=0, label='stalled solves (cycle cap)')
ax.set(xlabel='Evaluation, cumulative', ylabel='Flow solves', title='(c) Solves per evaluation (walks through intermediate designs)')
ax.legend(fontsize=8, frameon=False, loc='upper left', ncol=NCOL)

ax = axs[3]
limits = [(f'gap above the {tmin*1e3:g} mm floor', aqua, '-', lambda h: (h['minimum_surface_gap_m']-tmin)*1e3),
          ('leading edge to the box', blue, '-', lambda h: (h['minimum_x_m']-box['xmin'])*1e3),
          ('trailing edge to the box', blue, (0, (4, 2)), lambda h: (box['xmax']-h['maximum_x_m'])*1e3),
          ('clearance above the box floor', orange, '-', lambda h: (h['minimum_clearance_m']-box['ymin'])*1e3),
          ('top of the section to the box top', orange, (0, (4, 2)), lambda h: (box['ymax']-h['maximum_height_m'])*1e3)]
if te_base: limits.append((f'trailing-edge base above the {te_base*1e3:g} mm floor', aqua, (0, (4, 2)), lambda h: (h['te_base_thickness_m']-te_base)*1e3))
if te_angle: limits.append((f'trailing-edge gap above the {te_angle:g}° taper', aqua, ':', lambda h: h['minimum_te_taper_margin_m']*1e3))
kk = v('k', I); margins = {}
for lab, c, ls, f in limits:
    y = np.array([f(e) for e in I]); margins[lab] = y
    ax.plot(kk, y, c=c, ls=ls, lw=1.8, label=lab)
for i, r in enumerate(runs[1:], 1):
    first = [e for e in I if e['run'] == i]
    if first: ax.axvline(first[0]['k']-.5, c=muted, lw=.8, ls=':')
ax.set_yscale('symlog', linthresh=1); ax.set_ylim(-.3, 1e4); ax.set_xlim(-.5, kk[-1]*1.25+1)
ax.axhline(0, c=ink, lw=1.2)
ax.set(xlabel='PSQP iteration, cumulative', ylabel='Distance to the limit (mm)', title='(d) Constraints')
ax.legend(fontsize=7.5, frameon=False, loc='upper left', ncol=2)

ax = axs[4]
ref = json.loads(a.compare.read_text()) if a.compare.exists() else []
if ref:
    ax.scatter([-q['cl_c_m'] for q in ref], [q['cd_c_m'] for q in ref], c='#c4c4c4', s=10, lw=0, label='Oct 31 remeshed chain (drag capped)')
for g in range(len(grids)):
    S = [e for e in I if e['grid'] == g]
    ax.plot(down(S), c0*v('cd', S), c=GRID_COLORS[g], lw=1.4, marker='o', ms=3.2, mec='white', mew=.4, label=grid_name(g))
xr = np.array([.4, 5.])
for r in (10, 20, 40):
    ax.plot(xr, xr/r, c=grid_c, lw=.9, zorder=0)
xm = max(down(I).max(), 3.4)*1.05; ym = max(c0*v('cd', I).max(), .33)*1.08
for r in (10, 20, 40):
    x = min(xm*.98, ym*.95*r); ax.text(x, x/r, f'$-L/D$ = {r}', fontsize=7.5, color=muted, ha='right', va='bottom')
ax.set(xlim=(.4, xm), ylim=(0, ym), xlabel='Downforce $-C_Lc$', ylabel='$C_Dc$', title=f'(e) Drag polar (cap $C_Dc$ = {c0*settings["drag_cap"]:g})')
ax.legend(fontsize=8, frameon=False, loc='upper left', ncol=NCOL)

ax = axs[5]
if ref:
    ax.scatter([-q['cl_c_m'] for q in ref], [q['reversed_area_cm2'] for q in ref], c='#c4c4c4', s=10, lw=0, label='Oct 31 remeshed chain')
for g in range(len(grids)):
    S = [(e, runs[e['run']]['area'].get(e['evaluation'])) for e in I if e['grid'] == g]
    S = [(e, A) for e, A in S if A is not None]
    if S: ax.plot(down([e for e, _ in S]), [A for _, A in S], c=GRID_COLORS[g], lw=1.4, marker='o', ms=3.2, mec='white', mew=.4, label=grid_name(g))
ax.set(xlabel='Downforce $-C_Lc$', ylabel='Area with u < 0 (cm²)', title='(f) Recirculation zone')
ax.legend(fontsize=8, frameon=False, loc='upper left', ncol=NCOL)
for ax in axs:
    ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True); ax.set_title(ax.get_title(), fontsize=10.5)
fig.suptitle(f'NACA 0012 to the box: {len(runs)} runs, {len(I)-1} PSQP iterations, {len(E)} evaluations, '
             f'{E[-1]["t_min"]:.0f} min on {settings["mpi_ranks"]} ranks (XMAX {settings["optimizer_options"]["XMAX"]}, '
             f'true gap ≥ {tmin*1e3:g} mm)', fontsize=12.5, y=.995)
fig.text(.5, -.01, '$C_Lc = L\'/q_\\infty$, $C_Dc = D\'/q_\\infty$. Each run restarts PSQP from the previous run\'s last iterate; a remeshed '
         'run starts from that design on a fresh grid, with the previous converged flow as its initial state.',
         ha='center', fontsize=8.5, color=muted, wrap=True)
fig.tight_layout()
to_pt = 72/axs[3].figure.dpi
items = sorted((axs[3].transData.transform((kk[-1], y[-1]))[1]*to_pt, lab, y[-1]) for lab, y in margins.items())
placed = []
for yp, *_ in items: placed.append(max(yp, placed[-1]+10) if placed else yp)
for (yp, lab, y), yn in zip(items, placed):
    axs[3].annotate(lab.split(' to ')[0].split(' above ')[0], (kk[-1], y), xytext=(6, yn-yp), textcoords='offset points',
                    va='center', fontsize=7.5, color=ink)
save(fig, 'chain_history')

# ---------------------------------------------------------------- sections
surf = lambda e: np.loadtxt(runs[e['run']]['dir']/f"iterates/eval_{e['evaluation']:03d}/surface_geometry.csv", delimiter=',', skiprows=1)
naca = np.loadtxt(runs[0]['dir']/'baseline/surface_geometry.csv', delimiter=',', skiprows=1)
norm = Normalize(0, len(I)-1)
fig, ax = plt.subplots(figsize=(13, 5.8))
ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=1.1))
ax.text(box['xmin'], box['ymax']+.003, 'Allowed surface box', ha='left', va='bottom', color=muted, fontsize=9)
ax.axhline(0, c=ink, lw=2); ax.text(box['xmax']-.005, .004, 'Slip ground', ha='right', va='bottom', color=muted, fontsize=9)
ax.fill(naca[:, 0], naca[:, 1], color='#e6e6e6', lw=0, label='NACA 0012 baseline', zorder=1)
for e in list({id(q): q for q in I[::max(1, len(I)//60)]+[I[-1]]}.values()):
    s = surf(e); last = e is I[-1]
    ax.plot(s[:, 0], s[:, 1], c=ramp[-1] if last else cmap(norm(e['k'])), lw=1.7 if last else .6, zorder=3 if last else 2,
            label=(f"iteration {e['k']}: $-C_Lc$ {-c0*e['cl']:.3f}, $C_Dc$ {c0*e['cd']:.3f}, gap {e['minimum_surface_gap_m']*1e3:.2f} mm, "
                   f"incidence {e['chord_incidence_deg']:.1f}°") if last else None)
ax.set(xlim=(box['xmin']-.02, box['xmax']+.02), ylim=(-.012, box['ymax']+.02), xlabel='x (m)', ylabel='Height above ground (m)')
ax.set_aspect('equal'); ax.legend(loc='upper left', bbox_to_anchor=(.02, .92), fontsize=8.5, frameon=False)
ax.set_title('Section at the PSQP iterates, true scale (light = early)', fontsize=10.5)
cb = fig.colorbar(ScalarMappable(norm, cmap), ax=ax, shrink=.75, pad=.01, aspect=25); cb.set_label('PSQP iteration, cumulative')
save(fig, 'chain_sections')

summary = dict(
    runs=[dict(run=r['dir'].name, grid=grids[r['grid']], evaluations=len(r['H']), iterations=len(r['its']),
               failed=sum(h['cfd_failed'] for h in r['H']), stalled_solves=sum(h['stalled_solves'] for h in r['H']),
               solves=sum(h['solves'] for h in r['H']), wall_min=r['H'][-1]['elapsed_wall_time_s']/60,
               start_cl_c=c0*r['H'][0]['cl'], end_cl_c=c0*min(h['cl'] for h in r['H'] if not h['cfd_failed']))
          for r in runs],
    iterates=[dict(k=e['k'], run=e['run'], evaluation=e['evaluation'], t_min=e['t_min'], cl_c=c0*e['cl'], cd_c=c0*e['cd'],
                   gap_mm=e['minimum_surface_gap_m']*1e3, **{lab: float(margins[lab][e['k']]) for lab in margins},
                   incidence_deg=e['chord_incidence_deg'], recirculation_cm2=runs[e['run']]['area'].get(e['evaluation']))
              for e in I])
(out/'chain_summary.json').write_text(json.dumps(summary, indent=1)+'\n')
print(json.dumps(summary['runs'], indent=1))
