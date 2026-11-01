"""Plot one optimization run evaluation by evaluation, failed solves included.

plot_continuation.py follows PSQP major iterates only; this run's story is in the
line-search trials and failed solves between them. Written for the November 1 clean
run from the NACA 0012 baseline (XMAX 1.0, true-gap thickness floor). Under <run>/plots:
  run_history          downforce over wall time; solve time, residual reached and
                       line-search step of every evaluation; distance to each limit;
                       drag polar against the October 31 continuation chains
  iterate_diagnostics  recirculation area and warped-grid quality at each saved iterate
  section_history      the section at every PSQP iterate, true scale, in the box
and <run>/run_summary.json. PSQP iterates are the gradient evaluations matched to the
NIT lines of PSQP.out (earlier gradients belong to the derivative check before PSQP).
A line-search step is the trial's distance from its iterate over that of the first
trial after the iterate (PSQP's full step, already limited by XMAX).
Run with the MACH-Aero env active (check_warped_cells needs its CGNS library), after
plot_iterate_flows.py for the recirculation panel:
  python scripts/plot_clean_run.py optimizations/clean_xmax1.0_gap2mm_20261101
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_warped_cells import coords, hex_volumes

P = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('run', type=Path)
p.add_argument('--compare', type=Path, default=P/'optimizations/continuation_summary_20261031/continuation_summary.json',
               help='continuation_summary.json whose iterates are drawn behind the drag polar')
a = p.parse_args()
run = a.run.resolve(); out = run/'plots'; out.mkdir(exist_ok=True)
settings = json.loads((run/'settings.json').read_text())
c0, box = settings['reference_chord_m'], settings['box_m']
cmd = (run/'command.txt').read_text().split()
tmin = float(cmd[cmd.index('--min-local-thickness')+1])
def find(d, key):
    if isinstance(d, dict):
        if key in d: return d[key]
        for v in d.values():
            r = find(v, key)
            if r is not None: return r
target = find(settings, 'L2Convergence') or 1e-10

H = [json.loads(l) for l in (run/'history.jsonl').read_text().splitlines()]
assert all(h['evaluation'] == i for i, h in enumerate(H))
grads = list(dict.fromkeys(json.loads(l)['evaluation_count']-1 for l in (run/'gradients.jsonl').read_text().splitlines()))
F = [float(f) for f in re.findall(r'NIT=\s*\d+\s+NFV=\s*\d+\s+NFG=\s*\d+\s+F=\s*(\S+)', (run/'PSQP.out').read_text())]
k0 = next(k for k in range(len(grads)) if len(grads[k:]) >= len(F)
          and np.allclose([H[e]['cl'] for e in grads[k:k+len(F)]], F, atol=1e-5))
its = grads[k0:k0+len(F)]

n = len(H); ev = np.arange(n); pi = np.array(its); k = np.arange(len(its))
failed = np.array([h['cfd_failed'] for h in H])
outcome = np.where(failed, 2, np.where([h['nk_fallback'] for h in H], 1, 0))
down = -c0*np.array([h['cl'] for h in H]); cdc = c0*np.array([h['cd'] for h in H])
minutes = np.array([h['elapsed_wall_time_s'] for h in H])/60
solve = np.array([h['evaluation_wall_time_s'] for h in H])
resid = np.array([h['residual_ratio_initial'] for h in H])
X = np.array([np.concatenate([np.ravel(h['design'][key]) for key in sorted(h['design'])]) for h in H])
step = np.full(n, np.nan)
for e in range(its[0]+1, n):
    b = max(i for i in its if i < e)
    step[e] = np.linalg.norm(X[e]-X[b])/np.linalg.norm(X[b+1]-X[b])
full = [e for e in range(its[0]+1, n) if e-1 in its]
wall = next((e for e in full if all(failed[f] for f in full if f >= e)), None)  # every full step failed from here
wall_it = its.index(wall-1) if wall else None

ramp = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
cmap = LinearSegmentedColormap.from_list('early_to_late', ramp)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid = '#2b2b2b', '#6b6b6b', '#e4e4e4'
blue, orange, aqua = '#2a78d6', '#eb6834', '#1baf7a'
OUTCOMES = [('converged', blue, 'o'), ('converged on the NK retry', aqua, 'D'), ('failed: cap hit, retry failed too', orange, 'X')]
def save(fig, name):
    fig.savefig(out/(name+'.png'), dpi=170, bbox_inches='tight')
    fig.savefig(out/(name+'.svg'), bbox_inches='tight'); plt.close(fig)
def by_outcome(ax, x, y, sel):
    for j, (lab, col, m) in enumerate(OUTCOMES):
        s = sel & (outcome == j)
        ax.scatter(x[s], y[s], c=col, marker=m, s=26 if m == 'X' else 18, lw=0, zorder=3, label=lab)
def end_labels(ax, items, spacing=11):
    """Direct labels at line ends, pushed apart vertically (in points) so they do not overlap."""
    to_pt = 72/ax.figure.dpi
    pts = sorted((ax.transData.transform((x, y))[1]*to_pt, lab, x, y) for lab, x, y in items)
    placed = []
    for yp, *_ in pts: placed.append(max(yp, placed[-1]+spacing) if placed else yp)
    for (yp, lab, x, y), yn in zip(pts, placed):
        ax.annotate(lab, (x, y), xytext=(6, yn-yp), textcoords='offset points', va='center', fontsize=8, color=ink)
def mark_wall(ax, x):
    if x is None: return
    ax.axvline(x, c=muted, lw=.9, ls=':', zorder=1)
    ax.text(x, 1, ' full steps fail from here', transform=ax.get_xaxis_transform(), fontsize=7.5, color=muted, va='top')

# ---------------------------------------------------------------- run history
fig, axs = plt.subplots(2, 3, figsize=(16.5, 9.4)); axs = axs.ravel()
ax = axs[0]
trial = np.setdiff1d(ev[its[0]+1:], pi)
ax.plot(minutes[pi], down[pi], c=blue, lw=1.2, zorder=2)
ax.scatter(minutes[pi], down[pi], c=blue, s=18, lw=0, zorder=3, label='PSQP iterate')
ok, bad = trial[~failed[trial]], trial[failed[trial]]
ax.scatter(minutes[ok], down[ok], facecolor='white', edgecolor=blue, s=18, lw=1, zorder=3, label='line-search trial, converged')
ax.scatter(minutes[bad], down[bad], c=orange, marker='X', s=30, lw=0, zorder=4, label='line-search trial, solve failed')
if wall: ax.axvspan(minutes[wall]-solve[wall]/60, minutes[-1], color='#f3f3f3', lw=0, zorder=0)
mark_wall(ax, minutes[wall]-solve[wall]/60 if wall else None)
ax.annotate(f'{down[pi[-1]]:.3f}', (minutes[pi[-1]], down[pi[-1]]), xytext=(-4, 8), textcoords='offset points', ha='right', fontsize=8.5, color=ink)
ax.set(xlabel='Wall time since the start (min)', ylabel='Downforce $-C_Lc$', title='(a) Progress: every evaluation')
ax.legend(fontsize=8, frameon=False, loc='lower right')

ax = axs[1]
for j, (lab, col, m) in enumerate(OUTCOMES):
    s = outcome == j; ax.bar(ev[s], solve[s], width=.8, color=col, lw=0, label=lab)
ax.text(.02, .97, f'{failed.sum()} of {n} solves failed, using {solve[failed].sum()/60:.0f} of the {minutes[-1]:.0f} min;\n'
        f'all solves {solve.sum()/60:.0f} min, gradients and setup the rest', transform=ax.transAxes, va='top', fontsize=8.5, color=ink)
ax.set_ylim(0, solve.max()*1.3)
ax.set(xlabel='Evaluation', ylabel='Flow-solve wall time (s)', title='(b) Cost of each evaluation')
ax.legend(fontsize=8, frameon=False, loc='upper left', bbox_to_anchor=(0, .86))

ax = axs[2]
by_outcome(ax, ev, resid, ev >= 0)
ax.axhline(target, c=ink, lw=1, ls='--'); ax.text(0, target, f' target {target:g}', va='bottom', fontsize=8, color=ink)
ax.set_yscale('log')
ax.set(xlabel='Evaluation', ylabel='Final residual / initial residual', title='(c) How far each solve converged')

ax = axs[3]
by_outcome(ax, ev, step, np.isfinite(step))
ax.set_yscale('log', base=2); ax.set_yticks([1, .5, .25], ['1 (full)', '1/2', '1/4']); ax.minorticks_off()
ax.set_ylim(.18, 1.4)
ax.set(xlabel='Evaluation', ylabel='Trial step / full step', title='(d) Line-search step of each trial')

ax = axs[4]
limits = [(f'gap above the {tmin*1e3:g} mm floor', aqua, '-', lambda h: (h['minimum_surface_gap_m']-tmin)*1e3),
          (f'clearance above {box["ymin"]*1e3:g} mm', orange, '-', lambda h: (h['minimum_clearance_m']-box['ymin'])*1e3),
          ('leading edge to the box', blue, '-', lambda h: (h['minimum_x_m']-box['xmin'])*1e3),
          ('trailing edge to the box', blue, (0, (4, 2)), lambda h: (box['xmax']-h['maximum_x_m'])*1e3)]
margins = {}
for lab, col, ls, f in limits:
    y = np.array([f(H[e]) for e in its]); margins[lab] = y
    ax.plot(k, y, c=col, ls=ls, lw=1.8, label=lab)
ax.set_yscale('symlog', linthresh=1); ax.set_ylim(-.3, 300); ax.set_xlim(-.5, len(its)+7)
ax.axhline(0, c=ink, lw=1.2); mark_wall(ax, wall_it)
ax.set(xlabel='PSQP iteration', ylabel='Distance to the limit (mm)', title='(e) Constraints: none active at the end')
ax.legend(fontsize=8, frameon=False, loc='lower left', bbox_to_anchor=(0, .1))

ax = axs[5]
if a.compare.exists():
    chains = json.loads(a.compare.read_text())
    pts = np.array([(-q['clc_m'], q['cdc_m']) for c in chains.values() for q in c['iterates']])
    ax.scatter(pts[:, 0], pts[:, 1], c='#c4c4c4', s=10, lw=0, zorder=1,
               label='Oct 31 continuation iterates (drag capped,\nchord-normal thickness pairs)')
ax.plot(down[pi], cdc[pi], c=blue, lw=1.2, marker='o', ms=3.5, mec='white', mew=.4, zorder=3, label='this run, PSQP iterates')
xr = np.array([.4, 3.5])
for r in (10, 20, 40):
    ax.plot(xr, xr/r, c=grid, lw=.9, zorder=0)
    x = min(3.45, .3*r); ax.text(x, x/r, f'$-L/D$ = {r}', fontsize=7.5, color=muted, ha='right', va='bottom')
ax.set(xlim=(.4, 3.5), ylim=(0, .34), xlabel='Downforce $-C_Lc$', ylabel='$C_Dc$', title=f'(f) Drag polar: drag is free (cap $C_Dc$ = {c0*settings["drag_cap"]:g})')
ax.legend(fontsize=8, frameon=False, loc='upper left')
for ax in axs:
    ax.grid(c=grid, lw=.6); ax.set_axisbelow(True); ax.set_title(ax.get_title(), fontsize=10.5)
fig.suptitle(f'{run.name}: {len(its)-1} PSQP iterations from the NACA 0012 baseline, {n} evaluations, '
             f'{minutes[-1]:.0f} min on 64 ranks', fontsize=12.5, y=.995)
fig.text(.5, -.01, '$C_Lc = L\'/q_\\infty$ and $C_Dc = D\'/q_\\infty$. A failed solve stopped at the cycle cap and failed again on the '
         'NK retry from the last converged state; PSQP then got a large penalty and halved the step along the same direction.',
         ha='center', fontsize=8.5, color=muted, wrap=True)
fig.tight_layout()
end_labels(axs[4], [(lab, k[-1], y[-1]) for lab, y in margins.items()])
save(fig, 'run_history')

# ---------------------------------------------------------------- iterate diagnostics
def grid_stats(path):
    s = {}
    for name, Xg, blank in coords(path):
        v = hex_volumes(Xg); bad = v*np.sign(np.median(v)) <= 0
        live = (blank == 1) | (blank == -1)
        s[name.split('.')[0]] = dict(inverted=int(bad.sum()), inverted_live=int((bad & live).sum()), min_live=float(np.abs(v[live]).min()))
    return s
base = grid_stats(run/'baseline/flow_000_vol.cgns')
G = [grid_stats(run/f'iterates/eval_{e:03d}/flow_000_vol.cgns') for e in its]
zones = {'airfoil_near': ('near-body grid', ramp[-1]), 'domain': ('background grid', orange)}
flows = sorted(run.glob('iterate_flows_*.json'))
area = {}
if flows:
    for r in json.loads(flows[0].read_text()): area[int(Path(r['iterate']).name.split('_')[1])] = r['reversed_area_cm2']

fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
ax = axs[0]
# plot_iterate_flows keeps a repeated design once (the start design is evaluations 1 and 5)
ax.plot(k, [area[max((q for q in area if q <= e), default=e)] if area else np.nan for e in its], c=orange, lw=1.8, marker='o', ms=3.5)
ax.set(xlabel='PSQP iteration', ylabel='Area with u < 0 (cm²)', title='(a) Recirculation zone (computed cells)')
ax = axs[1]
inv = np.array([g['domain']['inverted'] for g in G]); inv_live = sum(g[z]['inverted_live'] for g in G for z in g)
ax.bar(k, inv, width=.8, color=muted, lw=0)
ax.text(.02, .97, f'all inside the section (flooded, never computed);\ninverted computed cells, both grids: {inv_live}; '
        f'near-body grid: {sum(g["airfoil_near"]["inverted"] for g in G)}', transform=ax.transAxes, va='top', fontsize=8.5, color=ink)
ax.set_ylim(0, max(inv.max(), 1)*1.35)
ax.set(xlabel='PSQP iteration', ylabel='Inverted cells', title='(b) Inverted background-grid cells')
ax = axs[2]
for z, (lab, col) in zones.items():
    y = [g[z]['min_live']/base[z]['min_live'] for g in G]
    ax.plot(k, y, c=col, lw=1.8, label=lab)
ax.set_yscale('log'); ax.set_ylim(.05, 2); ax.set_xlim(-.5, len(its)+6)
ax.set(xlabel='PSQP iteration', ylabel='Smallest computed cell / NACA grid\'s', title='(c) Smallest computed cell volume')
ax.legend(fontsize=8, frameon=False, loc='lower left')
for ax in axs:
    mark_wall(ax, wall_it); ax.grid(c=grid, lw=.6); ax.set_axisbelow(True); ax.set_title(ax.get_title(), fontsize=10.5)
fig.suptitle('Flow and warped grid at each saved iterate: what changed where full steps began to fail', fontsize=12, y=1.0)
fig.tight_layout()
end_labels(axs[2], [(lab, k[-1], G[-1][z]['min_live']/base[z]['min_live']) for z, (lab, _) in zones.items()])
save(fig, 'iterate_diagnostics')

# ---------------------------------------------------------------- sections
surf = lambda d: np.loadtxt(d/'surface_geometry.csv', delimiter=',', skiprows=1)
naca = surf(run/'baseline')
norm = Normalize(0, len(its)-1)
fig, ax = plt.subplots(figsize=(13, 5.6))
ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=muted, lw=1.1))
ax.text(box['xmax']-.005, box['ymax']-.008, 'Allowed surface box', ha='right', va='top', color=muted, fontsize=9)
ax.axhline(0, c=ink, lw=2); ax.text(box['xmax']-.005, .004, 'Slip ground', ha='right', va='bottom', color=muted, fontsize=9)
ax.fill(naca[:, 0], naca[:, 1], color='#e6e6e6', lw=0, label='NACA 0012 baseline', zorder=1)
for i, e in enumerate(its):
    s = surf(run/f'iterates/eval_{e:03d}')
    last = i == len(its)-1
    ax.plot(s[:, 0], s[:, 1], c=ramp[-1] if last else cmap(norm(i)), lw=1.6 if last else .7, zorder=3 if last else 2,
            label=f'iteration {i} (evaluation {e}): $-C_Lc$ {down[e]:.3f}, $C_Dc$ {cdc[e]:.3f}, '
                  f'gap {H[e]["minimum_surface_gap_m"]*1e3:.2f} mm' if last else None)
ax.set(xlim=(box['xmin']-.02, box['xmax']+.02), ylim=(-.012, box['ymax']+.012), xlabel='x (m)', ylabel='Height above ground (m)')
ax.set_aspect('equal'); ax.legend(loc='upper left', bbox_to_anchor=(.02, .97), fontsize=8.5, frameon=False)
ax.set_title(f'Section at every PSQP iteration, true scale (light = early); chord incidence {H[its[0]]["chord_incidence_deg"]:.1f}° → '
             f'{H[its[-1]]["chord_incidence_deg"]:.1f}°', fontsize=10.5)
cb = fig.colorbar(ScalarMappable(norm, cmap), ax=ax, shrink=.75, pad=.01, aspect=25); cb.set_label('PSQP iteration')
save(fig, 'section_history')

last = H[its[-1]]
summary = dict(
    psqp_iterations=len(its)-1, evaluations=n, wall_time_min=float(minutes[-1]),
    psqp_iterate_evaluations=its, pre_optimization_evaluations=int(its[0]),
    failed_solves=ev[failed].tolist(), rescued_by_nk_retry=ev[outcome == 1].tolist(),
    failed_solve_minutes=float(solve[failed].sum()/60), all_solve_minutes=float(solve.sum()/60),
    full_steps_fail_from_evaluation=wall, full_steps_fail_from_psqp_iteration=wall_it,
    trial_steps={int(e): float(step[e]) for e in ev if np.isfinite(step[e])},
    start=dict(evaluation=its[0], cl_c=float(-down[its[0]]), cd_c=float(cdc[its[0]])),
    best=dict(evaluation=its[-1], cl_c=float(-down[its[-1]]), cd_c=float(cdc[its[-1]]),
              true_gap_mm=last['minimum_surface_gap_m']*1e3, clearance_mm=last['minimum_clearance_m']*1e3,
              leading_edge_margin_mm=(last['minimum_x_m']-box['xmin'])*1e3, trailing_edge_margin_mm=(box['xmax']-last['maximum_x_m'])*1e3,
              chord_incidence_deg=last['chord_incidence_deg'], recirculation_cm2=area.get(its[-1])),
    grid=[dict(iteration=i, evaluation=e, **{f'{z}_{q}': g[z][q]/base[z]['min_live'] if q == 'min_live' else g[z][q]
                                             for z in g for q in g[z]}) for i, (e, g) in enumerate(zip(its, G))])
(run/'run_summary.json').write_text(json.dumps(summary, indent=1)+'\n')
print(json.dumps({q: v for q, v in summary.items() if q not in ('trial_steps', 'grid', 'psqp_iterate_evaluations')}, indent=1))
