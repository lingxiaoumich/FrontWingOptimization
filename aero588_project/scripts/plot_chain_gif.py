"""Animate a chain of runs, one frame per PSQP iterate, in four panels of equal size:
  (a) top left      total-pressure coefficient C_pt
  (b) top right     pressure coefficient C_p: maize at +1, white at 0, blue at -10; the colour bar is linear in C_p
  (c) bottom left   the section on a background grid that is the same in every run (50 mm in x, 25 mm in y), and a
                    table of the section area and the geometric limits
  (d) bottom right  downforce -C_L c against drag C_D c, by mesh (0 = first grid, then each remesh); a drag or downforce
                    limit dashed; a star at the best feasible iterate so far (most downforce, or least loss in a
                    --max-clc run), since PSQP's last iterate can sit outside a limit; a table of the iterate's
                    -C_L c, C_D c, outlet C_pt loss (and Cp min with a Cp limit) against their constraints
Tables: columns current and constraint (N/A when not set); bold when the value is within 1 % of its limit (active),
bold and underlined when it is past the limit by more than 0.1 % (infeasible).
Curves are blue, one shade per round of the chain, darker each round (r00 lightest; the shades are fixed by the
round's place in the chain, the 6th and later share the darkest). Axis ranges, colour scales, grids and labels
are the same for every run, so frames of different runs can be flipped between. The field panels reach just below
y = 0: grey ground in a ground run, flow in a free-air run. Numbers in the tables have 2 decimals; the outlet loss
is in % of the free-stream dynamic pressure.
The iterates are those of <out>/chain_summary.json (plot_run_chain.py); each needs its flow file (a frame
without one says so). Writes <out>/plots/chain_iterates.gif and chain_iterates_last.png (the last frame, sharper).
Run with the MACH-Aero env:
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
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MultipleLocator
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
free_air = settings.get('ground_boundary', '').startswith('none')
drag_cap = c0*settings['drag_cap'] if settings.get('drag_cap', 10) < 5 else None   # C_D c; --drag-cap 10 means none
loss_cap, clc_max, cp_lim = settings.get('loss_cap'), settings.get('max_clc_m'), settings.get('cp_min')
amin, tmin = settings.get('min_section_area_m2') or 0., settings.get('min_local_thickness_m')
te_base, te_angle = settings.get('te_min_thickness_m'), settings.get('te_min_angle_deg')

its = json.loads((out/'chain_summary.json').read_text())['iterates']
for q in its:
    q['dir'] = runs[q['run']]/'iterates'/f"eval_{q['evaluation']:03d}"
    q['xy'] = np.loadtxt(q['dir']/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]
    q['res'] = json.loads((q['dir']/'result.json').read_text())
    q['loss'], q['cp_min'], q['feasible'] = (q['res'].get(k) for k in ('cpt_loss_outlet', 'cp_min_ks', 'feasible'))
N = len(its)
t = np.array([q['t_min'] for q in its])
down, cd, g = -np.array([q['cl_c'] for q in its]), np.array([q['cd_c'] for q in its]), np.array([q['run'] for q in its])

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c, ground_c = '#2b2b2b', '#6b6b6b', '#e4e4e4', '#d4d4d4'
SHADES = ['#93bff2', '#4387dd', '#1f62bd', '#154a8f', '#0d3366', '#081f45']  # round 1, 2, ... of a chain
shade = lambda k: SHADES[min(k, len(SHADES)-1)]
CPT = dict(cmap=LinearSegmentedColormap.from_list('cpt', ['#0d366b', '#2a78d6', '#86b6ef', '#f4f8fd']),
           norm=Normalize(-1.2, 1.02), levels=np.linspace(-1.2, 1.02, 75), ticks=None)
# Michigan blue at Cp = -10, a mid blue at -4, white at 0, maize at +1; on a linear scale, so 0 sits 10/11 up the bar
cp_at = lambda v: (v+10)/11
CP = dict(cmap=LinearSegmentedColormap.from_list('cp', [(0, '#00274c'), (cp_at(-4), '#2f6fbf'), (cp_at(0), 'white'), (1, '#ffcb05')]),
          norm=Normalize(-10., 1.), levels=np.r_[np.linspace(-10, 0, 51), np.linspace(0, 1, 11)[1:]],
          ticks=[-10, -8, -6, -4, -2, 0, 1])

# Fixed layout (inches): four axes of the same size, each with a colour-bar slot on its right
W, H, AW, AH = 15.5, 8.6, 5.7, 2.85
COLS, ROWS = (.75, 8.55), (4.90, 1.20)
rect = lambda x, y, w, h: [x/W, y/H, w/W, h/H]
FIELD = dict(xlim=(-.15, .95), ylim=(-.08, .47))   # 1.10 x 0.55 m: the panel's 2:1 at true scale
SECTION = dict(xlim=(-.08, .62), ylim=(-.115, .235))  # 0.70 x 0.35 m; the band below y = 0 holds the table
POLAR = dict(xlim=(0, .26), ylim=(0, 2.7))         # every run of the campaign fits (runs 002-039)

def ground(ax, xlim):
    if free_air: ax.axhline(0, c=muted, lw=.8, ls=':', zorder=4)
    else:
        ax.add_patch(Rectangle((xlim[0], -1), xlim[1]-xlim[0], 1, color=ground_c, lw=0, zorder=4))
        ax.axhline(0, c=ink, lw=1.6, zorder=4)

def scene(ax, xy, xlim):
    ax.fill(*xy.T, color='#9a9a9a', lw=.6, ec=ink, zorder=5)
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=ink, lw=.8, zorder=6))
    ground(ax, xlim); ax.set_aspect('equal')

def field(fig, ax, cax, zones, key, style, xy, label, title):
    if zones is None: ax.text(.5, .75, 'flow not saved for this iterate', transform=ax.transAxes, ha='center', color=muted)
    else:
        for tag in ['background', 'near']:
            f = zones[tag]
            ax.tricontourf(tri(f, tag == 'near'), np.clip(f[key], -50, 2).ravel(), extend='both', zorder=2, rasterized=True,
                           cmap=style['cmap'], norm=style['norm'], levels=style['levels'])
    scene(ax, xy, FIELD['xlim']); ax.set(**FIELD, xlabel='x (m)', ylabel='y (m)')
    ax.set_title(title, loc='left', fontsize=10.5)
    fig.colorbar(ScalarMappable(style['norm'], style['cmap']), cax=cax, label=label, extend='both', ticks=style['ticks'])

def shape_grid(ax):
    """The same background grid in every run, for comparing shapes: 50 mm in x, 25 mm in y."""
    ax.xaxis.set_major_locator(MultipleLocator(.1)); ax.xaxis.set_minor_locator(MultipleLocator(.05))
    ax.yaxis.set_major_locator(MultipleLocator(.05)); ax.yaxis.set_minor_locator(MultipleLocator(.025))
    ax.grid(which='both', color=grid_c, lw=.6); ax.set_axisbelow(True)

def status(v, lim, sense, ref):
    """'active' when the value is within 1 % of `ref` (the limit's scale) of its limit, 'violated' when it is past the
    limit by more than 0.1 % of `ref`; None otherwise or without a limit. sense: 'min' (value >= limit) or 'max'."""
    if v is None or lim is None: return None
    margin = v-lim if sense == 'min' else lim-v
    return 'violated' if margin < -.001*ref else 'active' if margin <= .01*ref else None

def table(ax, rows, corner, xy, fs=9):
    """A small table under the header ('', 'current', 'constraint'): rows (label, value, limit, state), numbers
    right-aligned. state 'active': value and limit bold; 'violated': bold, value underlined. The block's `corner`
    ('lower left' or 'lower right') sits at axes fraction `xy`. Returns its extent in axes fraction."""
    r = ax.figure.canvas.get_renderer(); bb = ax.get_window_extent(r)
    def width(s, bold=False):
        if not s: return 0.
        tx = ax.text(0, 0, s, fontsize=fs, fontweight='bold' if bold else 'normal', transform=ax.transAxes)
        w = tx.get_window_extent(r).width/bb.width; tx.remove(); return w
    cells = [('', 'current', 'constraint')]+[row[:3] for row in rows]
    wl = max(width(c[0]) for c in cells); wv = max(width(c[1], True) for c in cells); wn = max(width(c[2], True) for c in cells)
    gap, dy = .025, fs*1.55/72/AH
    tw, th = wl+gap+wv+gap+wn, dy*len(cells)
    x0 = xy[0]-tw if corner == 'lower right' else xy[0]; y0 = xy[1]
    xv, xn = x0+wl+gap+wv, x0+tw
    for i, c in enumerate(cells):
        y = y0+th-(i+.5)*dy
        kw = dict(fontsize=fs, transform=ax.transAxes, va='center', zorder=11)
        if i == 0:
            for txt, x in ((c[1], xv), (c[2], xn)): ax.text(x, y, txt, ha='right', color=muted, **kw)
            continue
        state = rows[i-1][3]; wt = 'bold' if state else 'normal'
        ax.text(x0, y, c[0], ha='left', **kw)
        tv = ax.text(xv, y, c[1], ha='right', fontweight=wt, **kw)
        ax.text(xn, y, c[2], ha='right', fontweight=wt, color=ink if state else muted, **kw)
        if state == 'violated':
            e = tv.get_window_extent(r).transformed(ax.transAxes.inverted())
            ax.add_artist(Line2D([e.x0, e.x1], [e.y0-.006, e.y0-.006], color=ink, lw=1, transform=ax.transAxes, zorder=11))
    return x0, y0, x0+tw, y0+th

def backdrop(ax, x0, y0, x1, y1, pad=.012):
    ax.add_patch(FancyBboxPatch((x0, y0), x1-x0, y1-y0, boxstyle=f'round,pad={pad}', transform=ax.transAxes,
                                fc='white', ec=grid_c, lw=.8, zorder=10))

f2 = lambda v: f'{round(v, 2)+0:.2f}'  # no -0.00
lim_text = lambda lim, sense: 'N/A' if lim is None else ('≥ ' if sense == 'min' else '≤ ')+f2(lim)
def row(label, v, lim, sense, ref, shown=None):
    return (label, shown if shown is not None else ('N/A' if v is None else f2(v)), lim_text(lim, sense), status(v, lim, sense, ref))

def frame(n, dpi=90, path=None):
    q = its[n]; res = q['res']
    flow = q['dir']/'flow_000_vol.cgns'
    zones = None
    if flow.exists(): zones = read_zones(flow); add_cpt(zones, mach)
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    axs = {k: fig.add_axes(rect(COLS[c], ROWS[r], AW, AH)) for k, (r, c) in zip('abcd', [(0, 0), (0, 1), (1, 0), (1, 1)])}
    cax = {k: fig.add_axes(rect(COLS[c]+AW+.12, ROWS[0], .14, AH)) for k, c in zip('ab', (0, 1))}
    field(fig, axs['a'], cax['a'], zones, 'Cpt', CPT, q['xy'], '$C_{pt}$', '(a) Total pressure $C_{pt}$')
    field(fig, axs['b'], cax['b'], zones, 'CoefPressure', CP, q['xy'], '$C_p$', '(b) Pressure $C_p$')

    ax = axs['c']
    ax.plot(*np.r_[its[0]['xy'], its[0]['xy'][:1]].T, c=muted, lw=.8, ls=(0, (3, 2)), zorder=3, label='start section')
    scene(ax, q['xy'], SECTION['xlim'])
    ax.set(**SECTION, xlabel='x (m)', ylabel='y (m)'); shape_grid(ax)
    ax.legend(loc='upper left', fontsize=8.5, frameon=False, handlelength=1.6, borderaxespad=.3)
    ax.set_title('(c) Section and limits', loc='left', fontsize=10.5)
    mm = lambda key: None if res.get(key) is None else 1e3*res[key]
    area = None if res.get('section_area_m2') is None else 1e4*res['section_area_m2']
    geo = [row('Area (cm²)', area, 1e4*amin if amin else None, 'min', 1e4*amin or 1),
           row('Thinnest (mm)', mm('minimum_surface_gap_m'), 1e3*tmin if tmin else None, 'min', 1e3*(tmin or 1e-3)),
           row('TE base (mm)', mm('te_base_thickness_m'), 1e3*te_base if te_base else None, 'min', 1e3*(te_base or 1e-3)),
           row('TE taper (mm)', mm('minimum_te_taper_margin_m'), 0. if te_angle else None, 'min', 10.),
           row('Box margin (mm)', mm('box_minimum_margin_m'), 0., 'min', 10.),
           row('Lowest y (mm)', mm('minimum_clearance_m'), 1e3*box['ymin'], 'min', 1e3*box['ymin'])]
    e1 = table(ax, geo[:3], 'lower left', (.025, .035), fs=8.5)
    e2 = table(ax, geo[3:], 'lower left', (e1[2]+.06, .035), fs=8.5)
    backdrop(ax, e1[0], e1[1], e2[2], max(e1[3], e2[3]))

    ax = axs['d']
    for k in range(len(runs)):
        s = g == k; done = s & (np.arange(N) <= n)
        if not s.any(): continue
        ax.plot(cd[s], down[s], c=shade(k), lw=1.2, alpha=.3)
        ax.plot(cd[done], down[done], c=shade(k), lw=1.8, marker='o', ms=2.5, label=f'mesh {k}')
    ok = [i for i in range(n+1) if its[i]['feasible']]
    if ok:  # the best feasible iterate so far
        by_loss = clc_max is not None and all(its[i]['loss'] is not None for i in ok)
        b = min(ok, key=lambda i: its[i]['loss']) if by_loss else max(ok, key=lambda i: down[i])
        ax.plot(cd[b], down[b], '*', ms=14, mfc=ink, mec='white', mew=.8, zorder=6, label='best feasible')
    ax.plot(cd[n], down[n], 'o', ms=9, mfc='white', mec=ink, mew=1.5, zorder=5)
    if drag_cap is not None:
        ax.axvline(drag_cap, c=ink, lw=1, ls='--')
        ax.text(drag_cap, .97*POLAR['ylim'][1], f' drag limit {drag_cap:.2f}', ha='left', va='top', fontsize=8.5, color=ink)
    if clc_max is not None:
        ax.axhline(-clc_max, c=ink, lw=1, ls='--')
        ax.text(.27*POLAR['xlim'][1], -clc_max, f'downforce limit {-clc_max:.2f}', ha='center', va='bottom', fontsize=8.5, color=ink)
    loss = None if free_air or q['loss'] is None else 100*q['loss']
    rows = [row('$-C_Lc$', down[n], None if clc_max is None else -clc_max, 'min', abs(clc_max or 1)),
            row('$C_Dc$', cd[n], drag_cap, 'max', drag_cap or 1),
            row('$C_{pt}$ loss (%)', loss, None if loss_cap is None else 100*loss_cap, 'max', 100*(loss_cap or .01))]
    if cp_lim is not None: rows.append(row('$C_p$ min', q['cp_min'], cp_lim, 'min', abs(cp_lim)))
    e = table(ax, rows, 'lower right', (.975, .035))
    backdrop(ax, *e)
    ax.set(**POLAR, xlabel='Drag $C_Dc$', ylabel='Downforce $-C_Lc$')
    ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
    ax.legend(loc='upper right', fontsize=8.5, frameon=False, borderaxespad=.3)
    ax.set_title('(d) Downforce vs. drag', loc='left', fontsize=10.5)

    fig.suptitle(f"{out.parent.name}{' (free air)' if free_air else ''}  ·  mesh {q['run']}  ·  iteration {n} of {N-1}  ·  {t[n]:.0f} min",
                 fontsize=12, y=.985)
    path = path or Path(tmp)/f'{n:04d}.png'; fig.savefig(path, dpi=dpi); plt.close(fig)
    return path

if __name__ == '__main__':
    (out/'plots').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp, Pool(a.jobs) as pool:
        frames = [Image.open(f).convert('RGB') for f in pool.map(frame, range(N))]
        frame(N-1, dpi=150, path=out/'plots/chain_iterates_last.png')
    frames[0].save(out/'plots/chain_iterates.gif', save_all=True, append_images=frames[1:],
                   duration=[a.ms]*(N-1)+[3000], loop=0, optimize=True)
    print(out/'plots/chain_iterates.gif', N, 'frames', frames[0].size)
