"""Animate a chain of runs, one frame per PSQP iterate: four panels of equal size and a constraint column.
  (a) top left      total-pressure coefficient C_pt
  (b) top right     pressure coefficient C_p: maize at +1, white at 0, blue at -10; the colour bar is linear in C_p.
                    In its top-right quadrant, an inset of the surface C_p against x (C_p up, fixed range;
                    wall cells of the near grid), lower surface blue, upper gold, a C_p limit dashed
  (c) bottom left   the section, its FFD control points over the chain's start lattice (dashed), and a maize comb of
                    the objective's surface sensitivity, on a background grid that is the same in every run (50 mm
                    in x, 25 mm in y); a control point outside the view is a triangle on the edge, pointing to it.
                    The comb: teeth along the outward normal toward -df/dn, the way a surface move improves the
                    objective (C_L, or the outlet loss in a --max-clc run), on an asinh scale fixed for the chain
                    (see comb). Only design-variable gradients are saved (gradients.jsonl, evaluation_count = the
                    iterate's evaluation + 1), so the surface sensitivity is that gradient taken back through the
                    FFD's normal moves (pyGeo, as the driver sets it up; see surface_sensitivity): the true
                    sensitivity as the FFD sees it, not ADflow's own dI/dXs. x shape variables are left out of it
  (d) bottom right  downforce -C_L c against drag C_D c, by mesh (0 = first grid, then each remesh); a drag or downforce
                    limit dashed; a star at the best feasible iterate so far (most downforce, or least loss in a
                    --max-clc run), since PSQP's last iterate can sit outside a limit
  (e) right         every constraint of the campaign's driver (run_slip_optimization.py), imposed or not, as a
                    horizontal bar of the share of its limit used: limit/value for a floor, value/limit for a cap,
                    1 - gap/box size for a side of the regulation box, 1 - margin/(local thickness floor) for the
                    trailing-edge taper. 1 is on the limit (the tick on the row); past it is outside. Colour by state:
                      grey   inactive: feasible and off the limit
                      blue   on the limit, within ACTIVE (1 %) of its value, or ACTIVE_MM (0.1 mm) for the box sides
                             and the taper, whose limits are positions or zero
                      maize  infeasible by the optimizer's own test (the 'feasible' of result.json): an aerodynamic
                             constraint past its limit by more than 1e-5 of it, a geometric one by more than 1e-6 c0
                             (0.5 um; area 1e-6 c0^2)
                    A row without a limit, the objective or a constraint the run did not impose (label and value
                    grey), has a hollow bar of its value over a fixed scale (NO_LIMIT_SCALE) and no tick.
                    The objective's row is outlined and bold, with an arrow along the bar in its limit column:
                    'max ->' for downforce, '<- min' for the outlet loss of a --max-clc run (downforce a floor there).
                    The surface C_p min of a run without a C_p limit is the lowest wall-cell C_p of the flow file
                    (within 0.02 of ADflow's KS cpmin, which is what a --cp-min run constrains and shows).
Curves are blue, one shade per round of the chain, darker each round (r00 lightest; the shades are fixed by the
round's place in the chain, the 6th and later share the darkest). Axis ranges, colour scales, grids and labels
are the same for every run, so frames of different runs can be flipped between. The field panels reach just below
y = 0: grey ground in a ground run, flow in a free-air run. Numbers have 2 decimals (C_D c 3, the outlet C_pt loss
4), more when an infeasible value would otherwise print as its limit. The loss is absolute, a fraction of the
free-stream dynamic pressure, as --loss-cap takes it.
The control points are the round's lattice (source/ffd.xyz), turned about the incidence pivot by the iterate's
incidence change and moved by its shape variables, as pyGeo does (checked to 1e-16 m against pyGeo, free
incidence included).
The iterates are those of <out>/chain_summary.json (plot_run_chain.py); each needs its flow file (a frame
without one says so). Writes <out>/plots/chain_iterates.gif and chain_iterates_last.png (the last frame, sharper);
with --png, only the one frame --frame (default the last) to that file, for trying out the layout.
Run with the MACH-Aero env:
  python scripts/plot_chain_gif.py --out optimizations/run_NNN_<name>/summary optimizations/run_NNN_<name>/r00 ...
"""
import argparse
import json
import os
import re
import tempfile
# One BLAS thread per process: the frames run in parallel (--jobs), often beside a solver using most cores
for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(var, '1')
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MultipleLocator
from PIL import Image

from flow_fields import add_cpt, read_zones, tri

p = argparse.ArgumentParser()
p.add_argument('runs', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True)
p.add_argument('--ms', type=int, default=250, help='ms per frame; the last one stays 3 s')
p.add_argument('--jobs', type=int, default=8, help='frames drawn at once (keep cores free for a running solver)')
p.add_argument('--png', type=Path, help='draw only the frame --frame to this file (no GIF)')
p.add_argument('--frame', type=int, default=-1, help='with --png: the iterate, counted from 0 (negative from the end)')
p.add_argument('--sens-fit', choices=['spread', 'projection'], default='spread',
               help='how the surface sensitivity comes from the design-variable gradient (see surface_sensitivity)')
a = p.parse_args()
out = a.out.resolve(); runs = [r.resolve() for r in a.runs]
settings = json.loads((runs[0]/'settings.json').read_text())
mach, box, c0 = settings['flow_derived']['mach'], settings['box_m'], settings['reference_chord_m']
free_air = settings.get('ground_boundary', '').startswith('none')
drag_cap = c0*settings['drag_cap'] if settings.get('drag_cap', 10) < 5 else None   # C_D c; --drag-cap 10 means none
loss_cap, clc_max, cp_lim = settings.get('loss_cap'), settings.get('max_clc_m'), settings.get('cp_min')
suction = bool(settings.get('min_suction'))  # --max-clc --min-suction: the suction peak is the objective, not the loss
amin, tmin = settings.get('min_section_area_m2') or None, settings.get('min_local_thickness_m') or None
te_base, te_angle = settings.get('te_min_thickness_m') or None, settings.get('te_min_angle_deg') or None
tmax = settings.get('max_section_thickness_min_m', .06*c0)  # always imposed: 0.06 c0

def lattice(run):
    """The round's FFD lattice (its first span plane; the two are tied), x and y: (stations, 2, 2), and the pivot
    of the incidence rotation (from the grid's baseline.json, or as the driver places it on a first grid)."""
    lines = (run/'source/ffd.xyz').read_text().split('\n')
    ns = int(lines[1].split()[0])
    raw = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
    lat = np.stack([raw[k].reshape(ns, 2, 2, order='F')[:, :, 0] for k in range(2)], -1)
    m = re.search(r'--baseline (\S+)', (run/'command.txt').read_text())
    if m: pivot = json.loads(Path(m.group(1)).read_text())['pivot_m']
    else:
        xy0 = np.loadtxt(run/'source/airfoil_coordinates.csv', delimiter=',', skiprows=1)
        pivot = [.25*c0, float(xy0[136, 1]+.25*c0*np.sin(np.deg2rad(3.)))]
    return lat, np.array(pivot)

def control_points(lat, pivot, design):
    """pyGeo's order: the incidence turns the lattice about the pivot, then the shape variables move its points."""
    th = np.deg2rad(-(design['incidence'][0]+3.))
    pts = (lat-pivot)@np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]).T+pivot
    for name, axis in ([('shape', 1)] if 'shape' in design else [('shape_x', 0), ('shape_y', 1)]):
        pts[..., axis] += np.reshape(design[name], lat.shape[:2])
    return pts

def normals(xy):
    """Outward unit normals of the closed surface polygon at its vertices, and each vertex's share of the arc."""
    t = np.roll(xy, -1, 0)-np.roll(xy, 1, 0)
    area = np.sum(xy[:, 0]*np.roll(xy[:, 1], -1)-np.roll(xy[:, 0], -1)*xy[:, 1])
    n = np.c_[t[:, 1], -t[:, 0]]*np.sign(area)
    e = np.linalg.norm(np.roll(xy, -1, 0)-xy, axis=1)
    return n/np.linalg.norm(n, axis=1)[:, None], .5*(e+np.roll(e, 1))

class SurfaceJacobian:
    """d(surface x, y)/d(design variables) of one round, with pyGeo set up as the driver sets it up: the round's
    lattice, the incidence turning it about the pivot, the shape variables moving its points after the turn. The
    surface is linear in the shape variables and their columns do not depend on the incidence, so they are
    differenced once; the incidence column (free-incidence runs only) at each iterate."""
    def __init__(self, run, pivot, groups):
        from pygeo import DVGeometry
        from pyspline import Curve
        span = settings['span_m']
        self.geo = geo = DVGeometry(str(run/'source/ffd.xyz'))
        geo.addRefAxis('incidence_axis', curve=Curve(k=2, t=[0., 0., 1., 1.], coef=np.array([
            [*pivot, -.005], [*pivot, span+.005]])), axis=np.array([0., 0., 1.]))
        def incidence(val, g): g.rot_z['incidence_axis'].coef[:] = -(val[0]+3.)
        geo.addGlobalDV('incidence', [-3.], incidence)
        idx = geo.getLocalIndex(0); ns = idx.shape[0]
        for name in groups:
            axis = {'shape': 1, 'shape_x': 0, 'shape_y': 1}[name]
            geo.addShapeFunctionDV(name, [{int(idx[i, j, k]): np.eye(3)[axis] for k in range(2)}
                                         for i in range(ns) for j in range(2)])
        xy0 = np.loadtxt(run/'source/airfoil_coordinates.csv', delimiter=',', skiprows=1)[:, :2]
        geo.addPointSet(np.c_[xy0, np.full(len(xy0), span/2)], 'surface')
        x0 = {'incidence': np.array([-3.]), **{name: np.zeros(2*ns) for name in groups}}
        base = self.at(x0); cols = []
        for name in groups:
            for d in range(2*ns):
                x = {k: v.copy() for k, v in x0.items()}; x[name][d] = 1e-3
                cols.append((self.at(x)-base)/1e-3)
        self.shape = np.stack(cols, -1)
    def at(self, x):
        self.geo.setDesignVars(x); return self.geo.update('surface')[:, :2].copy()
    def __call__(self, design):
        """With the shape variables after the turn, the incidence column depends on the incidence alone."""
        if settings.get('incidence_fixed'): return self.shape
        h, inc = 1e-3, np.array(design['incidence'], float)
        up, dn = (self.at({'incidence': inc+s*h}) for s in (1, -1))
        return np.concatenate([((up-dn)/(2*h))[..., None], self.shape], -1)

def surface_sensitivity(xy, jac, grad, groups):
    """The objective's sensitivity to an outward normal move of the surface, per unit arc length, from the saved
    gradient df/dx through the normal moves phi_d = n . dX/dx_d of the design variables:
      'projection'  the smallest-norm density s with int s phi_d ds = df/dx_d for every d: the true surface
                    sensitivity projected on the FFD's moves; reproduces the gradient exactly, but undoing the
                    overlap of neighbouring B-spline modes rings at the station spacing
      'spread'      s = sum_d phi_d (df/dx_d)/int phi_d^2 ds: each variable's gradient spread over its own move;
                    smooth (the true sensitivity seen through the FFD's modes), close to the gradient, not exact"""
    n, w = normals(xy)
    Jn = np.einsum('kc,kcd->kd', n, jac)
    g = np.concatenate(([] if settings.get('incidence_fixed') else [np.ravel(grad['incidence'])])
                       + [np.ravel(grad[name]) for name in groups])
    if a.sens_fit == 'spread': return Jn@(g/(w[:, None]*Jn**2).sum(0))
    A = np.sqrt(w)[:, None]*Jn
    return np.linalg.lstsq(A.T, g, rcond=None)[0]/np.sqrt(w)

its = json.loads((out/'chain_summary.json').read_text())['iterates']
lattices, jacobians, gradients = {}, {}, {}
for q in its:
    k = q['run']
    q['dir'] = runs[k]/'iterates'/f"eval_{q['evaluation']:03d}"
    q['xy'] = np.loadtxt(q['dir']/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]
    q['res'] = json.loads((q['dir']/'result.json').read_text())
    q['loss'], q['cp_min'], q['feasible'] = (q['res'].get(key) for key in ('cpt_loss_outlet', 'cp_min_ks', 'feasible'))
    # The fit uses the y moves (and a free incidence): an x move slides a thin section's surface mostly along
    # itself, so its gradient is mostly the grid's response, which the fit would blow up by 1/slope
    design = q['res']['design']; groups = [name for name in ('shape', 'shape_y') if name in design]
    if k not in lattices:
        lattices[k] = lattice(runs[k])
        jacobians[k] = SurfaceJacobian(runs[k], lattices[k][1], groups)
        path = runs[k]/'gradients.jsonl'
        gradients[k] = {r['evaluation_count']: r for r in map(json.loads, path.read_text().splitlines())} if path.exists() else {}
    q['ffd'] = control_points(*lattices[k], design)
    grad = gradients[k].get(q['evaluation']+1)  # written just after the iterate is saved, when the count is one on
    q['sens'] = None if grad is None or grad['failed'] else surface_sensitivity(q['xy'], jacobians[k](design), grad['objective'], groups)
N = len(its)
sens = np.abs(np.concatenate([q['sens'] for q in its if q['sens'] is not None] or [np.ones(1)]))
SENS_0, SENS_REF = np.median(sens), np.percentile(sens, 99)  # the chain's comb scale: asinh knee, full tooth
t = np.array([q['t_min'] for q in its])
down, cd, g = -np.array([q['cl_c'] for q in its]), np.array([q['cd_c'] for q in its]), np.array([q['run'] for q in its])

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid_c, ground_c = '#2b2b2b', '#6b6b6b', '#e4e4e4', '#d4d4d4'
SHADES = ['#93bff2', '#4387dd', '#1f62bd', '#154a8f', '#0d3366', '#081f45']  # round 1, 2, ... of a chain
shade = lambda k: SHADES[min(k, len(SHADES)-1)]
FFD_C = '#2a78d6'
MAIZE, MAIZE_EDGE = '#ffcb05', '#8a6a00'
COMB_M, COMB_TOOTH_M = .035, .008  # tooth length at the chain's 99th-percentile |sensitivity|; arc between teeth
COMB_LW = .6                       # teeth and the curve through their tips
STATE = {'inactive': ('#bdbdbd', 'inactive'), 'active': ('#00274c', 'on limit (within tol.)'),
         'infeasible': ('#ffcb05', 'infeasible')}  # grey, Michigan blue, maize
ACTIVE, ACTIVE_MM = .01, .1
BAR_MAX = 1.2  # the bar track of panel (e), in share of the limit
# A row without a limit (the objective, or a constraint the run does not impose) has a hollow bar of its value over
# a fixed scale, the same in every run: panel (d)'s axes, panel (b)'s colour bar; the outlet C_pt loss (absolute,
# a fraction of the free-stream dynamic pressure) and the area above the campaign's largest (0.032, run_023;
# 224.8 cm2, run_030)
NO_LIMIT_SCALE = dict(down=2.7, drag=.26, loss=.04, cp=-10., area=250.)
CPT = dict(cmap=LinearSegmentedColormap.from_list('cpt', ['#0d366b', '#2a78d6', '#86b6ef', '#f4f8fd']),
           norm=Normalize(-1.2, 1.02), levels=np.linspace(-1.2, 1.02, 75), ticks=None)
# Michigan blue at Cp = -10, a mid blue at -4, white at 0, maize at +1; on a linear scale, so 0 sits 10/11 up the bar
cp_at = lambda v: (v+10)/11
CP = dict(cmap=LinearSegmentedColormap.from_list('cp', [(0, '#00274c'), (cp_at(-4), '#2f6fbf'), (cp_at(0), 'white'), (1, '#ffcb05')]),
          norm=Normalize(-10., 1.), levels=np.r_[np.linspace(-10, 0, 51), np.linspace(0, 1, 11)[1:]],
          ticks=[-10, -8, -6, -4, -2, 0, 1])

# Fixed layout (inches): four axes of the same size, each with a colour-bar slot on its right, and the constraint
# column (e) on the right over the height of both rows
W, H, AW, AH = 20.6, 8.6, 5.7, 2.85
COLS, ROWS = (.75, 8.55), (4.90, 1.20)
XE, WE = 15.8, 4.55
HE = ROWS[0]+AH-ROWS[1]
rect = lambda x, y, w, h: [x/W, y/H, w/W, h/H]
FIELD = dict(xlim=(-.15, .95), ylim=(-.08, .47))   # 1.10 x 0.55 m: the panel's 2:1 at true scale
SECTION = dict(xlim=(-.08, .62), ylim=(-.05, .30))   # 0.70 x 0.35 m
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

# C_p up; every 8th and the last iterate of runs 001-039 span -12.6 (run_001) to 1.3
CP_INSET = dict(xlim=(-.03, .54), ylim=(-13., 1.5))
def cp_inset(ax, zones):
    """Surface C_p against x in the top-right quadrant of a field panel, on a white backing: the wall cells of the
    near grid (their order is the surface's: cells 0-135 the upper surface, 136-271 the lower, the rest the base)."""
    if zones is None: return
    ax.add_patch(FancyBboxPatch((.5, .5), .5, .5, boxstyle='round,pad=0,rounding_size=.02', transform=ax.transAxes,
                                fc='white', ec=grid_c, lw=.8, alpha=.93, zorder=9))
    ins = ax.inset_axes([.585, .615, .395, .315], zorder=10)
    f = zones['near']; x = .5*(f['CoordinateX'][:-1, 0]+f['CoordinateX'][1:, 0]); cp = f['CoefPressure'][:, 0]
    for s, c, name in ((slice(136, 272), '#00274c', 'lower'), (slice(0, 136), '#d9a300', 'upper')):
        ins.plot(x[s], cp[s], c=c, lw=1.2, label=name)
    if cp_lim is not None:
        ins.axhline(cp_lim, c=ink, lw=.8, ls='--')
        ins.text(CP_INSET['xlim'][1], cp_lim, 'limit ', ha='right', va='top', fontsize=6.5, color=ink)
    ins.set(**CP_INSET); ins.set_facecolor('none')
    ins.tick_params(labelsize=6.5, length=2, pad=1.5)
    ins.set_xlabel('x (m)', fontsize=7, labelpad=1); ins.set_ylabel('surface $C_p$', fontsize=7, labelpad=1)
    ins.yaxis.set_major_locator(MultipleLocator(4)); ins.xaxis.set_major_locator(MultipleLocator(.1))
    ins.grid(c=grid_c, lw=.5); ins.set_axisbelow(True)
    ins.legend(loc='lower center', ncol=2, fontsize=6.5, frameon=False, handlelength=1.2, borderaxespad=.2,
               columnspacing=1.2)

def shape_grid(ax):
    """The same background grid in every run, for comparing shapes: 50 mm in x, 25 mm in y."""
    ax.xaxis.set_major_locator(MultipleLocator(.1)); ax.xaxis.set_minor_locator(MultipleLocator(.05))
    ax.yaxis.set_major_locator(MultipleLocator(.05)); ax.yaxis.set_minor_locator(MultipleLocator(.025))
    ax.grid(which='both', color=grid_c, lw=.6); ax.set_axisbelow(True)

def ffd(ax, pts, start):
    """The lattice (upper and lower rows, one line per station) and its control points over the outline of the
    chain's start lattice (dashed: its rows and end stations); a point outside the view is a triangle on the
    nearest edge, pointing to it."""
    ax.add_collection(LineCollection([start[:, 0], start[:, 1], start[0], start[-1]], color=FFD_C, lw=.8, alpha=.55,
                                     ls=(0, (4, 2)), zorder=3))
    kw = dict(color=FFD_C, lw=.7, alpha=.6, zorder=7)
    ax.add_collection(LineCollection([pts[:, 0], pts[:, 1], *pts], **kw))
    (x0, x1), (y0, y1) = SECTION['xlim'], SECTION['ylim']
    p = pts.reshape(-1, 2); inside = (p[:, 0] >= x0) & (p[:, 0] <= x1) & (p[:, 1] >= y0) & (p[:, 1] <= y1)
    ax.plot(*p[inside].T, 'o', ms=4, mfc=FFD_C, mec='white', mew=.6, zorder=8)
    dx, dy = .006*(x1-x0), .012*(y1-y0)
    for x, y in p[~inside]:
        over = np.array([x0-x, x-x1, (y0-y)*2, (y-y1)*2])  # y counts double: the view is 2:1
        side = int(np.argmax(over))
        cx, cy = np.clip(x, x0+dx, x1-dx), np.clip(y, y0+dy, y1-dy)
        ax.plot(cx, cy, '<>v^'[side], ms=6, mfc=FFD_C, mec='white', mew=.6, zorder=8, clip_on=False)

def comb(ax, xy, s):
    """The sensitivity as a comb on the surface: a tooth along the outward normal in the direction of -s (the way
    a move improves the objective; outward when pushing the surface out helps), of length asinh(|s|/SENS_0) over
    asinh(SENS_REF/SENS_0) x COMB_M, at most 1.25 x COMB_M: linear below the chain's median |s|, logarithmic
    above, since the trailing edge's sensitivity is some 100 times the rest. The tips are joined into a curve.
    Teeth every COMB_TOOTH_M of arc; the curve through every vertex."""
    n, w = normals(xy)
    length = np.clip(np.arcsinh(-s/SENS_0)/np.arcsinh(SENS_REF/SENS_0), -1.25, 1.25)*COMB_M
    tips = xy+length[:, None]*n
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    teeth = np.unique(np.searchsorted(arc, np.arange(0, arc[-1], COMB_TOOTH_M)))
    ax.add_collection(LineCollection(np.stack([xy[teeth], tips[teeth]], 1), colors=MAIZE_EDGE, lw=COMB_LW, zorder=7.5))
    ax.plot(*np.r_[tips, tips[:1]].T, c=MAIZE, lw=COMB_LW, zorder=7.6)

def fmt(v, lim, nd, infeasible):
    """v to nd decimals (no -0.00); more, up to 6, when an infeasible value would print as its limit."""
    s = lambda v, d: f'{round(v, d)+0:.{d}f}'
    if infeasible and lim is not None:
        while nd < 6 and s(v, nd) == s(lim, nd): nd += 1
    return s(v, nd)

def constraints(q, zones, n):
    """Every constraint of the driver, as rows (group, label, value, limit, sense, share, decimals, state)."""
    res = q['res']
    mm = lambda key: None if res.get(key) is None else 1e3*res[key]
    def row(group, label, v, lim, sense, share, nd, ftol, atol, objective=None, scale=None):
        """sense 'min' (v >= lim) or 'max'; ftol the optimizer's feasibility tolerance, atol the 'on limit' one;
        objective 'maximize' or 'minimize' on the objective's row (never itself bounded); scale the value of a
        full bar when there is no limit (NO_LIMIT_SCALE)."""
        if lim is None or v is None:
            state = 'off'
            share = None if v is None or scale is None else BAR_MAX*v/scale
        else:
            margin = v-lim if sense == 'min' else lim-v
            state = 'infeasible' if margin < -ftol else 'active' if margin <= atol else 'inactive'
            share = share(v, lim, margin)
        return dict(group=group, label=label, v=v, lim=lim, sense=sense, nd=nd, state=state, objective=objective,
                    share=share)
    floor = lambda v, lim, m: lim/v if v > 0 else np.inf
    cap = lambda v, lim, m: v/lim
    side = lambda size: (lambda v, lim, m: 1-m/size)
    rel = lambda lim: 1e-5*abs(lim or 1)        # aerodynamic: vals <= 1 + 1e-5, each scaled by its limit
    L, A = 1e-3*c0, 1e-6*c0**2*1e4              # geometric: 1e-6 of c0 (mm), 1e-6 c0^2 (cm^2)
    act = lambda lim: ACTIVE*abs(lim or 1)
    loss = None if free_air else q['loss']
    cp_v = q['cp_min'] if cp_lim is not None or suction else None if zones is None else float(zones['near']['CoefPressure'][:, 0].min())
    area = None if res.get('section_area_m2') is None else 1e4*res['section_area_m2']
    lim_dn = None if clc_max is None else -clc_max
    lim_a = None if amin is None else 1e4*amin
    bx = {k: 1e3*v for k, v in box.items()}; bw, bh = bx['xmax']-bx['xmin'], bx['ymax']-bx['ymin']
    ag, ge, bo = 'Aerodynamic', 'Section', 'Regulation box'
    min_loss = clc_max is not None and not suction  # the outlet loss is the objective and downforce a floor; else downforce
    S = NO_LIMIT_SCALE
    return [
        row(ag, 'Downforce $-C_Lc$', down[n], lim_dn, 'min', floor, 2, rel(lim_dn), act(lim_dn),
            None if clc_max is not None else 'maximize', S['down']),
        row(ag, 'Drag $C_Dc$', cd[n], drag_cap, 'max', cap, 3, rel(drag_cap), act(drag_cap), scale=S['drag']),
        row(ag, 'Outlet $C_{pt}$ loss', loss, loss_cap, 'max', cap, 4, rel(loss_cap), act(loss_cap),
            'minimize' if min_loss else None, S['loss']),
        # with --min-suction the objective: a shorter bar (C_p min nearer 0) is a lower suction peak
        row(ag, 'Surface $C_p$ min', cp_v, cp_lim, 'min', cap, 2, rel(cp_lim), act(cp_lim),
            'minimize' if suction else None, scale=S['cp']),
        row(ge, 'Area (cm²)', area, lim_a, 'min', floor, 2, A, act(lim_a), scale=S['area']),
        row(ge, 'Max thickness (mm)', mm('maximum_section_thickness_m'), 1e3*tmax, 'min', floor, 2, L, act(1e3*tmax)),
        row(ge, 'Local thickness (mm)', mm('minimum_surface_gap_m'), tmin and 1e3*tmin, 'min', floor, 2, L, act(tmin and 1e3*tmin)),
        row(ge, 'TE base (mm)', mm('te_base_thickness_m'), te_base and 1e3*te_base, 'min', floor, 2, L, act(te_base and 1e3*te_base)),
        row(ge, 'TE taper margin (mm)', mm('minimum_te_taper_margin_m'), 0. if te_angle else None, 'min',
            side(1e3*(tmin or .002)), 2, L, ACTIVE_MM),
        row(bo, 'Front, min x (mm)', mm('minimum_x_m'), bx['xmin'], 'min', side(bw), 2, L, ACTIVE_MM),
        row(bo, 'Rear, max x (mm)', mm('maximum_x_m'), bx['xmax'], 'max', side(bw), 2, L, ACTIVE_MM),
        row(bo, 'Floor, min y (mm)', mm('minimum_clearance_m'), bx['ymin'], 'min', side(bh), 2, L, ACTIVE_MM),
        row(bo, 'Top, max y (mm)', mm('maximum_height_m'), bx['ymax'], 'max', side(bh), 2, L, ACTIVE_MM)]

def constraint_panel(ax, rows):
    """Drawn in inches from the panel's top left: label | bar of the share of the limit used | current | limit."""
    ax.set(xlim=(0, WE), ylim=(HE, 0)); ax.axis('off')
    ax.set_title('(e) Objective and constraints', loc='left', fontsize=10.5)
    XB0, XV, XL = 1.5, 3.62, WE           # bar start; right edges of current, limit
    U = 1.45/BAR_MAX                       # inches per unit of share
    xs = lambda u: XB0+U*min(max(u, 0.), BAR_MAX)
    width = lambda t: t.get_window_extent(ax.figure.canvas.get_renderer()).width/ax.figure.dpi
    def key(x, y, patch, label):
        ax.add_patch(patch); return x+width(ax.text(x, y, label, va='center', fontsize=8.5))+.22
    x, y = 0., .14  # the key: the three states, then a limitless bar and the objective's outline
    for state in ('inactive', 'active', 'infeasible'):
        c, label = STATE[state]
        x = key(x+.27, y, Rectangle((x, y-.065), .2, .13, fc=c, ec=ink, lw=.5), label)
    x, y = 0., .40
    x = key(x+.27, y, Rectangle((x, y-.065), .2, .13, fc='white', ec=ink, lw=.7), 'no limit (fixed scale)')
    key(x+.45, y, FancyBboxPatch((x, y-.09), .36, .18, boxstyle='round,pad=0,rounding_size=.05', fc='none', ec=ink,
                                 lw=.9), 'objective (arrow: improving direction)')
    y = .78; kw = dict(va='center', fontsize=8.5, color=muted)
    ax.text(XV, y, 'current', ha='right', **kw); ax.text(XL, y, 'limit', ha='right', **kw)
    ax.text(xs(1), y, '1', ha='center', **kw); ax.text(XB0, y, '0', ha='center', **kw)
    group = None
    for r in rows:
        if r['group'] != group:
            group = r['group']; y += .38
            ax.text(0, y, group.upper(), va='center', fontsize=8, fontweight='bold', color=muted)
        y += .335
        obj = r['objective']; off = r['state'] == 'off' and not obj
        txt = dict(va='center', color=muted if off else ink, fontweight='bold' if obj else 'normal', zorder=7)
        if obj:  # the objective: outlined across the row, bold
            ax.add_patch(FancyBboxPatch((-.05, y-.14), WE+.1, .28, boxstyle='round,pad=0,rounding_size=.06',
                                        fc='white', ec=ink, lw=.9, zorder=6, clip_on=False))
        ax.text(0, y, r['label'], fontsize=9, **txt)
        if r['share'] is not None:
            ax.add_patch(Rectangle((XB0, y-.1), U*BAR_MAX, .2, fc='#f2f2f2', lw=0, zorder=6.5))
            fill = 'white' if r['state'] == 'off' else STATE[r['state']][0]
            ax.add_patch(Rectangle((XB0, y-.1), xs(r['share'])-XB0, .2, fc=fill, ec=ink,
                                   lw=.7 if r['state'] == 'off' else .5, zorder=6.6))
            if r['share'] > BAR_MAX: ax.plot(xs(BAR_MAX)+.07, y, '>', ms=5, color=ink, clip_on=False, zorder=6.7)
        if r['state'] != 'off':  # the limit, at 1
            ax.plot([xs(1)]*2, [y-.15, y+.15], c=ink, lw=1., solid_capstyle='butt', zorder=6.8)
        v = 'N/A' if r['v'] is None else fmt(r['v'], r['lim'], r['nd'], r['state'] == 'infeasible')
        ax.text(XV, y, v, ha='right', fontsize=9, **txt)
        if obj: lim = 'max →' if obj == 'maximize' else '← min'  # along the bar: longer or shorter is better
        elif r['lim'] is None: lim = '–'
        else: lim = ('≥ ' if r['sense'] == 'min' else '≤ ')+fmt(r['lim'], None, r['nd'], False)
        ax.text(XL, y, lim, ha='right', fontsize=10 if obj else 9, **txt)
    S = NO_LIMIT_SCALE
    ax.text(0, y+.34, 'Bar: share of the limit used, 1 at the limit: limit/value for a floor, value/limit for a cap,\n'
            '1 − gap/box size for a box side, 1 − margin/local-thickness floor for the TE taper.\n'
            f"No limit: value over a full bar of $-C_Lc$ {S['down']}, $C_Dc$ {S['drag']}, $C_{{pt}}$ loss {S['loss']}, "
            f"$C_p$ {S['cp']:.0f}, area {S['area']:.0f} cm².", fontsize=7.5, color=muted, va='top', linespacing=1.4)

def frame(n, dpi=90, path=None):
    q = its[n]
    flow = q['dir']/'flow_000_vol.cgns'
    zones = None
    if flow.exists(): zones = read_zones(flow); add_cpt(zones, mach)
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    axs = {k: fig.add_axes(rect(COLS[c], ROWS[r], AW, AH)) for k, (r, c) in zip('abcd', [(0, 0), (0, 1), (1, 0), (1, 1)])}
    cax = {k: fig.add_axes(rect(COLS[c]+AW+.12, ROWS[0], .14, AH)) for k, c in zip('ab', (0, 1))}
    field(fig, axs['a'], cax['a'], zones, 'Cpt', CPT, q['xy'], '$C_{pt}$', '(a) Total pressure $C_{pt}$')
    field(fig, axs['b'], cax['b'], zones, 'CoefPressure', CP, q['xy'], '$C_p$', '(b) Pressure $C_p$')
    cp_inset(axs['b'], zones)

    ax = axs['c']
    ax.plot(*np.r_[its[0]['xy'], its[0]['xy'][:1]].T, c=muted, lw=.8, ls=(0, (3, 2)), zorder=3, label='start section')
    scene(ax, q['xy'], SECTION['xlim'])
    ffd(ax, q['ffd'], its[0]['ffd'])
    if q['sens'] is not None: comb(ax, q['xy'], q['sens'])
    ax.plot([], [], c=FFD_C, lw=.8, alpha=.55, ls=(0, (4, 2)), label='start FFD box')
    ax.plot([], [], 'o-', c=FFD_C, lw=.7, ms=4, mec='white', mew=.6, label='FFD control points')
    ax.plot([], [], c=MAIZE, lw=1.5, label='objective sensitivity $-\\partial f/\\partial n$ (asinh)' if q['sens'] is not None else 'no gradient saved')
    ax.set(**SECTION, xlabel='x (m)', ylabel='y (m)'); shape_grid(ax)
    ax.legend(loc='upper left', fontsize=8.5, frameon=False, handlelength=1.6, borderaxespad=.3)
    ax.set_title('(c) Section, FFD control points and sensitivity', loc='left', fontsize=10.5)

    ax = axs['d']
    for k in range(len(runs)):
        s = g == k; done = s & (np.arange(N) <= n)
        if not s.any(): continue
        ax.plot(cd[s], down[s], c=shade(k), lw=1.2, alpha=.3)
        ax.plot(cd[done], down[done], c=shade(k), lw=1.8, marker='o', ms=2.5, label=f'mesh {k}')
    ok = [i for i in range(n+1) if its[i]['feasible']]
    if ok:  # the best feasible iterate so far
        by_loss = clc_max is not None and not suction and all(its[i]['loss'] is not None for i in ok)
        by_cp = suction and all(its[i]['cp_min'] is not None for i in ok)
        b = (max(ok, key=lambda i: its[i]['cp_min']) if by_cp else min(ok, key=lambda i: its[i]['loss']) if by_loss
             else max(ok, key=lambda i: down[i]))
        ax.plot(cd[b], down[b], '*', ms=14, mfc=ink, mec='white', mew=.8, zorder=6, label='best feasible')
    ax.plot(cd[n], down[n], 'o', ms=9, mfc='white', mec=ink, mew=1.5, zorder=5)
    if drag_cap is not None:
        ax.axvline(drag_cap, c=ink, lw=1, ls='--')
        ax.text(drag_cap, .97*POLAR['ylim'][1], f' drag limit {drag_cap:.2f}', ha='left', va='top', fontsize=8.5, color=ink)
    if clc_max is not None:
        ax.axhline(-clc_max, c=ink, lw=1, ls='--')
        ax.text(.27*POLAR['xlim'][1], -clc_max, f'downforce limit {-clc_max:.2f}', ha='center', va='bottom', fontsize=8.5, color=ink)
    ax.set(**POLAR, xlabel='Drag $C_Dc$', ylabel='Downforce $-C_Lc$')
    ax.grid(c=grid_c, lw=.6); ax.set_axisbelow(True)
    ax.legend(loc='upper right', fontsize=8.5, frameon=False, borderaxespad=.3)
    ax.set_title('(d) Downforce vs. drag', loc='left', fontsize=10.5)

    constraint_panel(fig.add_axes(rect(XE, ROWS[1], WE, HE)), constraints(q, zones, n))

    fig.suptitle(f"{out.parent.name}{' (free air)' if free_air else ''}  ·  mesh {q['run']}  ·  iteration {n} of {N-1}  ·  {t[n]:.0f} min",
                 fontsize=12, y=.985)
    path = path or Path(tmp)/f'{n:04d}.png'; fig.savefig(path, dpi=dpi); plt.close(fig)
    return path

if __name__ == '__main__':
    if a.png:
        print(frame(a.frame % N, dpi=150, path=a.png)); raise SystemExit
    (out/'plots').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp, Pool(a.jobs) as pool:
        frames = [Image.open(f).convert('RGB') for f in pool.map(frame, range(N))]
        frame(N-1, dpi=150, path=out/'plots/chain_iterates_last.png')
    frames[0].save(out/'plots/chain_iterates.gif', save_all=True, append_images=frames[1:],
                   duration=[a.ms]*(N-1)+[3000], loop=0, optimize=True)
    print(out/'plots/chain_iterates.gif', N, 'frames', frames[0].size)
