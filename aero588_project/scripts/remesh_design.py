"""Build a fresh overset grid around an optimized design and re-base its parameterization.

The warped grid fails (negative volumes near the trailing edge) once the design has
moved far from the NACA baseline. This regrows the near-field block with pyHyp from
the design's actual 289-vertex surface, with the settings of the original case
(cases/naca0012_c050_h005_u200_am3/run_case.py: 81 layers, 2 um first cell, 25 mm
march), and rebuilds the identical Cartesian background. The vertex order is kept, so
the local-thickness pairs and chord definition in run_slip_optimization.py still apply.

The FFD lattice is written in its deformed position at the design, so starting the
optimizer with --baseline <out>/baseline.json and all-zero shape variables (incidence
-3, i.e. no rotation) reproduces the design; later moves are relative to it. The
incidence pivot stays the original one. Run serially with the MACH-Aero environment:
  python scripts/remesh_design.py --design D.json --out DIR [--baseline B.json]
(--baseline when D belongs to an earlier remeshed parameterization).
"""
import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

P = Path(__file__).resolve().parents[1]
h = argparse.ArgumentParser()
h.add_argument('--design', type=Path, required=True)
h.add_argument('--out', type=Path, required=True)
h.add_argument('--baseline', type=Path)
h.add_argument('--hyp-options', default='{}', help='JSON overrides of the pyHyp marching options')
# Grid family for a lower box floor (see `mesh` below); each default is that of --baseline, else the original case's.
h.add_argument('--march-dist', type=float, help='m; near-grid depth, pyHyp marchDist (original 0.025)')
h.add_argument('--ground-spacing', type=float, help='m; background cell height from the ground up to --fine-height, '
               'then growing 10 %% a cell to the original 2.5 mm (original 0.0025: uniform)')
h.add_argument('--fine-height', type=float, help='m; top of the --ground-spacing layer (original 0)')
h.add_argument('--ground-gap', type=float, help='m; lowest height of the near grid: its grid lines that reach lower are '
               'shortened about their wall point, smoothly from line to line (original 0: off)')
h.add_argument('--free-air', action='store_true', default=None,
               help='no ground: the background reaches 10 chords below y = 0 as well, to a far field there, with 2.5 mm '
               'cells from y = -0.15 m (the ground runs use it with run_slip_optimization.py --no-ground)')
own = h.parse_args()
MESH0 = {'march_dist_m': .025, 'ground_spacing_m': .0025, 'fine_height_m': 0., 'ground_gap_m': 0., 'free_air': False}
previous = json.loads(own.baseline.read_text()).get('mesh', {}) if own.baseline else {}
mesh = {k: v if v is not None else previous.get(k, MESH0[k])
        for k, v in zip(MESH0, [own.march_dist, own.ground_spacing, own.fine_height, own.ground_gap, own.free_air])}
if mesh['free_air'] and (mesh['ground_spacing_m'] < .0025-1e-12 or mesh['ground_gap_m'] > 0):
    h.error('--free-air has no ground: it takes neither --ground-spacing nor --ground-gap')
out = Path(os.path.abspath(own.out)); out.mkdir(parents=True, exist_ok=True)  # not resolved: keep a linked run folder's path

# The optimizer's own geometry setup, in a scratch root (it writes setup files there).
scratch = Path(tempfile.mkdtemp())
shutil.copytree(P/'setups/slip_ground_xy5_cd005/source', scratch/'source')
if own.baseline: shutil.copy(own.baseline.parent/'airfoil_coordinates.csv', scratch/'source')
sys.path.insert(0, str(scratch/'source'))
raw = json.loads(own.design.read_text())
directions = 'y' if 'shape' in raw else 'xy'  # 'shape' holds the moves of --shape-directions y
nstation = len(raw['shape' if directions == 'y' else 'shape_x'])//2  # upper/lower pair per station
sys.argv = ['x', '--root', str(scratch), '--ffd-stations', str(nstation), '--shape-directions', directions] + \
           (['--baseline', os.path.abspath(own.baseline)] if own.baseline else [])
S = {}
exec(compile((P/'scripts/run_slip_optimization.py').read_text().split('options={')[0], 'setup', 'exec'), S)
geo, pivot, span, c0 = S['geo'], S['pivot'], S['span'], S['c0']
design = S['load_design'](own.design)
geo.setDesignVars(design)
pts = geo.update('surface_check')

# Deformed lattice: rotate the undeformed one about the pivot, then add the shape moves (pyGeo's order).
lines = (scratch/'source/ffd.xyz').read_text().split('\n')
dims = [int(v) for v in lines[1].split()]
vals = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
L0 = np.stack([v.reshape(dims, order='F') for v in vals], -1)  # (station, lower/upper, span, xyz)
t = np.deg2rad(-(design['incidence'][0]+3.))
R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
L = L0.copy()
L[..., :2] = (L0[..., :2]-pivot[:2])@R.T+pivot[:2]
for name, ax in [('shape_x', 0), ('shape_y', 1), ('shape', 1)]:
    if name in design: L[..., ax] += np.reshape(design[name], dims[:2])[:, :, None]
with (out/'ffd.xyz').open('w') as f:
    f.write(f'1\n{dims[0]} {dims[1]} {dims[2]}\n')
    for axis in range(3): np.savetxt(f, L[..., axis].ravel(order='F'), fmt='%.16e')
# The next run rebuilds the surface from this lattice; a strongly deformed lattice can lose
# vertices (they no longer project inside it), which gives a self-crossing start shape.
from pygeo import DVGeometry
check = DVGeometry(str(out/'ffd.xyz')); check.addPointSet(pts, 'p')
embed_err = float(np.abs(check.update('p')-pts).max())
if embed_err > 1e-9:
    shutil.rmtree(scratch)
    raise SystemExit(f'the design surface does not re-embed in the deformed lattice ({embed_err:.1e} m)')

xy = pts[:, :2]
np.savetxt(out/'airfoil_coordinates.csv', xy, delimiter=',', header='x_m,y_m', comments='')
xyz = np.empty((len(xy), 2, 1, 3)); xyz[..., 0] = xy[:, 0, None, None]; xyz[..., 1] = xy[:, 1, None, None]
xyz[:, 0, :, 2], xyz[:, 1, :, 2] = 0, span
with (out/'airfoil_surface.xyz').open('w') as f:
    f.write(f'1\n{len(xy)} 2 1\n')
    for axis in range(3): np.savetxt(f, xyz[..., axis].ravel(order='F'), fmt='%.15e')

from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import Grid, Block, Boco, readGrid, combineGrids
hyp_options = {
    'unattachedEdgesAreSymmetry': False, 'outerFaceBC': 'overset', 'autoConnect': True,
    'BC': {1: {'jLow': 'zSymm', 'jHigh': 'zSymm'}}, 'families': 'airfoil',
    'N': 81, 's0': 2e-6, 'marchDist': mesh['march_dist_m'],
    'ps0': -1., 'pGridRatio': -1., 'cMax': 3., 'epsE': 1., 'epsI': 2.,
    'theta': 3., 'volCoef': .25, 'volBlend': .0001, 'volSmoothIter': 100}
hyp_options.update(json.loads(own.hyp_options))
hyp = pyHyp(options={'inputFile': str(out/'airfoil_surface.xyz'), **hyp_options})
hyp.run()
hyp.writeCGNS(str(out/'airfoil_near.cgns'))

def squash(coords, gap, width=12):
    """Lift the near grid off the ground where it reaches below `gap` (index 0 around the section, 1 span,
    2 layers): on line i the downward offsets from its wall point are scaled by s_i, the factor that puts
    the line's lowest point at `gap`, opened (cyclic min filter, then mean filter, both +-width lines),
    which never exceeds a line's own need. Only the vertical part shrinks, so lines that point forward or
    back (nose, trailing edge) keep their shape while the downward ones under the wing are compressed."""
    wall = coords[:, :, :1]
    yw, low = wall[:, 0, 0, 1], coords[..., 1].min(axis=(1, 2))
    need = np.where(low < gap, (yw-gap)/np.maximum(yw-low, 1e-30), 1.)[:-1]  # the last line is the first one
    win = lambda v, f: np.array([f(np.take(v, range(i-width, i+width+1), mode='wrap')) for i in range(len(v))])
    s = win(win(need, np.min), np.mean)
    s = np.r_[s, s[0]][:, None, None]
    out = coords.copy(); dy = coords[..., 1]-wall[..., 1]
    out[..., 1] = wall[..., 1]+np.where(dy < 0, s*dy, dy)
    return out, float(s.min())

def extend(start, stop, step, ratio=1.08):
    values = [start]; direction = np.sign(stop-start)
    while abs(stop-values[-1]) > step:
        values.append(values[-1]+direction*step); step *= ratio
    values.append(stop)
    return np.array(values)
chord = c0  # background spacing of the original case, which this grid must match
xx = np.r_[extend(-.15, -10*chord, .005)[::-1][:-1], np.linspace(-.15, .85, 201), extend(.85, 15*chord, .005)[1:]]
def ground_layer(spacing, top, coarse=.0025, ratio=1.1, end=.35):
    """y from the ground to `end`: `spacing` up to `top`, growing by `ratio` to `coarse`, then uniform."""
    if spacing >= coarse-1e-12: return np.linspace(0, end, int(round(end/coarse))+1)
    y = list(np.arange(0, top+spacing/2, spacing)); s = spacing
    while s*ratio < coarse: s *= ratio; y.append(y[-1]+s)
    return np.r_[y[:-1], np.linspace(y[-1], end, int(round((end-y[-1])/coarse))+1)]
yy = np.r_[ground_layer(mesh['ground_spacing_m'], mesh['fine_height_m']), extend(.35, 10*chord, .0025)[1:]]
if mesh['free_air']:  # the same lines above y = 0, then 2.5 mm down to -0.15 m, growing 8 % a cell to -10 chords
    yy = np.r_[extend(-.15, -10*chord, .0025)[::-1][:-1], np.linspace(-.15, 0, 61)[:-1], yy]
coordinates = np.stack(np.meshgrid(xx, yy, [0, span], indexing='ij'), axis=-1)
bdims = list(coordinates.shape[:3])
block = Block('background', bdims, coordinates)
def boundary(name, kind, axis, index, family):
    point_range = np.array([[1, bdims[0]], [1, bdims[1]], [1, bdims[2]]]); point_range[axis] = index
    block.addBoco(Boco(name, kind, point_range, family))
boundary('inlet', 'bcfarfield', 0, 1, 'farfield'); boundary('outlet', 'bcfarfield', 0, bdims[0], 'farfield')
boundary('top', 'bcfarfield', 1, bdims[1], 'farfield')
if mesh['free_air']: boundary('bottom', 'bcfarfield', 1, 1, 'farfield')
else: boundary('ground', 'bcsymmetryplane', 1, 1, 'ground')
boundary('span0', 'bcsymmetryplane', 2, 1, 'symmetry'); boundary('span1', 'bcsymmetryplane', 2, 2, 'symmetry')
background = Grid(); background.blocks = [block]
near = readGrid(str(out/'airfoil_near.cgns'))
squashed = 1.
if mesh['ground_gap_m'] > 0:
    near.blocks[0].coords, squashed = squash(near.blocks[0].coords, mesh['ground_gap_m'])
    near.writeToCGNS(str(out/'airfoil_near.cgns'))
grid = combineGrids([near, background])
if not mesh['free_air'] and grid.blocks[0].coords[..., 1].min() <= 0:  # the near grid must stay above the slip ground
    shutil.rmtree(scratch)
    raise SystemExit(f"the near grid reaches the ground (lowest y {grid.blocks[0].coords[..., 1].min():.4f} m)")
grid.connect()
grid.writeToCGNS(str(out/'ground.cgns'))
checks = []
for b in grid.blocks:
    p = b.coords
    a_, d_, e_ = p[1:, :-1, :-1]-p[:-1, :-1, :-1], p[:-1, 1:, :-1]-p[:-1, :-1, :-1], p[:-1, :-1, 1:]-p[:-1, :-1, :-1]
    jac = np.einsum('...i,...i->...', np.cross(a_, d_), e_)
    q = [p[:-1, :-1, :-1], p[1:, :-1, :-1], p[1:, 1:, :-1], p[:-1, 1:, :-1], p[:-1, :-1, 1:], p[1:, :-1, 1:], p[1:, 1:, 1:], p[:-1, 1:, 1:]]
    vol = sum(np.einsum('...i,...i', np.cross(q[b_]-q[a0], q[c_]-q[a0]), q[d0]-q[a0])/6
              for a0, b_, c_, d0 in [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)])
    sgn = np.sign(np.median(vol))
    # One-corner Jacobians are stricter than the cell volumes ADflow checks; a sharp,
    # upswept trailing edge can give a valid cell with one negative corner.
    checks.append({'block': b.name, 'vertices': list(p.shape[:3]), 'cells': int(jac.size),
                   'min_corner_jacobian_m3': float(jac.min()), 'nonpositive_corners': int((jac*np.sign(np.median(jac)) <= 0).sum()),
                   'nonpositive_cell_volumes': int((vol*sgn <= 0).sum()),
                   'min_cell_volume_over_median': float(np.abs(vol).min()/np.median(np.abs(vol)))})
info = {'design': os.path.abspath(own.design), 'previous_baseline': os.path.abspath(own.baseline) if own.baseline else None,
        'pivot_m': [float(pivot[0]), float(pivot[1])], 'lattice_reembedding_error_m': embed_err,
        'pyhyp': {k: v for k, v in hyp_options.items() if k != 'BC'}, 'mesh': mesh, 'near_grid_shortest_line_scale': squashed,
        'background': 'free air: the original case\'s lines above y = 0, 2.5 mm cells from -0.15 m, growing 8 % a cell to '
        'a far field at y = -10 chords (no ground)' if mesh['free_air'] else
        'identical to the original case' if mesh['ground_spacing_m'] >= .0025-1e-12 else
        f"the original case's, but {mesh['ground_spacing_m']*1e3:g} mm cells from the ground to {mesh['fine_height_m']*1e3:g} mm, "
        'growing 10 % a cell to 2.5 mm',
        'surface': {'x_min_m': float(xy[:, 0].min()), 'x_max_m': float(xy[:, 0].max()),
                    'y_min_m': float(xy[:, 1].min()), 'y_max_m': float(xy[:, 1].max())},
        'blocks': checks, 'total_cells_before_blanking': sum(c['cells'] for c in checks)}
(out/'baseline.json').write_text(json.dumps(info, indent=2)+'\n')
print(json.dumps(info, indent=2))
shutil.rmtree(scratch)
# run_chain.py takes a non-zero exit as a failed grid
assert all(c['nonpositive_cell_volumes'] == 0 for c in checks), 'Inverted cells in the new grid'
