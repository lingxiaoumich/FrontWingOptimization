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
own = h.parse_args()
out = own.out.resolve(); out.mkdir(parents=True, exist_ok=True)

# The optimizer's own geometry setup, in a scratch root (it writes setup files there).
scratch = Path(tempfile.mkdtemp())
shutil.copytree(P/'setups/slip_ground_xy5_cd005/source', scratch/'source')
if own.baseline:
    for f in ['airfoil_coordinates.csv']: shutil.copy(own.baseline.parent/f, scratch/'source'/f)
sys.path.insert(0, str(scratch/'source'))
sys.argv = ['x', '--root', str(scratch), '--ffd-stations', '5', '--shape-directions', 'xy'] + \
           (['--baseline', str(own.baseline.resolve())] if own.baseline else [])
S = {}
exec(compile((P/'scripts/run_slip_optimization.py').read_text().split('options={')[0], 'setup', 'exec'), S)
geo, pivot, span, c0 = S['geo'], S['pivot'], S['span'], S['c0']
design = S['load_design'](own.design)
geo.setDesignVars(design)
pts = geo.update('surface_check')

# Deformed lattice: rotate the undeformed one about the pivot, then add the shape moves
# (pyGeo's order; checked below by re-embedding the surface in the new lattice).
lines = (scratch/'source/ffd.xyz').read_text().split('\n')
dims = [int(v) for v in lines[1].split()]
vals = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
L0 = np.stack([v.reshape(dims, order='F') for v in vals], -1)  # (station, lower/upper, span, xyz)
t = np.deg2rad(-(design['incidence'][0]+3.))
R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
L = L0.copy()
L[..., :2] = (L0[..., :2]-pivot[:2])@R.T+pivot[:2]
for name, ax in [('shape_x', 0), ('shape_y', 1)]:
    L[..., ax] += np.reshape(design[name], dims[:2])[:, :, None]
with (out/'ffd.xyz').open('w') as f:
    f.write(f'1\n{dims[0]} {dims[1]} {dims[2]}\n')
    for axis in range(3): np.savetxt(f, L[..., axis].ravel(order='F'), fmt='%.16e')
from pygeo import DVGeometry
check = DVGeometry(str(out/'ffd.xyz')); check.addPointSet(pts, 'p')
embed_err = float(np.abs(check.update('p')-pts).max())
assert embed_err < 1e-9, f'design surface does not re-embed in the deformed lattice ({embed_err:.1e} m)'

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
    'N': 81, 's0': 2e-6, 'marchDist': .025,
    'ps0': -1., 'pGridRatio': -1., 'cMax': 3., 'epsE': 1., 'epsI': 2.,
    'theta': 3., 'volCoef': .25, 'volBlend': .0001, 'volSmoothIter': 100}
hyp_options.update(json.loads(own.hyp_options))
hyp = pyHyp(options={'inputFile': str(out/'airfoil_surface.xyz'), **hyp_options})
hyp.run()
hyp.writeCGNS(str(out/'airfoil_near.cgns'))

def extend(start, stop, step, ratio=1.08):
    values = [start]; direction = np.sign(stop-start)
    while abs(stop-values[-1]) > step:
        values.append(values[-1]+direction*step); step *= ratio
    values.append(stop)
    return np.array(values)
chord = c0  # background spacing of the original case, which this grid must match
xx = np.r_[extend(-.15, -10*chord, .005)[::-1][:-1], np.linspace(-.15, .85, 201), extend(.85, 15*chord, .005)[1:]]
yy = np.r_[np.linspace(0, .35, 141), extend(.35, 10*chord, .0025)[1:]]
coordinates = np.stack(np.meshgrid(xx, yy, [0, span], indexing='ij'), axis=-1)
bdims = list(coordinates.shape[:3])
block = Block('background', bdims, coordinates)
def boundary(name, kind, axis, index, family):
    point_range = np.array([[1, bdims[0]], [1, bdims[1]], [1, bdims[2]]]); point_range[axis] = index
    block.addBoco(Boco(name, kind, point_range, family))
boundary('inlet', 'bcfarfield', 0, 1, 'farfield'); boundary('outlet', 'bcfarfield', 0, bdims[0], 'farfield')
boundary('top', 'bcfarfield', 1, bdims[1], 'farfield'); boundary('ground', 'bcsymmetryplane', 1, 1, 'ground')
boundary('span0', 'bcsymmetryplane', 2, 1, 'symmetry'); boundary('span1', 'bcsymmetryplane', 2, 2, 'symmetry')
background = Grid(); background.blocks = [block]
grid = combineGrids([readGrid(str(out/'airfoil_near.cgns')), background])
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
                   'min_cell_volume_over_median': float(np.abs(vol).min()/np.median(np.abs(vol))),
                   'lowest_grid_y_m': float(p[..., 1].min())})
info = {'design': str(own.design.resolve()), 'previous_baseline': str(own.baseline.resolve()) if own.baseline else None,
        'pivot_m': [float(pivot[0]), float(pivot[1])], 'lattice_reembedding_error_m': embed_err,
        'pyhyp': {k: v for k, v in hyp_options.items() if k != 'BC'}, 'background': 'identical to the original case',
        'surface': {'x_min_m': float(xy[:, 0].min()), 'x_max_m': float(xy[:, 0].max()),
                    'y_min_m': float(xy[:, 1].min()), 'y_max_m': float(xy[:, 1].max())},
        'blocks': checks, 'total_cells_before_blanking': sum(c['cells'] for c in checks)}
(out/'baseline.json').write_text(json.dumps(info, indent=2)+'\n')
print(json.dumps(info, indent=2))
assert all(c['nonpositive_cell_volumes'] == 0 for c in checks), 'Inverted cells in the new grid'
assert float(grid.blocks[0].coords[..., 1].min()) > 0, 'Near grid crosses ground'
shutil.rmtree(scratch)
