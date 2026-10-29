"""Export the actual pyGeo lattice at the saved baseline and final designs.

Run in the same MDO Lab environment used by run_slip_optimization.py.
Reconstructs its global incidence rotation and span-tied shape functions, and
checks the reconstructed surface against the surface saved by the CFD run.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from pygeo import DVGeometry
from pyspline import Curve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    root = parser.parse_args().root
    out = root / 'ffd_geometry'
    out.mkdir(exist_ok=True)
    xy = np.loadtxt(root / 'source/airfoil_coordinates.csv', delimiter=',', skiprows=1)
    pivot = [.125, float(xy[136, 1] + .125*np.sin(np.deg2rad(3.))), 0.]
    geo = DVGeometry(str(root / 'source/ffd.xyz'))
    geo.addRefAxis('incidence_axis', curve=Curve(
        k=2, t=[0., 0., 1., 1.], coef=np.array([
            [pivot[0], pivot[1], -.005], [pivot[0], pivot[1], .055]])),
        axis=np.array([0., 0., 1.]))

    def incidence(value, geometry):
        geometry.rot_z['incidence_axis'].coef[:] = -(value[0] + 3.)

    geo.addGlobalDV('incidence', [-3.], incidence, scale=1.)
    indices = geo.getLocalIndex(0)
    shapes = [{int(indices[i, j, k]): np.array([0., 1., 0.]) for k in range(2)}
              for i in range(10) for j in range(2)]
    geo.addShapeFunctionDV('shape', shapes, scale=100.)
    geo.addPointSet(np.c_[xy, np.full(len(xy), .025)], 'surface_check')
    arrays, checks = {}, {}
    for label in ['baseline', 'final']:
        design = json.loads((root / label / 'design.json').read_text())
        geo.setDesignVars({key: np.array(value) for key, value in design.items()})
        surface = np.real(geo.update('surface_check'))
        saved = np.loadtxt(root / label / 'surface_geometry.csv', delimiter=',', skiprows=1)
        error = float(np.max(np.abs(surface - saved)))
        assert error < 1e-10, (label, error)
        lattice = np.real(geo.FFD.coef[indices]).copy()
        assert lattice.shape == (10, 2, 2, 3)
        assert np.max(np.abs(lattice[:, :, 0, :2] - lattice[:, :, 1, :2])) < 1e-12
        arrays[label] = lattice
        checks[label + '_surface_max_error_m'] = error
        rows = [[i, j, k, *lattice[i, j, k]]
                for i in range(10) for j in range(2) for k in range(2)]
        np.savetxt(out / (label + '_control_points.csv'), rows, delimiter=',',
                   header='i,j,k,x_m,y_m,z_m', comments='',
                   fmt=['%d', '%d', '%d', '%.16e', '%.16e', '%.16e'])
        with (out / (label + '_ffd.xyz')).open('w') as stream:
            stream.write('1\n10 2 2\n')
            for axis in range(3):
                np.savetxt(stream, lattice[..., axis].ravel(order='F'), fmt='%.16e')
    np.savez(out / 'control_points.npz', **arrays)
    checks.update(control_points_3d=40, projected_span_pairs=20,
                  source='pyGeo FFD coefficients after applying saved physical design variables')
    (out / 'verification.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
