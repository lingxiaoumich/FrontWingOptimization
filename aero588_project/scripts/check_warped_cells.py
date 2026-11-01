"""Find inverted cells in the warped volume grid saved with a flow solution.

  python scripts/check_warped_cells.py <flow_000_vol.cgns> [...]
Reads every structured zone's coordinates through the CGNS library and computes
each hexahedron's volume (sum of six tetrahedra, so it is also signed for
non-planar faces). A cell whose sign differs from its zone's majority is inverted.
Inverted cells are split by their overset blanking (the 'Iblank' cell field ADflow
writes: 1 computed, -1 interpolated, 0 hole, -2 flooded, -3 flood seed); only
non-live ones (iblank <= 0 except -1) are harmless. Also reports the smallest live cell.
"""
import ctypes as ct
import sys
from pathlib import Path

import numpy as np

lib = ct.CDLL(str(Path(sys.prefix)/'lib/libcgns.so'))
cgsize = ct.c_int64

def coords(path):
    fn = ct.c_int()
    if lib.cg_open(str(path).encode(), 0, ct.byref(fn)): raise RuntimeError(lib.cg_get_error())
    nz = ct.c_int(); lib.cg_nzones(fn, 1, ct.byref(nz))
    zones = []
    nsol = ct.c_int()
    for z in range(1, nz.value+1):
        name = ct.create_string_buffer(64); size = (cgsize*9)()
        lib.cg_zone_read(fn, 1, z, name, size)
        n = [size[0], size[1], size[2]]
        xyz = []
        for c in [b'CoordinateX', b'CoordinateY', b'CoordinateZ']:
            a = np.zeros(n[0]*n[1]*n[2]); lo = (cgsize*3)(1, 1, 1); hi = (cgsize*3)(*n)
            if lib.cg_coord_read(fn, 1, z, c, 4, lo, hi, a.ctypes.data_as(ct.c_void_p)): raise RuntimeError(c)
            xyz.append(a.reshape(n[::-1]).transpose(2, 1, 0))  # (i, j, k)
        lib.cg_nsols(fn, 1, z, ct.byref(nsol)); blank = None
        for sidx in range(1, nsol.value+1):
            nf = ct.c_int(); lib.cg_nfields(fn, 1, z, sidx, ct.byref(nf))
            for f in range(1, nf.value+1):
                dt = ct.c_int(); fname = ct.create_string_buffer(64)
                lib.cg_field_info(fn, 1, z, sidx, f, ct.byref(dt), fname)
                if fname.value.lower() in (b'blank', b'iblank'):
                    m = [x-1 for x in n]; b = np.zeros(m[0]*m[1]*m[2])
                    lo = (cgsize*3)(1, 1, 1); hi = (cgsize*3)(*m)
                    lib.cg_field_read(fn, 1, z, sidx, fname, 4, lo, hi, b.ctypes.data_as(ct.c_void_p))
                    blank = b.reshape(m[::-1]).transpose(2, 1, 0)
        zones.append((name.value.decode(), np.stack(xyz, -1), blank))
    lib.cg_close(fn)
    return zones

def hex_volumes(X):
    p = [X[:-1, :-1, :-1], X[1:, :-1, :-1], X[1:, 1:, :-1], X[:-1, 1:, :-1],
         X[:-1, :-1, 1:], X[1:, :-1, 1:], X[1:, 1:, 1:], X[:-1, 1:, 1:]]
    tets = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    return sum(np.einsum('...i,...i', np.cross(p[b]-p[a], p[c]-p[a]), p[d]-p[a])/6 for a, b, c, d in tets)

for path in sys.argv[1:]:
    print(f'== {path}')
    for name, X, blank in coords(path):
        v = hex_volumes(X)
        s = np.sign(np.median(v)); bad = v*s <= 0
        line = f'  {name}: {X.shape[:3]} nodes, {v.size} cells, inverted {bad.sum()}'
        if blank is not None and bad.any():
            kinds = {int(k): int(((blank == k) & bad).sum()) for k in np.unique(blank[bad])}
            label = {1: 'computed', 0: 'blanked', -1: 'interpolated'}
            line += ' (by blanking: '+', '.join(f'{label.get(k, k)} {c}' for k, c in kinds.items())+')'
        if blank is not None:
            live = blank != 0
            live &= blank > -2  # computed (1) and interpolated (-1) cells
            r = np.abs(v[live])/np.median(np.abs(v))
            line += f'; live cells: min |V|/median {r.min():.1e}, ratio < 0.01: {(r < .01).sum()}'
        print(line)
        if bad.any():
            c = (X[:-1, :-1, :-1]+X[1:, 1:, 1:])[bad]/2
            print(f'     inverted cells span x {c[:, 0].min():.3f}–{c[:, 0].max():.3f} m, y {c[:, 1].min():.3f}–{c[:, 1].max():.3f} m')
