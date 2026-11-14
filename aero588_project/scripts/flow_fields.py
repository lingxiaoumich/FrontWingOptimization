"""Read the flow fields of an ADflow volume solution and triangulate them for plotting."""
import ctypes as ct
import sys
from pathlib import Path

import numpy as np
import matplotlib.tri as mtri

lib = ct.CDLL(str(Path(sys.prefix)/'lib/libcgns.so')); cgsize = ct.c_int64
def read_zones(path):
    """Vertex x, y and cell-centred fields of each zone's first z plane."""
    fn = ct.c_int()
    if lib.cg_open(str(path).encode(), 0, ct.byref(fn)): raise RuntimeError(path)
    zones = {}
    for z in (1, 2):
        name = ct.create_string_buffer(64); size = (cgsize*9)(); lib.cg_zone_read(fn, 1, z, name, size)
        n, m = list(size[:3]), list(size[3:6])
        near = b'airfoil' in name.value
        plane = (lambda v: v[:, 0, :]) if near else (lambda v: v[:, :, 0])  # near grid is (i, z, j)
        def read(field, dims, coord):
            buf = np.empty(int(np.prod(dims))); lo = (cgsize*3)(1, 1, 1); hi = (cgsize*3)(*dims)
            args = (fn, 1, z, field, 4, lo, hi, buf.ctypes.data_as(ct.c_void_p))
            code = lib.cg_coord_read(*args) if coord else lib.cg_field_read(fn, 1, z, 1, *args[3:])
            if code: raise RuntimeError(field)
            return plane(buf.reshape(dims, order='F'))
        f = {k: read(k.encode(), n, True) for k in ['CoordinateX', 'CoordinateY']}
        f.update({k: read(k.encode(), m, False) for k in ['CoefPressure', 'Mach', 'VelocityX', 'VelocityY', 'Iblank']})
        zones['near' if near else 'background'] = f
    lib.cg_close(fn)
    return zones

def add_cpt(zones, mach, g=1.4):
    """Cpt = 1 - (p0,inf - p0)/q_inf: 1 in the freestream, about Cp in dead air."""
    p0_inf = (1+.2*mach**2)**3.5  # over p_inf; q_inf/p_inf = g M^2/2
    for f in zones.values():
        pp = 1+f['CoefPressure']*g*mach**2/2
        f['Cpt'] = 1-(p0_inf-pp*(1+.2*f['Mach']**2)**3.5)/(g*mach**2/2)

def centers(x): return .25*(x[:-1, :-1]+x[1:, :-1]+x[1:, 1:]+x[:-1, 1:])
def tri(f, near):
    x, y = centers(f['CoordinateX']), centers(f['CoordinateY'])
    ids = np.arange(x.size).reshape(x.shape)
    quads = [[ids[:-1, :-1], ids[1:, :-1], ids[1:, 1:], ids[:-1, 1:]]]
    if near: quads.append([ids[-1, :-1], ids[0, :-1], ids[0, 1:], ids[-1, 1:]])  # O-grid seam
    t = np.concatenate([np.r_[np.c_[q[0].ravel(), q[1].ravel(), q[2].ravel()], np.c_[q[0].ravel(), q[2].ravel(), q[3].ravel()]] for q in quads])
    live = np.isin(f['Iblank'], [1, -1]).ravel()
    return mtri.Triangulation(x.ravel(), y.ravel(), t[live[t].all(axis=1)])
