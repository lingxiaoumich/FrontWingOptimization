"""KKT residual of a stored PSQP point, independent of PSQP's own test (2026-11-27, runs 036-039).

From the optimizer database (optimization.hst): the gradient of the objective and the Jacobians of the constraints
within `tol` of a bound (the active set), in the variables as PSQP sees them (shape variables x dv_scale 100).
The multipliers come from a sign-constrained least-squares fit (scipy nnls): grad f = sum lam_i grad c_i, lam >= 0 for
an upper bound, <= 0 for a lower one. The residual is the largest component of grad f - J_A^T lam. PSQP's own test
(exit 4) uses its QP multipliers and the same max-norm with TOLG 1e-4, so the two agree to about that level but not
exactly: nnls minimizes the 2-norm of the residual, not the max-norm.
  python kkt_residual.py ROUND [n_last] [act_tol]          (MACH-Aero env: needs pyoptsparse)
"""
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

DV_SCALE = 100.
LUMPED = ('lift', 'cp_min', 'max_thickness', 'section_area')

def bounds_of(root):
    """{constraint group: [(lower, upper)]} from the round's problem.txt."""
    sec = (Path(root)/'problem.txt').read_text().split('Constraints')[1]
    out = {}
    for m in re.finditer(r'^\s*\d+\s+(\S+)\s+i\s+(\S+)\s+\S+\s+(\S+)', sec, re.M):
        name = m.group(1) if m.group(1) in LUMPED else re.sub(r'_\d+$', '', m.group(1))
        out.setdefault(name, []).append((float(m.group(2)), float(m.group(3))))
    return out

def residual(funcs, sens, bounds, tol=1e-3):
    """(max |residual|, [(constraint, multiplier)] sorted by size, number of active constraints) at one point."""
    g = np.asarray(sens['objective']['shape']).ravel()/DV_SCALE
    rows, names = [], []
    for nm, bl in bounds.items():
        if nm not in funcs: continue
        v = np.atleast_1d(funcs[nm]); J = np.atleast_2d(np.asarray(sens[nm]['shape']))/DV_SCALE
        if J.shape[0] != len(v): J = J.reshape(len(v), -1)
        for i, (lo, up) in enumerate(bl[:len(v)]):
            sc = max(1., abs(lo) if lo > -1e19 else 0, abs(up) if up < 1e19 else 0)
            if up < 1e19 and v[i] >= up-tol*sc: rows.append(-J[i]); names.append(f'{nm}[{i}] upper')
            if lo > -1e19 and v[i] <= lo+tol*sc: rows.append(J[i]); names.append(f'{nm}[{i}] lower')
    if rows:
        A = np.array(rows).T; lam, _ = nnls(A, g); r = g-A@lam
    else: lam, r = np.array([]), g
    mult = [(names[i], float(lam[i])) for i in np.argsort(-lam) if lam[i] > 1e-8]
    return float(np.abs(r).max()), mult, len(names)

def last_points(root, n=1, tol=1e-3):
    """[(call, objective, max residual, multipliers, n active)] for the last n gradient points of a round."""
    from pyoptsparse import History
    root = Path(root); b = bounds_of(root)
    h = History(str(root/'optimization.hst'), flag='r')
    keys = [k for k in sorted(h.getCallCounters(), key=int) if 'funcsSens' in h.read(k)]
    out = []
    for k in keys[-n:]:
        c = h.read(k); res, mult, na = residual(c['funcs'], c['funcsSens'], b, tol)
        out.append((int(k), float(np.atleast_1d(c['funcs']['objective'])[0]), res, mult, na))
    return out

if __name__ == '__main__':
    root = Path(sys.argv[1]); n = int(sys.argv[2]) if len(sys.argv) > 2 else 3; tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    for k, obj, res, mult, na in last_points(root, n, tol):
        print(f'call {k:4d} objective {obj:.7f} KKT residual {res:.2e} ({na} active, {len(mult)} with a multiplier)')
    print('multipliers at the last point:', [(a, round(b, 4)) for a, b in mult[:12]])
