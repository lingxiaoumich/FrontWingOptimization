"""Finite-difference step study of the adjoint gradient at one design.
Directions: each function's own adjoint gradient in the shape variables, and incidence.
Every perturbed solve restarts from the base design's converged state.
  mpirun -np 64 python fd_study.py --script S/run_slip_optimization.py --at D.json --out R.json <run args>"""
import argparse, copy, json, sys
from pathlib import Path
import numpy as np
h_ = argparse.ArgumentParser(add_help=False)
h_.add_argument('--script', type=Path, required=True)
h_.add_argument('--at', type=Path, required=True)
h_.add_argument('--out', type=Path, required=True)
h_.add_argument('--steps', default='1e-5,3e-5,1e-4,3e-4', help='m for shape directions; x100 deg for incidence')
own, rest = h_.parse_known_args()
sys.argv = [str(own.script)]+rest
sys.path.insert(0, str(own.script.parent))
exec(compile(own.script.read_text().split('# Baseline always starts from the same NACA shape')[0], str(own.script), 'exec'))
names = ['objective', 'drag']
x0 = load_design(own.at)
f0, fail = objfun(x0); assert not fail
base = solver.getStates().copy()
g, adjfail = sensfun(x0, f0); assert not adjfail
def solve_from_base(x):
    global anchor_state, anchor_x
    solver.setStates(base); anchor_state, anchor_x = base.copy(), copy.deepcopy(x0)
    f, fl = objfun(x)
    return {n: float(f[n]) for n in names}, bool(fl)
dirs = {}
for n in names:
    u = {k: np.asarray(g[n][k], float).ravel() for k in shape_groups}
    nrm = np.sqrt(sum((v**2).sum() for v in u.values()))
    dirs[n+'_shape'] = ({k: v/nrm for k, v in u.items()}, 'm', 1.)
dirs['incidence'] = ({'incidence': np.array([1.])}, 'deg', 100.)
out = {'design': str(own.at), 'f0': {n: float(f0[n]) for n in names}, 'directions': {}}
for dname, (u, unit, mult) in dirs.items():
    adj = {n: float(sum(np.asarray(g[n][k], float).ravel()@u[k] for k in u)) for n in names}
    rows = []
    for s in own.steps.split(','):
        step = float(s)*mult; r = {'step': step}
        for sgn, tag in [(1, 'plus'), (-1, 'minus')]:
            x = copy.deepcopy(x0)
            for k in u: x[k] = x0[k]+sgn*step*u[k]
            r[tag], r[tag+'_failed'] = solve_from_base(x)
        r['fd'] = {n: (r['plus'][n]-r['minus'][n])/(2*step) for n in names}
        r['second_difference_over_first'] = {n: abs(r['plus'][n]+r['minus'][n]-2*out['f0'][n])/max(abs(r['plus'][n]-r['minus'][n]), 1e-30) for n in names}
        rows.append(r)
        if comm.rank == 0: print('FD_STUDY', dname, json.dumps(r), flush=True)
    out['directions'][dname] = {'unit': unit, 'adjoint': adj, 'rows': rows}
# Path dependence: re-solve the base design from the last perturbed state.
f, fl = objfun(x0)
out['base_resolved_from_other_state'] = {n: float(f[n]) for n in names}
dump(own.out, out)
if comm.rank == 0: print('FD_STUDY_DONE', flush=True)
