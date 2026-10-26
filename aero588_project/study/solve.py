"""Solve one case with one solver variant and keep the full convergence history.

    mpirun -np 8 python study/solve.py --case study/grids/<case> --solver repo [--l2 1e-12]

Writes study/results/solve/<case>__<solver>.json: forces, wall time, iterations, residual
reduction, failure flags and the per-iteration history (residual, C_L, C_D, CPU time).
"""

import argparse
import json
import time

import numpy as np
from mpi4py import MPI

from common import FUNCS, SOLVERS, STUDY, adflow_options, add_functions, aero_problem, flow_conditions, load_meta

comm = MPI.COMM_WORLD
p = argparse.ArgumentParser()
p.add_argument("--case", required=True)
p.add_argument("--solver", default="repo", choices=list(SOLVERS))
p.add_argument("--l2", type=float, default=1e-12, help="target residual relative to freestream")
p.add_argument("--ncycles", type=int, default=20000)
p.add_argument("--tag", default="")
args = p.parse_args()

from adflow import ADFLOW

meta = load_meta(args.case)
outdir = STUDY / "results" / "solve"
if comm.rank == 0:
    outdir.mkdir(parents=True, exist_ok=True)
comm.Barrier()

t0 = time.perf_counter()
solver = ADFLOW(options=adflow_options(meta, outdir, args.solver, l2=args.l2, ncycles=args.ncycles), comm=comm)
ap = aero_problem(meta, name="fwo")
add_functions(solver, ap)
t_setup = time.perf_counter() - t0
t0 = time.perf_counter()
solver(ap)
t_solve = time.perf_counter() - t0
funcs = {}
solver.evalFunctions(ap, funcs)
initial, start, final = solver.getResNorms()
try:
    hist = solver.getConvergenceHistory()
except UnicodeDecodeError:  # ADflow decodes the iteration-type strings as utf-8; they can hold junk bytes
    hist = {}

if comm.rank == 0:
    fc = flow_conditions()
    q = 0.5 * fc["rho"] * fc["V"] ** 2 * meta["chord"]
    res = {f: float(funcs[f"fwo_{f}_airfoil"]) for f in FUNCS}
    h = {k: (np.asarray(v).tolist() if not isinstance(v, list) else [str(x) for x in v]) for k, v in hist.items()}
    row = dict(
        case=meta["name"], level=meta["level"], alpha=meta["alpha"], nCells=meta["nCells"], solver=args.solver,
        tag=args.tag, nProcs=comm.size, l2_target=args.l2, setup_s=t_setup, solve_s=t_solve,
        residual_reduction=float(final / initial), solve_failed=bool(ap.solveFailed), fatal_fail=bool(ap.fatalFail),
        **res, downforce_to_drag=-res["cl"] / res["cd"],
        downforce_N_per_m=-res["cl"] * q, drag_N_per_m=res["cd"] * q, mach=fc["mach"], reynolds=fc["reynolds"],
        history=h,
    )
    tag = f"_{args.tag}" if args.tag else ""
    (outdir / f"{meta['name']}__{args.solver}{tag}.json").write_text(json.dumps(row) + "\n")
    print(json.dumps({k: v for k, v in row.items() if k != "history"}, indent=2))
