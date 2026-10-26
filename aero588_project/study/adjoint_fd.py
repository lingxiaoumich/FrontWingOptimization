"""Adjoint versus central finite differences for the design variables of the optimization
(docs/shape_control_and_optimization.md): incidence, ride height, FFD shape variables.

    mpirun -np 8 python study/adjoint_fd.py --case study/grids/<case> [--solver ank_uncoupled]
        [--shape 1:0 4:0 4:1 8:0 8:1] [--steps 1e-3 1e-4 1e-5 1e-6]

Shape variables are named station:row (row 0 = lower, 1 = upper control point). FD steps are
in each variable's own unit (deg for incidence, m for ride height and shape); the incidence
steps are the listed steps x 100 (1e-1 ... 1e-4 deg). Every FD solve restarts from the current
state and converges to the same --l2 relative to the freestream residual.

Writes study/results/adjoint/<case>__<solver>__<overset>.json.
"""

import argparse
import json
import time

import numpy as np
from mpi4py import MPI

from common import CASE, SOLVERS, STUDY, adflow_options, add_functions, aero_problem, load_meta

comm = MPI.COMM_WORLD
p = argparse.ArgumentParser()
p.add_argument("--case", required=True)
p.add_argument("--solver", default="repo", choices=list(SOLVERS))
p.add_argument("--overset", default="frozen", choices=["frozen", "full", "fast"])
p.add_argument("--l2", type=float, default=1e-12)
p.add_argument("--shape", nargs="*", default=["1:0", "4:0", "4:1", "8:0", "8:1"])
p.add_argument("--steps", type=float, nargs="+", default=[1e-3, 1e-4, 1e-5, 1e-6])
p.add_argument("--globals", nargs="*", default=["incidence", "rideheight"])
p.add_argument("--checkOnly", action="store_true", help="stop after the parameterization checks")
p.add_argument("--alpha0", type=float, default=0.0,
               help="ADflow angle of attack, deg. 0 puts the freestream parallel to the top far-field face, whose "
                    "inflow/outflow switch then sits on the baseline and adds a jump to every FD; 0.05 moves it off")
p.add_argument("--exactWallDist", action="store_true", help="useApproxWallDistance = False")
p.add_argument("--tag", default="")
args = p.parse_args()

from adflow import ADFLOW
from idwarp import USMesh
from pygeo import DVGeometry
from pyspline import Curve

meta = load_meta(args.case)
span = CASE["span"]
outdir = STUDY / "results" / "adjoint"
if comm.rank == 0:
    outdir.mkdir(parents=True, exist_ok=True)
comm.Barrier()

# ---------------------------------------------------------------- parameterization
DVGeo = DVGeometry(meta["ffdFile"])
xc, yc = meta["pivot"]
DVGeo.addRefAxis("c4", curve=Curve(x=np.array([xc, xc]), y=np.array([yc, yc]),
                                   z=np.array([-0.25 * span, 1.25 * span]), k=2), axis="x")


def incidence(val, geo):
    """Change of incidence, deg. Repo sign: negative = nose down, i.e. rotate the box by -val."""
    for i in range(geo.rot_z["c4"].coef.shape[0]):
        geo.rot_z["c4"].coef[i] = -val[0]


def rideheight(val, geo):
    """Vertical translation of the box, m."""
    C = geo.extractCoef("c4")
    C[:, 1] += val[0]
    geo.restoreCoef(C, "c4")


DVGeo.addGlobalDV("incidence", value=0.0, func=incidence, lower=-4.0 - meta["alpha"], upper=4.0 - meta["alpha"])
DVGeo.addGlobalDV("rideheight", value=0.0, func=rideheight, lower=-0.05, upper=0.10)
DVGeo.addSpanwiseLocalDV("shape", "k", axis="y", lower=-0.05 * CASE["chord"], upper=0.05 * CASE["chord"])
nShape = DVGeo.getNDV() - 2

opts = adflow_options(meta, outdir, args.solver, l2=args.l2, overset=args.overset)
if args.exactWallDist:
    opts["useApproxWallDistance"] = False
solver = ADFLOW(options=opts, comm=comm)
ap = aero_problem(meta, name="fwo")
ap.alpha = args.alpha0
add_functions(solver, ap, funcs=("cl", "cd"))
FUNCS = ["cl_airfoil", "cd_airfoil"]
solver.setDVGeo(DVGeo)
solver.setMesh(USMesh(options={"gridFile": meta["gridFile"]}, comm=comm))

# ---------------------------------------------------------------- parameterization self-checks
pts = solver.getSurfaceCoordinates(solver.designFamilyGroup)
DVGeo.addPointSet(pts, "check")


def gmax(a):
    return comm.allreduce(float(np.max(a)) if a.size else -np.inf, op=MPI.MAX)


def gmin(a):
    return comm.allreduce(float(np.min(a)) if a.size else np.inf, op=MPI.MIN)


x0 = DVGeo.getValues()
embed_err = gmax(np.abs(DVGeo.update("check") - pts))

# shape DV index -> (station i, row j) of the control point it moves
coef0 = DVGeo.FFD.coef.copy()
lidx = DVGeo.getLocalIndex(0)
dvmap = {}
for d in range(nShape):
    x = {k: np.array(v, dtype=float).copy() for k, v in x0.items()}
    x["shape"][d] = 1e-3
    DVGeo.setDesignVars(x)
    DVGeo.update("check")
    moved = np.where(np.abs(DVGeo.FFD.coef[:, 1] - coef0[:, 1]) > 1e-9)[0]
    ijk = [tuple(int(v) for v in np.argwhere(lidx == m)[0]) for m in moved]
    dvmap[f"{ijk[0][0]}:{ijk[0][1]}"] = d
DVGeo.setDesignVars(x0)

# sign check: -1 deg of incidence must lower the LE (nose down) relative to the TE
x = {k: np.array(v, dtype=float).copy() for k, v in x0.items()}
x["incidence"][0] = -1.0
DVGeo.setDesignVars(x)
moved_pts = DVGeo.update("check")
DVGeo.setDesignVars(x0)
DVGeo.update("check")
xle, xte = gmin(pts[:, 0]), gmax(pts[:, 0])
near_le = np.abs(pts[:, 0] - xle) < 1e-3
near_te = np.abs(pts[:, 0] - xte) < 1e-3
dy_le = comm.allreduce(float(np.sum(moved_pts[near_le, 1] - pts[near_le, 1])), op=MPI.SUM) / max(
    comm.allreduce(int(near_le.sum()), op=MPI.SUM), 1)
dy_te = comm.allreduce(float(np.sum(moved_pts[near_te, 1] - pts[near_te, 1])), op=MPI.SUM) / max(
    comm.allreduce(int(near_te.sum()), op=MPI.SUM), 1)
checks = dict(embedding_error_m=embed_err, nShape=nShape, dvmap=dvmap,
              incidence_minus1deg_dyLE_m=dy_le, incidence_minus1deg_dyTE_m=dy_te,
              sign_ok=bool(dy_le < 0 < dy_te))
if comm.rank == 0:
    print("parameterization checks:", json.dumps(checks), flush=True)
if args.checkOnly:
    raise SystemExit(0)


# ---------------------------------------------------------------- flow, adjoint, FD
def evaluate(x):
    DVGeo.setDesignVars(x)
    t = time.perf_counter()
    solver(ap)
    f = {}
    solver.evalFunctions(ap, f, evalFuncs=FUNCS)
    solver.checkSolutionFailure(ap, f)
    out = {k: float(f[f"fwo_{k}"]) for k in FUNCS}
    out["ratio"] = out["cl_airfoil"] / out["cd_airfoil"]
    out["fail"] = bool(f.get("fail", False))
    out["time_s"] = time.perf_counter() - t
    return out


t = time.perf_counter()
f0 = evaluate(x0)
t_flow = time.perf_counter() - t
t = time.perf_counter()
sens = {}
solver.evalFunctionsSens(ap, sens, evalFuncs=FUNCS)
solver.checkAdjointFailure(ap, sens)
t_adj = time.perf_counter() - t


def adjoint(dv, i):
    g = {k: float(np.ravel(sens[f"fwo_{k}"][dv])[i]) for k in FUNCS}
    g["ratio"] = (g["cl_airfoil"] * f0["cd_airfoil"] - f0["cl_airfoil"] * g["cd_airfoil"]) / f0["cd_airfoil"] ** 2
    return g


variables = [(g, 0) for g in args.globals] + [("shape", dvmap[s]) for s in args.shape]
labels = {("shape", dvmap[s]): f"shape[{s}]" for s in args.shape}
rows = []
for dv, i in variables:
    label = labels.get((dv, i), dv)
    adj = adjoint(dv, i)
    for h0 in args.steps:
        h = h0 * 100 if dv == "incidence" else h0
        fp, fm = None, None
        for sgn in (+1, -1):
            x = {k: np.array(v, dtype=float).copy() for k, v in x0.items()}
            x[dv][i] += sgn * h
            r = evaluate(x)
            fp, fm = (r, fm) if sgn > 0 else (fp, r)
        for fn in ("cl_airfoil", "cd_airfoil", "ratio"):
            fd = (fp[fn] - fm[fn]) / (2 * h)
            rows.append(dict(variable=label, function=fn, step=h, adjoint=adj[fn], fd=fd,
                             rel_diff=abs(adj[fn] - fd) / max(abs(fd), 1e-30),
                             fail=fp["fail"] or fm["fail"], solve_s=0.5 * (fp["time_s"] + fm["time_s"])))
        if comm.rank == 0:
            r = rows[-3]
            print(f"{label:14s} h={h:8.1e}  dCL adj {r['adjoint']: .6e}  fd {r['fd']: .6e}  rel {r['rel_diff']:.2e}", flush=True)
    evaluate(x0)  # back to the baseline state before the next variable

if comm.rank == 0:
    res = dict(case=meta["name"], level=meta["level"], alpha=meta["alpha"], nCells=meta["nCells"], solver=args.solver,
               overset=args.overset, l2=args.l2, nProcs=comm.size, alpha0=args.alpha0, exactWallDist=args.exactWallDist,
               checks=checks, baseline=f0, flow_s=t_flow, adjoint_s=t_adj, nAdjoints=len(FUNCS), rows=rows)
    tag = f"__{args.tag}" if args.tag else ""
    (outdir / f"{meta['name']}__{args.solver}__{args.overset}{tag}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(f"flow {t_flow:.1f} s, {len(FUNCS)} adjoints {t_adj:.1f} s on {comm.size} ranks")
