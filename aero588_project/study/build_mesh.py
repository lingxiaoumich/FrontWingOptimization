"""Build one level of the mesh family for the repo design case (overset, repo BCs).

    mpirun -np 4 python study/build_mesh.py --level L2 --alpha 0

Writes study/grids/<case>/: <case>.cgns, near.cgns, surface.xyz, ffd.xyz, case.json.
"""

import argparse
import json
import time

import numpy as np
from mpi4py import MPI

from common import CASE, LEVELS, MARCH, STUDY, background_axes, case_name, section, write_ffd, write_surface_xyz

comm = MPI.COMM_WORLD
p = argparse.ArgumentParser()
p.add_argument("--level", default="L2", choices=list(LEVELS))
p.add_argument("--alpha", type=float, default=0.0, help="incidence, deg; negative = nose down (repo)")
p.add_argument("--nzNear", type=int, default=1, help="spanwise cells of the O-mesh (1 = repo)")
args = p.parse_args()

from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import Block, Boco, Grid, combineGrids, readGrid

started = time.perf_counter()
lv = LEVELS[args.level]
chord, span = CASE["chord"], CASE["span"]
name = case_name(args.level, args.alpha) + ("" if args.nzNear == 1 else f"_nz{args.nzNear}")
out = STUDY / "grids" / name
if comm.rank == 0:
    out.mkdir(parents=True, exist_ok=True)
comm.Barrier()

_, (xw, yw), pivot = section(args.level, args.alpha)
if comm.rank == 0:
    write_surface_xyz(out / "surface.xyz", xw, yw, span, nz=args.nzNear + 1)
    write_ffd(out / "ffd.xyz", args.alpha, chord, span, pivot[1])
comm.Barrier()

hyp = pyHyp(options={
    "inputFile": str(out / "surface.xyz"), "unattachedEdgesAreSymmetry": False, "outerFaceBC": "overset",
    "autoConnect": True, "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}}, "families": "airfoil",
    "N": lv["N"], "s0": lv["s0"], "marchDist": MARCH,
    "ps0": -1.0, "pGridRatio": -1.0, "cMax": 3.0, "epsE": 1.0, "epsI": 2.0,
    "theta": 3.0, "volCoef": 0.25, "volBlend": 0.0001, "volSmoothIter": 100,
})
hyp.run()
hyp.writeCGNS(str(out / "near.cgns"))

if comm.rank == 0:
    xx, yy = background_axes(args.level)
    coords = np.stack(np.meshgrid(xx, yy, [0.0, span], indexing="ij"), axis=-1)
    dims = list(coords.shape[:3])
    block = Block("background", dims, coords)

    def boundary(bname, kind, axis, index, family):
        pr = np.array([[1, dims[0]], [1, dims[1]], [1, dims[2]]])
        pr[axis] = index
        block.addBoco(Boco(bname, kind, pr, family))

    boundary("inlet", "bcfarfield", 0, 1, "farfield")
    boundary("outlet", "bcfarfield", 0, dims[0], "farfield")
    boundary("top", "bcfarfield", 1, dims[1], "farfield")
    boundary("ground", "bcsymmetryplane", 1, 1, "ground")
    boundary("span0", "bcsymmetryplane", 2, 1, "symmetry")
    boundary("span1", "bcsymmetryplane", 2, 2, "symmetry")
    background = Grid()
    background.blocks = [block]
    grid = combineGrids([readGrid(str(out / "near.cgns")), background])
    grid.connect()
    gridFile = out / f"{name}.cgns"
    grid.writeToCGNS(str(gridFile))

    checks = []
    for b in grid.blocks:
        P = b.coords
        a = P[1:, :-1, :-1] - P[:-1, :-1, :-1]
        d = P[:-1, 1:, :-1] - P[:-1, :-1, :-1]
        e = P[:-1, :-1, 1:] - P[:-1, :-1, :-1]
        jac = np.einsum("...i,...i->...", np.cross(a, d), e)
        checks.append({"block": b.name, "vertices": list(P.shape[:3]), "cells": int(jac.size),
                       "nonpositive_corners": int((jac <= 0).sum()), "lowest_grid_y_m": float(P[..., 1].min())})
    assert all(c["nonpositive_corners"] == 0 for c in checks), "Invalid corner Jacobian"
    assert yw.min() >= CASE["clearance"] - 1e-9, "Insufficient surface clearance"
    near = [c for c in checks if c["vertices"] == [lv["nSurf"] + lv["nTE"] + 1, args.nzNear + 1, lv["N"]]]
    assert len(near) == 1 and near[0]["lowest_grid_y_m"] > 0, "Near grid crosses ground"
    i_le, i_te = int(np.argmin(xw)), 0
    le, te = (xw[i_le], yw[i_le]), (xw[i_te], 0.5 * (yw[0] + yw[lv["nSurf"]]))
    meta = dict(
        name=name, level=args.level, alpha=args.alpha, nzNear=args.nzNear, **CASE, march=MARCH, **lv,
        gridFile=str(gridFile), ffdFile=str(out / "ffd.xyz"), pivot=list(pivot),
        lowest_surface_y=float(yw.min()),
        LE_to_TE_angle_deg=float(np.degrees(np.arctan2(te[1] - le[1], te[0] - le[0]))),
        nCells=sum(c["cells"] for c in checks), blocks=checks, build_s=time.perf_counter() - started,
    )
    (out / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("name", "nCells", "lowest_surface_y", "LE_to_TE_angle_deg")}))
