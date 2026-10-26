"""Shared set-up for the solver, mesh and derivative studies of the repo's design case.

Everything that defines the case lives here: flow conditions, the mesh family, the geometry
recipe, the boundary conditions and the ADflow solver variants. The geometry, boundary
conditions and sign convention are those of scripts/run_ground_case.py:

  * NACA 0012, standard finite trailing edge (-0.1015) closed by a short vertical face;
  * negative incidence = nose down = LE below TE, rotation about c/4, freestream horizontal;
  * the section is lifted so its lowest point is exactly `clearance` above the ground;
  * background block: inlet/outlet/top bcfarfield, ground bcsymmetryplane (ideal slip),
    span faces bcsymmetryplane; near-field O-mesh one cell wide, outer face overset.

Mesh family: L2 is the repo's baseline mesh (77k cells). A level with refinement factor f
divides every spacing (surface points, wall spacing, layers, background core) by f and takes
the background stretching ratio to the power 1/f, so the levels are geometrically similar.
Levels step by sqrt(2): L2 (f = 1), L15 (sqrt 2), L1 (2), L05 (2 sqrt 2). The coarser L3 (f = 1/2)
cannot be hole-cut: with the near-field thickness fixed, halving the layers squares pyHyp's
growth ratio, the near-field cells outgrow the background cells only ~4 layers inside the
overset boundary, and the fringe on the ground side is left without donors (28 orphans; a
two-cell-wide O-mesh gives 126).
"""

import json
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent

# repo design case (README / case_config.json)
CASE = dict(chord=0.7, clearance=0.2, speed_kmh=100.0, T=288.15, P=101325.0, span=0.07)

def level(f):
    """Repo baseline mesh (f = 1) refined by factor f in every direction."""
    return dict(f=f, nSurf=2 * int(round(136 * f)), nTE=int(round(16 * f)), N=int(round(80 * f)) + 1,
                s0=5.0e-6 / f, dx=0.01 / f, dy=0.005 / f, r=1.08 ** (1 / f))


LEVELS = {"L3": level(0.5), "L2": level(1.0), "L15": level(2**0.5), "L1": level(2.0), "L05": level(2**1.5)}
MARCH = 0.11  # near-field thickness, m (repo); stays below the 0.20 m clearance
NFFD = 10  # FFD chordwise stations


def flow_conditions():
    c = CASE
    V = c["speed_kmh"] / 3.6
    rho = c["P"] / (287.058 * c["T"])
    mu = 1.716e-5 * (c["T"] / 273.15) ** 1.5 * (273.15 + 110.4) / (c["T"] + 110.4)
    mach = V / np.sqrt(1.4 * 287.058 * c["T"])
    return dict(V=V, rho=rho, mu=mu, mach=float(mach), reynolds=float(rho * V * c["chord"] / mu))


def case_name(level, alpha):
    return f"c070_h020_a{alpha:+.1f}_{level}"


def section(level, alpha, chord=CASE["chord"], clearance=CASE["clearance"]):
    """Closed surface loop TE(upper) -> LE -> TE(lower) -> TE face, placed in the world frame.
    Returns body-frame coordinates, world-frame coordinates and the placed c/4 pivot."""
    lv = LEVELS[level]
    theta = np.linspace(0, 2 * np.pi, lv["nSurf"] + 1)
    x = (1 + np.cos(theta)) / 2
    yt = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    y = np.where(theta <= np.pi, yt, -yt)
    te = np.linspace(y[-1], y[0], lv["nTE"] + 1)[1:]
    xb = np.r_[x, np.ones(te.size)] * chord
    yb = np.r_[y, te] * chord
    xw, yw, pivot = place(xb, yb, alpha, chord, clearance)
    return (xb, yb), (xw, yw), pivot


def place(xb, yb, alpha, chord=CASE["chord"], clearance=CASE["clearance"], offset=None):
    """Body frame -> world frame: rotate by -alpha about (c/4, 0), then shift vertically so the
    lowest point of the section sits at `clearance` (or by a given offset)."""
    a = np.deg2rad(-alpha)
    pivot = 0.25 * chord
    xr = np.cos(a) * (xb - pivot) - np.sin(a) * yb + pivot
    yr = np.sin(a) * (xb - pivot) + np.cos(a) * yb
    if offset is None:
        sec = section_body(chord)
        offset = clearance - float(rotate_only(*sec, alpha, chord)[1].min())
    return xr, yr + offset, (pivot, offset)


def section_body(chord, n=2000):
    """Dense body-frame section, only used to find the lowest point consistently at all levels."""
    theta = np.linspace(0, 2 * np.pi, n + 1)
    x = (1 + np.cos(theta)) / 2
    yt = 0.6 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    return x * chord, np.where(theta <= np.pi, yt, -yt) * chord


def rotate_only(xb, yb, alpha, chord):
    a = np.deg2rad(-alpha)
    p = 0.25 * chord
    return np.cos(a) * (xb - p) - np.sin(a) * yb + p, np.sin(a) * (xb - p) + np.cos(a) * yb


def write_surface_xyz(path, xw, yw, span, nz=2):
    """Surface strip for pyHyp on nz spanwise planes (nz = 2: the repo's one-cell-wide O-mesh)."""
    xyz = np.empty((xw.size, nz, 1, 3))
    xyz[..., 0], xyz[..., 1] = xw[:, None, None], yw[:, None, None]
    xyz[..., 2] = np.linspace(0.0, span, nz)[None, :, None]
    with open(path, "w") as f:
        f.write(f"1\n{xw.size} {nz} 1\n")
        for d in range(3):
            np.savetxt(f, xyz[..., d].ravel(order="F"), fmt="%.15e")


def write_ffd(path, alpha, chord, span, offset, nffd=NFFD, xmargin=0.001, ymargin=0.02):
    """Fitted FFD box (nffd x 2 x 2), built in the body frame around the half-thickness
    distribution and placed with the same rotation/shift as the section."""
    xs = np.linspace(-xmargin, 1 + xmargin, nffd)
    xc = np.clip(xs, 0, 1)
    yt = 0.6 * (0.2969 * np.sqrt(xc) - 0.126 * xc - 0.3516 * xc**2 + 0.2843 * xc**3 - 0.1015 * xc**4)
    X = np.zeros((nffd, 2, 2, 3))
    for k, z in enumerate((-0.25 * span, 1.25 * span)):
        for j, sgn in enumerate((-1, 1)):
            xb, yb = xs * chord, sgn * (yt + ymargin) * chord
            xw, yw, _ = place(xb, yb, alpha, chord, offset=offset)
            X[:, j, k, 0], X[:, j, k, 1], X[:, j, k, 2] = xw, yw, z
    with open(path, "w") as f:
        f.write(f"1\n{nffd} 2 2\n")
        for d in range(3):
            np.savetxt(f, X[..., d].ravel(order="F"), fmt="%.15e")


def extend(start, stop, step, ratio):
    """Repo stretching: steps of `step` growing by `ratio` from start until stop is reached."""
    values = [start]
    direction = np.sign(stop - start)
    while abs(stop - values[-1]) > step:
        values.append(values[-1] + direction * step)
        step *= ratio
    values.append(stop)
    return np.array(values)


def background_axes(level, chord=CASE["chord"]):
    """Repo baseline background (core x in [-0.5c, 2c], y in [0, 0.65 m]; far field 10c upstream,
    15c downstream, 10c above), with the level's core spacing and stretching ratio."""
    lv = LEVELS[level]
    x0, x1, y1 = -0.5 * chord, 2.0 * chord, 0.65
    nx, ny = int(round((x1 - x0) / lv["dx"])), int(round(y1 / lv["dy"]))
    xx = np.r_[extend(x0, -10 * chord, lv["dx"], lv["r"])[::-1][:-1], np.linspace(x0, x1, nx + 1),
               extend(x1, 15 * chord, lv["dx"], lv["r"])[1:]]
    yy = np.r_[np.linspace(0, y1, ny + 1), extend(y1, 10 * chord, lv["dy"], lv["r"])[1:]]
    return xx, yy


# ---------------------------------------------------------------- ADflow
REPO_OPTIONS = {  # scripts/run_ground_case.py
    "equationType": "RANS", "equationMode": "steady", "liftIndex": 2, "turbulenceModel": "SA",
    "turbulenceOrder": "first order", "discretization": "central plus scalar dissipation",
    "useWallFunctions": False, "useANKSolver": True, "useNKSolver": True,
    "ANKSwitchTol": 10.0, "NKSwitchTol": 1e-6,
    "ANKCoupledSwitchTol": 1e-2, "ANKSecondOrdSwitchTol": 1e-2,
    "ANKADPC": True, "NKADPC": True, "NKViscPC": True,
    "NKSubspaceSize": 120, "NKPCILUFill": 3, "turbResScale": 10000.0,
    "MGCycle": "sg", "useZipperMesh": False, "nearWallDist": 0.02,
}

# Solver variants for the ANK test. All other options are the repo's.
SOLVERS = {
    "repo": {},  # ANK -> coupled ANK at 1e-2 -> NK at 1e-6
    "ank_uncoupled": {"ANKCoupledSwitchTol": 1e-16},  # segregated ANK all the way to NK
    "ank_off": {"useANKSolver": False},  # DADI smoother -> NK at 1e-6
    # NK stalls (steps 0.01-0.1) on the finer meshes after switching at 1e-6: let segregated ANK
    # carry the solution two decades deeper and give NK a larger Krylov space and stronger ILU
    "ank_deep": {"ANKCoupledSwitchTol": 1e-16, "NKSwitchTol": 1e-8, "NKSubspaceSize": 400,
                 "NKPCILUFill": 4, "NKASMOverlap": 3},
}

ADJOINT_OPTIONS = {
    "adjointSolver": "GMRES", "adjointL2Convergence": 1e-11, "ADPC": True, "adjointMaxIter": 5000,
    "adjointSubspaceSize": 200, "ILUFill": 3, "ASMOverlap": 3, "outerPreconIts": 3,
    "frozenTurbulence": False, "restartADjoint": False,
}

FUNCS = ("cl", "cd", "cdp", "cdv", "cmz")


def adflow_options(meta, outdir, solver="repo", l2=1e-12, ncycles=20000, overset="frozen"):
    o = dict(REPO_OPTIONS)
    o.update(SOLVERS[solver])
    o.update(ADJOINT_OPTIONS)
    o.update({
        "gridFile": meta["gridFile"], "outputDirectory": str(outdir),
        "L2Convergence": l2, "nCycles": ncycles,
        "monitorVariables": ["cpu", "resrho", "resturb", "cl", "cd"],
        "writeVolumeSolution": False, "writeSurfaceSolution": False,
        "oversetUpdateMode": overset,
    })
    return o


def aero_problem(meta, name="fwo"):
    from baseclasses import AeroProblem

    fc = flow_conditions()
    c = CASE
    return AeroProblem(name=name, mach=fc["mach"], reynolds=fc["reynolds"], reynoldsLength=c["chord"],
                       T=c["T"], alpha=0.0, areaRef=c["chord"] * c["span"], chordRef=c["chord"],
                       xRef=meta["pivot"][0], yRef=meta["pivot"][1], zRef=0.5 * c["span"])


def add_functions(solver, ap, funcs=FUNCS):
    names = [solver.addFunction(f, "airfoil", name=f"{f}_airfoil") for f in funcs]
    ap.evalFuncs = names
    return names


def load_meta(case_dir):
    return json.loads((Path(case_dir) / "case.json").read_text())
