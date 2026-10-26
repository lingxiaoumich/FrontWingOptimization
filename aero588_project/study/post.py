"""Tables and figures for the solver test, the mesh study and the adjoint/FD study.

    python study/post.py            (any Python with numpy + matplotlib; no MACH-Aero needed)

Reads study/results/{solve,adjoint}/*.json and the ADflow iteration tables in study/logs,
writes study/figs/*.png and study/results/summary.json, and prints markdown tables.
"""

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STUDY = Path(__file__).resolve().parent
FIGS = STUDY / "figs"
FIGS.mkdir(exist_ok=True)
LEVEL_ORDER = ["L3", "L2", "L15", "L1", "L05"]
SOLVER_LABEL = {"repo": "repo: ANK + coupled ANK (1e-2) + NK", "ank_uncoupled": "ANK, coupled phase off + NK",
                "ank_off": "ANK off: DADI + NK"}


def parse_log(path):
    """ADflow iteration table -> dict of columns (iter, iterTot, type, cfl, step, lin, cpu, totalRes, ...)."""
    cols, rows = None, []
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("#  Grid") and "Iter" in line:
            names = [c.strip().lower() for c in line.strip("# |\n").split("|")]
            # the three 'Iter' columns are iteration, total (linear) iterations and solver type
            renamed, seen = [], 0
            for c in names:
                if c == "iter":
                    c = ("iter", "itertot", "type")[seen]
                    seen += 1
                renamed.append({"res rho": "resrho", "res nuturb": "resturb", "c_lift": "cl", "c_drag": "cd"}.get(c, c))
            cols = renamed
            continue
        if cols and re.match(r"^\s+1\s+\d+\s+\d+\s+\S+", line):
            tok = line.split()
            if len(tok) == len(cols):
                rows.append(tok)
    if not rows:
        return {}
    out = {}
    for j, c in enumerate(cols):
        vals = [r[j] for r in rows]
        try:
            out[c] = np.array([float(v) if v != "----" else np.nan for v in vals])
        except ValueError:
            out[c] = vals
    return out


def load_solves():
    runs = []
    for f in sorted((STUDY / "results" / "solve").glob("*__*.json")):
        d = json.loads(f.read_text())
        d.pop("history", None)
        log = STUDY / "logs" / f"{d['case']}__{d['solver']}.log"
        d["log"] = parse_log(log) if log.exists() else {}
        runs.append(d)
    return runs


def time_to(run, tol):
    """Wall time (from the CPU-time monitor) and nonlinear iterations until totalRes/totalRes_0 <= tol."""
    h = run["log"]
    if not h or "totalres" not in h:
        return None, None
    r = h["totalres"] / h["totalres"][0]
    idx = np.where(r <= tol)[0]
    if not idx.size:
        return None, None
    i = int(idx[0])
    t = h.get("wall")
    return (float(t[i]) if t is not None else None), int(h["iter"][i])


def solver_table(runs):
    print("\n### Solver test\n")
    print("| case | cells | solver | converged | final res. red. | wall to 1e-8 (s) | wall to 1e-12 (s) | iters | C_L | C_D |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(runs, key=lambda r: (r["alpha"], LEVEL_ORDER.index(r["level"]), r["solver"])):
        t8, _ = time_to(r, 1e-8)
        t12, _ = time_to(r, 1e-12)
        it = int(r["log"]["iter"][-1]) if r["log"] else None
        ok = not (r["solve_failed"] or r["fatal_fail"])
        f = lambda v: f"{v:.0f}" if v is not None else "—"
        print(f"| {r['case']} | {r['nCells']} | {r['solver']} | {'yes' if ok else 'NO'} | {r['residual_reduction']:.1e} | "
              f"{f(t8)} | {f(t12)} | {it} | {r['cl']:.6f} | {r['cd']:.6f} |")
    # residual histories, one panel per case
    cases = sorted({r["case"] for r in runs}, key=lambda c: c)
    if not cases:
        return
    fig, axes = plt.subplots(1, len(cases), figsize=(4.2 * len(cases), 3.6), squeeze=False)
    for ax, c in zip(axes[0], cases):
        for r in [r for r in runs if r["case"] == c]:
            h = r["log"]
            if not h:
                continue
            t = h.get("wall")
            ax.semilogy(t, h["totalres"] / h["totalres"][0], label=r["solver"])
        ax.set_title(c, fontsize=9)
        ax.set_xlabel("wall time (s)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0][0].set_ylabel("total residual / initial")
    axes[0][0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "solver_test_residuals.png", dpi=150)


def richardson(h, f):
    """Three levels (coarse -> fine) with constant ratio: observed order, extrapolated value, GCI_fine."""
    h, f = np.asarray(h), np.asarray(f)
    r = h[0] / h[1]
    e21, e32 = f[2] - f[1], f[1] - f[0]
    if e21 == 0 or e32 / e21 <= 0:
        return dict(p=np.nan, f_ext=np.nan, gci_fine=np.nan, monotone=False)
    p = np.log(e32 / e21) / np.log(r)
    if p <= 0:  # differences grow with refinement: not in the asymptotic range, no extrapolation
        return dict(p=float(p), f_ext=np.nan, gci_fine=np.nan, monotone=False)
    f_ext = f[2] + e21 / (r**p - 1)
    gci = 1.25 * abs(e21 / f[2]) / (r**p - 1)
    return dict(p=float(p), f_ext=float(f_ext), gci_fine=float(gci), monotone=True)


def mesh_table(runs, solver_pref=("ank_uncoupled", "repo", "ank_off")):
    print("\n### Mesh study\n")
    summary = {}
    for a in sorted({r["alpha"] for r in runs}):
        best = {}
        for r in runs:
            if r["alpha"] != a or r["solve_failed"] or r["fatal_fail"]:
                continue
            k = r["level"]
            if k not in best or solver_pref.index(r["solver"]) < solver_pref.index(best[k]["solver"]):
                best[k] = r
        lv = [l for l in LEVEL_ORDER if l in best]
        if not lv:
            continue
        n = np.array([best[l]["nCells"] for l in lv], float)
        h = np.sqrt(n[0] / n)  # relative spacing, coarsest = 1
        print(f"\nα = {a:+.1f}°\n")
        print("| level | cells | h/h_L2 | C_L | C_D | C_Dp | C_Dv | downforce/drag | solver |")
        print("|---|---|---|---|---|---|---|---|---|")
        for l, hh in zip(lv, h):
            r = best[l]
            print(f"| {l} | {r['nCells']} | {hh:.3f} | {r['cl']:.6f} | {r['cd']:.6f} | {r['cdp']:.6f} | {r['cdv']:.6f} | "
                  f"{r['downforce_to_drag']:.3f} | {r['solver']} |")
        res = {}
        if len(lv) >= 3:
            print("\n| quantity | observed order p | extrapolated | GCI (finest) | error of L2 vs extrapolated |")
            print("|---|---|---|---|---|")
            for q in ("cl", "cd", "cdp", "cdv", "downforce_to_drag"):
                f = [best[l][q] for l in lv]
                rr = richardson(h[-3:], f[-3:])
                err_l2 = abs(best[lv[0]][q] - rr["f_ext"]) / abs(rr["f_ext"]) if rr["monotone"] else np.nan
                res[q] = dict(levels=lv, values=f, **rr, err_L2=err_l2)
                print(f"| {q} | {rr['p']:.2f} | {rr['f_ext']:.6f} | {100 * rr['gci_fine']:.2f} % | {100 * err_l2:.2f} % |")
        summary[f"{a:+.1f}"] = res
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
        for ax, q, lab in zip(axes, ("cl", "cd", "downforce_to_drag"), ("$C_L$", "$C_D$", "downforce / drag")):
            ax.plot(h, [best[l][q] for l in lv], "o-")
            for l, hh in zip(lv, h):
                ax.annotate(l, (hh, best[l][q]), textcoords="offset points", xytext=(4, 4), fontsize=8)
            if q in res and res[q]["monotone"]:
                ax.axhline(res[q]["f_ext"], ls="--", c="gray", lw=1, label="Richardson")
                ax.legend(fontsize=8)
            ax.set_xlabel("h / h$_{L2}$")
            ax.set_title(lab)
            ax.set_xlim(left=0)
            ax.grid(alpha=0.3)
        fig.suptitle(f"Mesh study, α = {a:+.1f}° (c = 0.7 m, h = 0.2 m, 100 km/h)", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIGS / f"mesh_study_a{a:+.1f}.png", dpi=150)
    return summary


def fd_without_jump(rs):
    """Central differences across a discontinuity of size delta at the baseline read
    FD(h) = g + delta/(2h) + b h^2. Least squares over the steps gives the derivative g of the
    smooth branch (what the adjoint linearizes) and the jump delta."""
    h = np.array([r["step"] for r in rs])
    A = np.column_stack([np.ones_like(h), 1 / (2 * h), h**2])
    g, delta, _ = np.linalg.lstsq(A, np.array([r["fd"] for r in rs]), rcond=None)[0]
    return float(g), float(delta)


def adjoint_tables():
    print("\n### Adjoint vs finite differences\n")
    out = {}
    for f in sorted((STUDY / "results" / "adjoint").glob("*.json")):
        d = json.loads(f.read_text())
        out[f.stem] = d
        print(f"\n{f.stem}: flow {d['flow_s']:.0f} s, {d['nAdjoints']} adjoints {d['adjoint_s']:.0f} s on {d['nProcs']} ranks; "
              f"checks {json.dumps({k: v for k, v in d['checks'].items() if k != 'dvmap'})}\n")
        print("| variable | function | adjoint | FD (best step) | best step | rel. diff (best) | FD, jump removed | rel. diff | jump δ | rel. diff at each step |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        rows = d["rows"]
        keys = []
        for r in rows:
            if (r["variable"], r["function"]) not in keys:
                keys.append((r["variable"], r["function"]))
        for v, fn in keys:
            rs = [r for r in rows if r["variable"] == v and r["function"] == fn]
            b = min(rs, key=lambda r: r["rel_diff"])
            steps = ", ".join(f"{r['step']:.0e}: {r['rel_diff']:.1e}" for r in rs)
            g, delta = fd_without_jump(rs)
            rel = abs(b["adjoint"] - g) / abs(g)
            print(f"| {v} | {fn} | {b['adjoint']:.6e} | {b['fd']:.6e} | {b['step']:.0e} | {b['rel_diff']:.1e} | "
                  f"{g:.6e} | {rel:.1e} | {delta:+.1e} | {steps} |")
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), sharey=True)
        for ax, fn in zip(axes, ("cl_airfoil", "cd_airfoil", "ratio")):
            for v in dict.fromkeys(r["variable"] for r in rows):
                rs = [r for r in rows if r["variable"] == v and r["function"] == fn]
                st = np.array([r["step"] for r in rs])
                st = st / 100 if v == "incidence" else st  # plot incidence on the same 1e-3..1e-6 scale
                ax.loglog(st, [max(r["rel_diff"], 1e-16) for r in rs], "o-", label=v)
            ax.set_title({"cl_airfoil": "$C_L$", "cd_airfoil": "$C_D$", "ratio": "$C_L/C_D$"}[fn])
            ax.set_xlabel("FD step (m; incidence: deg/100)")
            ax.invert_xaxis()
            ax.grid(True, which="both", alpha=0.3)
        axes[0].set_ylabel("|adjoint − FD| / |FD|")
        axes[0].legend(fontsize=7)
        fig.suptitle(f.stem, fontsize=9)
        fig.tight_layout()
        fig.savefig(FIGS / f"adjoint_fd_{f.stem}.png", dpi=150)
    return out


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8")  # Windows consoles default to cp1252
    runs = load_solves()
    solver_table(runs)
    ms = mesh_table(runs)
    adjoint_tables()
    (STUDY / "results" / "summary.json").write_text(json.dumps({"mesh_study": ms}, indent=1, default=float) + "\n")
