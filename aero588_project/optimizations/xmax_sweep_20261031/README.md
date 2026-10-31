# PSQP step-limit (XMAX) sweep, October 31

**Stopped early, by decision, after the first step size beyond the default.** The question was whether larger PSQP steps make faster progress than the October 29 run, where every one of the 50 steps hit the limit `XMAX = 0.05` (scaled variables: incidence in degrees, shape displacements in cm). The answer so far: with the main solver settings, **the flow solver, not PSQP, is what limits the step size**.

Forces are reported as C_L·c = L′/q∞ and C_D·c = D′/q∞ (m), so the reference chord drops out.

| Run | Result |
|---|---|
| `xmax_0.05/` (control) | 10 iterations, all steps exactly at the limit, all solves converged in 40–60 s. C_L·c −0.74456 → −0.80357 m (+7.9 % downforce), C_D·c 0.01346 → 0.01471 m (59 % of the 0.025 m cap). Chord 0.511 m, clearance 0.0529 m, max thickness 0.0561 m, chord incidence −4.01°. |
| `xmax_0.2/` | The first trial step (length 0.2, i.e. up to 2 mm of control-point movement or 0.2° of incidence) was still solving after 30.7 min, against under a minute for a 0.05 step. Stopped without an accepted iterate. |
| 0.5, 1.0 | Not run. |
| `failed_ncycles4000_xmax_0.05/` | First attempt. The driver re-solves the NACA baseline before restarting, then jumped from the baseline flow to the restart design; with `nCycles = 4000` (which counts NK linear iterations too) that solve stopped at a 9e-7 residual ratio after 562 s, was flagged failed, and PSQP quit at iteration 0. Kept as a record. |

Other observations:

- The restart resets PSQP's Hessian approximation. Its first step points in a different direction from the last October 29 steps (cosine 0.24) and gains 0.0090 m per step at the start, against 0.0021 m at the end of October 29. After about five steps the direction locks in again (successive cosines ≥ 0.997).
- Incidence takes 53 % of the first restart step and 24 % of the tenth.
- This laptop (WSL2 + Docker, 4 ranks) is about 3× slower than the Mac that ran October 29: about 50 s per flow solve and 130 s for the two adjoints, against 17 s and 42 s.
- The baseline reproduces the partner's to 1e-13 in C_L with ADflow 2.13.1 (theirs 2.12.2), so results are directly comparable.

## Next step, if the sweep is resumed

Our solver study (branch `study/solver-mesh-adjoint`) diagnosed this pattern as a Newton–Krylov stall after an early hand-off (`NKSwitchTol = 1e-4` in the driver). It found `ank_uncoupled` fastest on this mesh, and `ank_deep` (NK at 1e-8, Krylov subspace 400, ILU fill 4) robust everywhere. Re-running all XMAX values, 0.05 included, with one of those would separate the solver limit from the optimizer's step limit.

## Reproduction

```sh
# inside the container, aero588_project mounted at /project
CHECKS=/project/optimizations/xmax_sweep_20261031/failed_ncycles4000_xmax_0.05/source/aerodynamic_derivative_check.json \
NCYCLES=20000 bash /project/scripts/run_xmax_sweep.sh 10 0.05 0.2 0.5 1.0
python scripts/plot_xmax_sweep.py --sweep optimizations/xmax_sweep_20261031
```

The sweep uses three driver options added on branch `opt1031` (defaults reproduce the main setup): `--xmax`, `--ncycles`, and `--restart-flow`, which starts from the converged October 29 final flow at the restart design instead of solving the NACA baseline first. Flow solutions (`*.cgns`, `*.plt`) and `optimization.hst` stay local and are not committed.
