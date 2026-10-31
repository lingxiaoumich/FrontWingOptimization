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

## Single-step test with `ank_deep` (`step0.2_ank_deep/`)

The exact trial shape that stalled under the main settings (`xmax_0.2/incoming_design.json`: ≤1.07 mm of control-point movement, −0.105° incidence) was re-solved from the same converged restart flow with `--solver ank_deep` (NK hand-off at 1e-8 instead of 1e-4, Krylov subspace 400, ILU fill 4) and live residual output.

- **Converged in 377 s** to a residual ratio of 8.6e-11: 132 segregated-ANK iterations from a starting residual of 1.0e3 down to 1.7e-2, then **2 NK iterations** to 1e-10.
- Under the main settings the NK hand-off fires at 1.8e2 absolute, which this trace crosses at ANK iteration 3, while C_L is still swinging between −1.54 and −1.58. The flow is four decades short of where Newton converges; that is the stall.
- C_L·c = −0.77936 m, C_D·c = 0.014451 m. The 0.05 control step along the same direction gave ΔC_L·c = −0.00900; four times that predicts −0.78056. The 0.2 step delivered 97 % of the linear extrapolation, so the aerodynamics are linear across it and the optimizer would have accepted it.
- Mesh cleared: 0 orphans, no negative-volume or warping warnings.
- Cost per unit progress: control = 4 × (≈50 s flow + ≈130 s adjoints) ≈ 720 s; this = 377 s + 130 s ≈ 510 s.

## Next step

Re-run the sweep (0.05, 0.2, 0.5, 1.0) with `--solver ank_deep` throughout, so the comparison measures PSQP's step limit rather than the flow solver. The main-settings arm of the 0.2 test was not re-run with live output, so the hand-off diagnosis above is read off the `ank_deep` trace rather than observed directly.

## Reproduction

```sh
# inside the container, aero588_project mounted at /project
CHECKS=/project/optimizations/xmax_sweep_20261031/failed_ncycles4000_xmax_0.05/source/aerodynamic_derivative_check.json \
NCYCLES=20000 bash /project/scripts/run_xmax_sweep.sh 10 0.05 0.2 0.5 1.0
python scripts/plot_xmax_sweep.py --sweep optimizations/xmax_sweep_20261031
```

Single-step test: `bash /project/scripts/run_step_test.sh <design.json> ank_deep <run name> [minutes]`.

The sweep uses four driver options added on branch `opt1031` (defaults reproduce the main setup): `--xmax`, `--ncycles`, `--solver`, and `--restart-flow`, which starts from the converged October 29 final flow at the restart design instead of solving the NACA baseline first. Flow solutions (`*.cgns`, `*.plt`) and `optimization.hst` stay local and are not committed.
