# Flow and adjoint solver settings: speed against accuracy, October 31

Under `ank_deep`, one optimizer iteration took 4–5 min at 32 ranks, and about 80 % of that was the flow solve. In the profiled solves, about three quarters of the flow iterations were spent in the slow ANK phase, because `ank_deep` hands over to Newton–Krylov only at `NKSwitchTol = 1e-8`. This test re-ran three real steps from the October 31 continuation under six settings. Each run measures the flow solve and both adjoints, and compares lift, drag and the full 21-variable gradients against the current setting.

Steps (scaled length; incidence in degrees, shape in cm), each solved from the converged flow at its start design:

| Case | Start → step design | Length |
|---|---|---:|
| S | 2 mm chain XMAX 2.0, stage 00 eval 0 → eval 3 | 1.50 |
| M | 2 mm chain XMAX 5.0, stage 00 eval 0 → eval 3 | 2.31 |
| L | 2 mm chain XMAX 2.0, stage 00 eval 0 → eval 4 | 3.38 (beyond any accepted step) |

## Results (32 ranks, flow solve + adjoint pair)

| Setting | Per iteration | Speed-up | CL, CD vs reference | Gradient vs reference | Converged |
|---|---:|---:|---:|---:|---|
| reference: `ank_deep` (NK switch 1e-8, L2 1e-10, adjoint 1e-7) | 252–300 s | 1.0× | — | — | 3/3 |
| **NK switch 1e-6** | **92–107 s** | **2.7–3.0×** | ≤ 1.1e-7 | ≤ 3.1e-6 | **3/3** |
| NK switch 1e-5 | 78–91 s | 2.8–3.7× | ≤ 1.4e-7 | ≤ 1.5e-6 | **2/3** (S stalled; 491 s, failed) |
| L2 convergence 1e-8 | 185–206 s | 1.3–1.6× | ≤ 6.3e-6 | ≤ 1.3e-5 | 3/3 |
| NK switch 1e-6 + L2 1e-8 | 79–99 s | 3.0–3.2× | ≤ 5.1e-5 | ≤ 2.2e-4 | 3/3 |
| … + adjoint 1e-6 | 82–99 s | 3.0–3.3× | ≤ 5.1e-5 | ≤ 2.2e-4 | 3/3 |

Full table: `comparison.txt` (per case) and `comparison.json`. Relative errors are taken against the reference; the gradient error is the norm of the difference over the norm of the reference gradient.

- **An earlier NK handover is where the time is.** At 1e-6 the flow solve drops from 206–255 s to 45–55 s, and the answers do not change at any level the optimizer can see. At 1e-5 it is faster again, but it stalls on one step in three. The solver study made `ank_deep` to avoid exactly that kind of stall.
- **Loosening convergence buys little.** L2 1e-8 alone is 1.3–1.6× faster. On top of the 1e-6 handover it adds about 10 % speed and makes the errors about 100× larger. The adjoint (about 45 s for both functions) is no faster at 1e-6.
- For scale: the gradient's own linear-model error on steps of length 2 is 2–17 % (actual change against adjoint prediction). Convergence errors of 1e-6 are irrelevant to the optimizer at this stage. They would matter only near a converged optimum, where steps change CL by about 1e-5.

## Fallback for stalled solves

`run_slip_optimization.py --nk-fallback-tol 1e-8`: if a flow solve stalls (`solveFailed` without `fatalFail`), the same design is re-solved once from the last converged state with the later NK handover. Only after that is a failure reported. It was checked on the case that stalls (S with NK switch 1e-5). The stall took 491 s, the retry then converged, and the result matches the reference to 7e-9 in CL and 2e-7 in the gradient (`runs/S_nk5_fallback8`). A stall therefore costs about 12 min instead of a backtracked step. The record's `nk_fallback` field marks evaluations that needed it.

## Adopted

From the evening of October 31, the 2 mm continuation chains run with `--nk-switch-tol 1e-6 --nk-fallback-tol 1e-8`. L2 and adjoint convergence keep their defaults (1e-10, 1e-7). Expected iteration time is about 1.7 min instead of 4.5 min.

## Reproduction

```sh
source ~/software/mdolab/env.sh
# one run (see scripts/test_solver_settings.py); cases and settings as in the tables above
mpirun -np 32 python scripts/test_solver_settings.py --step-to designs/S_to.json --root runs/S_nk6 \
    --ffd-stations 5 --shape-directions xy --drag-cap 0.1 --min-local-thickness 0.002 --solver ank_deep \
    --ncycles 20000 --reuse-checks --restart-design designs/S_from.json \
    --restart-flow ../continuation_tmin2mm_xmax2.0_20261031/stage_00/baseline/flow_000_vol.cgns --nk-switch-tol 1e-6
python scripts/compare_solver_settings.py optimizations/solver_tolerance_test_20261031
```

Each run needs `source/` staged from `setups/slip_ground_xy5_cd005/source` plus the derivative-check JSON, and an OpenMPI rankfile for core pinning. Flow files are not committed.
