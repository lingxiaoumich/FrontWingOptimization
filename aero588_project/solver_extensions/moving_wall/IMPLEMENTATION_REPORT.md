# Moving-ground implementation report

Implemented on October 26, 2026, in an isolated `aero588/adflow-moving-wall:2.12.2` image. The original MDOLab image, previous CFD results, and proposal are preserved.

## Confirmed validation

- Empty moving-floor test: 55.5555555556 m/s wall and inlet, maximum fluid-velocity error 1.421e-14 m/s, drag coefficient -2.134e-18. The residual starts and remains at roundoff; see README for the absolute acceptance criterion and solver flag.
- Stationary no-slip control: normally converged, positive drag coefficient 0.00097827, maximum velocity deficit 54.5227 m/s.
- MPI tests use two ranks. Setting another wall speed and refreshing velocities preserves the requested value.
- Coordinate-Jacobian forward/reverse dot-product relative error 3.060e-14. Finite-difference relative error 3.883e-06. These validate the fixed-speed spatial derivative path, not speed-design-variable derivatives.

## Exact changes

The full upstream diff is in `adflow-moving-wall.patch`. `prescribedWallMotion.F90` stores physical wall velocities by family. `patch_adflow.py` adds the Python `setWallVelocity` and `getWallVelocityStats` APIs, F2PY bindings, and overrides in flow/forward/reverse/fast-primal wall-motion routines. `Dockerfile` rebuilds and installs the real library. `validate_moving_wall.py` and `inspect_wall_resolution.py` provide reproducible boundary and mesh checks.

The project `scripts/run_ground_case.py` now selects a viscous road, generates road-normal stretching, prescribes the speed, uses an exact velocity/pressure/temperature inflow, restricts slices and forces to the airfoil, and records the actual ground/airfoil speeds. The separate moving-ground case stores its configuration and results; no optimization has been run.

## Usage

```python
solver.setWallVelocity([200.0 / 3.6, 0.0, 0.0], groupName="ground")
solver(aero_problem)
print(solver.getWallVelocityStats("ground"))
```

The CGNS ground must be `bcwallviscous`; the setter alone does not convert a symmetry plane into a viscous wall.

## Limits

Supported: steady, real, fixed-speed, tangential wall motion on a stationary plane. The wall speed is prescribed data. The complex-step library and derivatives with respect to the wall/inflow speed have not been implemented. Shape/mesh derivatives are checked above. Grid independence and full optimization derivative validation remain separate tasks.

## Airfoil check

A separate 0.50 m NACA 0012 case at 200 km/h, 5 cm clearance and -3 degrees converged in 189.1 s on four MPI ranks, with residual reduction 9.9000e-09. The actual ground boundary speed is (55.5555555556,0,0) m/s and the airfoil speed is exactly (0,0,0). Airfoil-only CL=-1.14820011, CD=0.02360433.

The converged road grid has ground maximum y-plus 5.454 and airfoil maximum y-plus 0.470. This establishes implementation functionality, not final wall/grid resolution. A 15 micrometre road-grid check brought provisional ground maximum y-plus to 1.088, but did not converge to the requested tolerance. It is saved separately as an unconverged refinement attempt; a stronger-NK continuation was stopped and is not represented as a completed solution.

Delivered case: ../../cases/naca0012_c050_h005_u200_am3_moving_ground/ . Its refinement_check_15um subdirectory must not be used as a converged aerodynamic result. Further road-grid/solver-convergence work is needed before optimization.

`scripts/plot_current_run.py` also labels moving and stationary no-slip cases correctly rather than labeling every ground as slip.
