# Moving no-slip ground for ADflow 2.12.2

This local extension prescribes a constant, **physical** velocity in metres per
second on selected viscous wall families. The original MDOLab image is unchanged.
The patched image is `aero588/adflow-moving-wall:2.12.2`.

## What changed

- `prescribedWallMotion.F90` is a new Fortran module storing velocities by family
  ID. It also reads the actual internal wall-velocity data for diagnostics.
- `patch_adflow.py` applies a version-checked patch to the installed ADflow source
  tree and writes an exact unified diff, saved here as `adflow-moving-wall.patch`.
- `Dockerfile` builds the real solver and installs the rebuilt Python package in
  a separate Docker image. Rebuilding the source alone is insufficient because
  the base image imports its installed package, rather than the repository copy.
- `validate_moving_wall.py` runs the moving-floor test, stationary-floor control,
  coordinate-Jacobian forward/reverse and finite-difference checks, and verifies
  that speed changes and solver refreshes preserve the prescribed velocity.
- The project script `scripts/run_ground_case.py` now accepts `ground_model`
  values `slip`, `stationary_no_slip`, and `moving_no_slip`. Existing case files
  without this key retain the slip-ground setup. Viscous-ground meshes use
  stretched spacing near the road, configured independently of airfoil spacing.

The exact upstream files added or modified are:

1. `src/modules/prescribedWallMotion.F90` (new)
2. `src/build/fileList`
3. `src/f2py/adflow.pyf`
4. `src/solver/solverUtils.F90`
5. `src/adjoint/outputForward/solverUtils_d.f90`
6. `src/adjoint/outputReverse/solverUtils_b.f90`
7. `src/adjoint/outputReverseFast/solverUtils_fast_b.f90`
8. `adflow/pyADflow.py`

## How the condition works

ADflow's existing viscous-wall boundary computes halo velocities as
`u_halo = -u_interior + 2*u_wall`. The extension supplies the existing `uSlip`
wall-velocity field; it does not replace the viscous flux or turbulence wall BC.
Despite its historical name, `uSlip` is the solid-wall velocity used by the
no-slip condition, not a free-slip boundary condition.

For a prescribed family, the finest-level wall values are overwritten after
each wall-motion update with `velocity_mps / uRef`. ADflow's existing restriction
propagates them to coarse multigrid levels. Other families and grid velocities
are unchanged. The forward tangent is zeroed on prescribed walls; reverse seeds
are discarded there before propagating the original wall-motion derivatives.
Every generated primal copy is patched as well. The prescribed speed is fixed
data rather than an optimization variable.

## Usage

The CGNS road boundary must be viscous (`bcwallviscous`), not symmetry. In Python:

```python
solver = ADFLOW(options=options, comm=comm)
solver.setWallVelocity([200.0 / 3.6, 0.0, 0.0], groupName="ground")
solver(aero_problem)
print(solver.getWallVelocityStats("ground"))
print(solver.getWallVelocityStats("airfoil"))
```

Use `AeroProblem(V=200.0/3.6, P=101325.0, T=288.15, ...)` to match the inlet
speed exactly. Independently computing Mach with a rounded gas constant can
produce a small mismatch between the intended wall speed and actual inflow.
Keep the freestream horizontal and rotate the airfoil for incidence changes.
Evaluate lift/drag on the airfoil family only, excluding ground forces.
The standard `cl`/`cd` iteration monitors sum all walls, including the road;
the case result uses the separately registered airfoil-only functions instead.

Call the setter on all MPI ranks with identical data. Diagnostic min/max speeds
are MPI-global and dimensional. The diagnostic entry count includes wall halos;
it is not the count of unique physical mesh faces.

Build from the project root:

```sh
docker build -t aero588/adflow-moving-wall:2.12.2 solver_extensions/moving_wall
```

For the supplied separate airfoil case, run `build` and then `solve` using
`scripts/run_ground_case.py --case-root <case-directory>` in that image. Use a
local staging directory for Docker mounts; copy results back to Google Drive.

## Supported scope

- Steady, real-valued ADflow 2.12.2, including MPI and its real discrete-adjoint
  spatial derivatives. The complex-step library has not been rebuilt.
- Constant tangential velocity on a stationary plane, with fixed freestream
  conditions. The caller must ensure tangency; this is not a moving/deforming
  mesh implementation.
- Shape derivatives treat the wall velocity as constant. Derivatives of a
  speed-matching condition with respect to flight speed are not implemented.
- Forward/reverse sources are patched explicitly. If Tapenade is rerun or the
  upstream solver changes, regenerate/reapply this extension and repeat tests.

## Validation

The independent tests and logs are saved under `validation/`.
The separate `validation/airfoil_100um` records are the converged airfoil boundary
check, with 100 micrometre road spacing and 2 micrometre airfoil spacing. Its road
y-plus exceeded one. A supplemental 15 micrometre road-grid check improved wall
resolution but did not converge; it is labeled separately and is not a final
aerodynamic solution. The supplied main case is the converged 100 micrometre run.

The moving empty tunnel starts at a mean-flow residual of order `1e-13`. A
relative reduction of eight orders cannot be obtained from that roundoff level;
ADflow may therefore set `solveFailed` even though the uniform solution is
already stationary to machine precision. This specific test accepts an absolute
residual below `1e-10`, a relative velocity error below `1e-6`, and negligible
drag, independently of that flag. The stationary control must converge normally.
This exception is **not** used for the airfoil case's convergence acceptance.

The saved airfoil mesh/solver checks remain preliminary; they do not establish
grid independence or validate an optimized design. See the implementation report
for the actual airfoil convergence status and wall-resolution checks.
