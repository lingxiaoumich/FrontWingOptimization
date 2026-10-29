# Racing-car front-wing section optimization

This project studies a NACA 0012 airfoil section in ground effect using MACH-Aero: ADflow for steady RANS and discrete adjoints, pyGeo for free-form deformation (FFD), IDWarp for mesh deformation, and pyOptSparse for optimization. It is a quasi-two-dimensional section study; finite-span effects, endplates, wheels, and the vehicle body are excluded.

The current optimization uses a 0.50 m initial chord, 200 km/h horizontal freestream, an initial incidence of −3°, and ideal slip ground. It minimizes upward-positive CL to maximize downforce. Incidence is free. The latest FFD lattice has five upper/lower pairs along the chord, extruded across two tied span planes: 20 physical control points, 10 x-displacement variables, 10 y-displacement variables, and one incidence variable.

## Completed optimization runs

Both runs completed 50 PSQP major iterations at the iteration limit, **before meeting the gradient stopping criterion**. The results are feasible improvements, not converged optima. Both use the original 88,548-cell mesh; grid independence has not been established.

| Run | FFD shape variables | CD cap | Final CL | Final CD | Downforce increase | Independent full-connectivity check |
|---|---|---:|---:|---:|---:|---:|
| [October 27](optimizations/slip_ground_box_20261027/README.md) | 20 y variables, 10 upper/lower pairs | Initial CD ≈ 0.0229614 | −1.32757 | 0.0229556 | 15.26% | 15.13% |
| [October 29](optimizations/slip_ground_xy5_cd005_20261029/README.md) | 10 x + 10 y variables, 5 upper/lower pairs | 0.05 | −1.48912 | 0.0269161 | 29.28% | 29.05% |

Each run also has one incidence variable. The common baseline has CL ≈ −1.15181 and CD ≈ 0.0229614. Coefficients use the fixed initial chord and reference area 0.025 m². Since the parameterization and drag cap changed together, this comparison does not isolate the benefit of x movement.

The latest run took 44.27 minutes in the optimizer and 45.96 minutes including baseline checks and the final re-solve, using four MPI ranks. All 57 CFD evaluations converged; all 52 optimizer/final evaluations were feasible. One preoptimization finite-difference perturbation slightly crosses the design clearance constraint and is identified in `final_checks.json`.

![Latest geometry and FFD comparison](optimizations/slip_ground_xy5_cd005_20261029/plots/ffd_before_after_comparison.png)

## Current constraints and geometry definition

- Drag coefficient CD ≤ 0.05.
- Actual CFD surface coordinates remain within −0.125 ≤ x ≤ 0.625 m and 0.05 ≤ y ≤ 0.35 m. The box is 1.5 initial chords wide and 0.6 initial chords tall; its bottom enforces a clearance of 0.1 initial chords.
- Maximum section thickness is at least 0.06 initial chords = 0.03 m. This is a lower bound on the **thickest section**, not a local thickness requirement everywhere.
- The box constrains the actual airfoil surface, rather than the FFD control points. There are no explicit incidence or x/y displacement bounds and no leading/trailing-edge fixing equalities. Invalid boundary topology and failed CFD trials are rejected as numerical failures.

With x movement enabled, corresponding upper/lower material points need not share a chordwise location. The current [thickness implementation](scripts/airfoil_geometry.py) therefore intersects chord-normal lines with the actual closed, piecewise-linear CFD surface and finds the longest continuous interior segment. Exterior gaps in concave sections are excluded. Its active-feature derivative is chained through pyGeo; a maximum can have a derivative kink when the active location changes. The original leading/trailing-edge markers define the deformed chord direction. Reference dimensions remain tied to the initial 0.5 m chord.

## Code, inputs, and results

- [Latest prepared setup](setups/slip_ground_xy5_cd005/README.md): parameterization, inputs, units, geometry tests, and reproduction instructions.
- [Optimization driver](scripts/run_slip_optimization.py): defaults to five FFD pairs, global x/y movement, CD cap 0.05, PSQP, and 50 major iterations.
- [Latest results folder](optimizations/slip_ground_xy5_cd005_20261029/README.md): settings, comparisons, source snapshots, optimizer history, mesh/flow solutions, and independent final CFD validation.
- [Latest comparison plots](optimizations/slip_ground_xy5_cd005_20261029/plots): mesh/FFD overlays, x/y control-point movement, airfoil geometry, Cp curves/contours, total pressure/loss, loads, and optimization history, in PNG and SVG.
- [Original slip-ground CFD case](cases/naca0012_c050_h005_u200_am3/README.md): NACA 0012 baseline at the current dimensions and operating conditions.

The optimization deforms the whole overset mesh with USMesh and uses frozen interpolation connectivity/weights for consistent primal/adjoint evaluations. Each selected final design was also evaluated with rebuilt overset connectivity. That check is not a grid-independence study. Solver residuals, convergence flags, surface constraints, grid-corner orientations, and overset orphans are documented with the saved results. Small numerical total-pressure overshoots are retained in the plots.

For reproduction, use a fresh run directory and the input/source snapshot for that run. The latest run used `mdolab/public@sha256:bfb64541108e7f819b7abe0f2f119ca08f1bbd7a5ea9ebe7377e90ad8f0d9333`. Its [execution record](optimizations/slip_ground_xy5_cd005_20261029/source/execution.json) contains the command and environment. Docker uses local staging because direct Google Drive mounts are unreliable. Do not hot-start the new x/y configuration from the old y-only optimizer history; their variable groups and shape spaces differ.

## Moving no-slip ground extension

A separate [ADflow moving-wall extension](solver_extensions/moving_wall/README.md) prescribes constant physical wall velocity while retaining the viscous no-slip boundary condition. It includes the exact solver patch, a Docker build, MPI checks, and real-valued spatial-derivative validation.

The [moving-ground airfoil case](cases/naca0012_c050_h005_u200_am3_moving_ground/README.md) is a converged boundary-condition validation, not an optimization result. The road grid still requires refinement/convergence work. Its supplemental 15 micrometre refinement attempt is explicitly **unconverged** and must not be used as a final aerodynamic result. The two optimization runs above use slip ground.

## Earlier project material

The top-level `case_config.json`, `geometry/`, `meshes/`, `results/`, and `runs/baseline_a0/` preserve the original 0.70 m chord, 100 km/h, 0.20 m clearance baseline. Those legacy settings are separate from the current 0.50 m / 200 km/h optimization. Proposal sources and graphics are under `output/pdf/` and `proposal latex/`; setup attempts and numerical-recovery history are retained with their run records.

MACH-Aero documentation: <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/>.
