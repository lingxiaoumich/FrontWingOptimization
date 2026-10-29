# Slip-ground wing optimization

The 50-iteration PSQP run completed at the iteration limit (exit 11), before the gradient stopping criterion. The selected design is feasible and improves downforce by 15.2591% in the optimization model. An independent full-connectivity CFD check also converges, satisfies the constraints, and gives a 15.1260% increase. These are preliminary results on one mesh, not a converged optimum.

This folder contains an adjoint-based PSQP optimization after SLSQP setup trials of the quasi-2D front-wing section. The objective is minimum upward-positive CL (maximum downforce), subject to the initial wing's drag coefficient, a surface containment box, and a lower limit on maximum section thickness.

- Reference chord: 0.5 m; speed: 200 km/h; extrusion span: 0.05 m.
- Ground: ideal slip/symmetry at y = 0. The airfoil is a stationary no-slip wall.
- Surface box: −0.125 ≤ x ≤ 0.625 m and 0.05 ≤ y ≤ 0.35 m. This is 1.5 times the initial chord wide and 0.6 times the initial chord tall. Its lower edge enforces 0.1 times the initial chord in ground clearance.
- Maximum section thickness ≥ 0.06 times the initial chord = 0.03 m. The maximum is evaluated over paired upper/lower material stations in the chord-normal direction.
- Twenty FFD shape variables plus one incidence variable. The two spanwise FFD planes are tied. There is no explicit incidence bound, no shape-displacement bound, and no leading/trailing-edge fixing equality.
- The reference chord and coefficient reference area stay fixed even if the actual chord changes.

## Files

`settings.txt` records CFD, mesh, geometry, and optimizer settings; `settings.json` is the machine-readable version. `comparison_results.txt` compares the initial and chosen optimized design; `comparison_results.json` contains their detailed results.

`baseline/` contains the initial solution. `final/` contains the optimizer's returned solution. `best_feasible/` contains the lowest-CL converged evaluation satisfying all requested bounds within the documented numerical feasibility tolerance. The comparison report identifies which folder it uses.

Each solution folder contains the design variables, actual surface coordinates, volume and surface CGNS solution files, extracted plotting arrays, Cp/Cf/y+ CSV data, and a result JSON. `plots/` contains PNG and SVG comparisons of geometry, mesh, surface Cp, static pressure contours, total pressure, pressure loss, performance, and evaluation history.

`plots/ffd_mesh_comparison` shows the baseline and final FFD lattices over their actual CFD meshes. `plots/ffd_before_after_comparison` overlays the airfoils and lattices and plots the control-point height changes. Both figures are available as PNG and SVG. The 40 three-dimensional control points project into 20 visible spanwise pairs. Mesh lines are subsampled for clarity; airfoil outlines use all saved surface vertices. These figures use the final design after 50 iterations, before optimizer convergence.

`ffd_geometry/` contains both actual pyGeo control-point arrays as CSV, Plot3D XYZ, and NPZ. `verification.json` records zero difference between the reconstructed pyGeo surfaces and the saved baseline/final surface coordinates. The lattice is an FFD enclosure, not the allowed surface box; the box constrains the airfoil surface rather than its control points. To regenerate these data, run `source/export_optimization_ffd.py --root RUN_DIRECTORY` in the MDO Lab environment, then `source/plot_ffd_mesh_comparison.py --root RUN_DIRECTORY` with NumPy and Matplotlib.

`solution.txt` prints optimizer-scaled variables; physical FFD displacements in meters and incidence in degrees are in each `design.json`. `source/replay_settings.json` supplies only the original fixed drag cap; active CFD controls are in the top-level `settings.txt`.

`optimization.hst` is the pyOptSparse binary history. `history.jsonl` records every completed CFD evaluation, including derivative checks and line-search trials. `history_summary.csv` is a compact table. `gradients.jsonl` records aerodynamic sensitivities and adjoint failure flags. `PSQP.out`, `solution.txt`, `optimizer_exit.json`, and `run_complete.json` document termination. `run.log` contains Python/MPI output. Native Fortran CFD output is also in `PSQP.out` because this run used Fortran output unit 6. Both logs are preserved.

`restarts/` preserves interrupted frozen-connectivity attempts. The continued run replays its completed optimizer evaluations and retains its drag cap. PSQP uses XMAX=0.05 in scaled coordinates, MET=1 (BFGS), and Powell curvature correction; this controls each numerical step and does not bound final incidence. The production CFD switches earlier to Newton–Krylov, uses a larger Krylov subspace and stronger preconditioning, disables adaptive linear tolerances, and penalizes failed CFD trials explicitly.

`full_connectivity_validation/` holds the final independent full-connectivity CFD run, restarted from the saved optimized flow. It rebuilds hole cutting and interpolation, without computing adjoints. Its settings and result files document the difference from the optimization model.

`source/` includes the input mesh, restart solution, FFD file, Python drivers, software versions/commits, and derivative-check results. The initial surface originates from `cases/naca0012_c050_h005_u200_am3`. Both baseline and optimized solutions use the original 88,548-cell mesh with global IDWarp deformation and frozen overset connectivity. An extended-overlap mesh was tested during setup but rejected because it did not cure full-update fragility. The baseline is solved at the current tighter tolerance before setting the drag cap.

## Reproduction

Use the recorded `mdolab/public:u24-gcc-ompi-stable` image digest. Copy the input `source/` into a fresh run directory mounted in the container. Do not overwrite this completed run when reproducing it. Inside the container:

```sh
source /home/mdolabuser/.bashrc_mdolab
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mpirun -np 4 python /project/run_slip_optimization.py --root /project/new_run --optimizer PSQP --restart-design /project/new_run/source/start_design.json
python /project/extract_optimization_fields.py --root /project/new_run
```

Run `plot_optimization_comparison.py --root /absolute/path/to/new_run` on a Python installation with NumPy and Matplotlib. Set `MPLCONFIGDIR` to a writable cache directory if needed. The mesh and restart paths are resolved from the supplied run root.

## Interpretation

These are preliminary results on one mesh, not a grid-independent design. ADflow uses fully turbulent SA RANS in this run. The pressure-loss plot retains small numerical total-pressure overshoots rather than clipping them away. The whole overset grid is warped together using USMesh and frozen interpolation connectivity/weights, as described in the installed ADflow documentation for this use case. Adjoint sensitivities use the same discrete model. Full recutting and fast weight updates were tested during setup and rejected for robustness problems. Frozen interpolation introduces a numerical approximation under deformation, in addition to the usual mesh/discretization errors. Geometry validity and solver failure flags are checked separately from the requested physical constraints. A successful local optimizer exit does not establish a global optimum.

`export_repair.txt` records a repaired post-optimization JSON serialization error. The exact final history variables were checked against the saved converged design; the repair used the exact last iterate and its already converged CFD solution. The executed pre-fix driver is retained, and the supplied driver uses `dataclasses.asdict` for the pyOptSparse exit record.
