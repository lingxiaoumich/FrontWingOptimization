# Five-pair x/y slip-ground optimization

The PSQP run returned: **iteration limit reached; stationarity not reached**. This is a preliminary single-grid result. The selected design is `final/`, with 29.285% more downforce than the baseline. An independent CFD run with full overset connectivity gives 29.045% more downforce and reports feasibility = True.

The model uses five upper/lower FFD pairs, tied across two span planes: 10 independent x displacements, 10 y displacements, and free incidence (21 variables). CD ≤ 0.05. Initial chord is 0.5 m, freestream is 200 km/h, and the ground is an ideal slip wall. The actual airfoil surface is constrained to x = [−0.125, 0.625] m and y = [0.05, 0.35] m. Maximum section thickness must be at least 0.03 m, calculated from actual-surface chord-normal intersections.

## Results and provenance

- `comparison_results.txt` and `.json`: baseline/selected force coefficients, dimensional loads, clearance, chord, thickness, incidence, solver status, and independent validation.
- `settings.txt` and `.json`: geometry, CFD, mesh, optimizer settings, units, and resolved solver defaults.
- `plots/`: PNG/SVG comparisons of geometry, mesh, Cp, pressure contours/loss, loads, iteration history, and FFD movement in x and y.
- `plots/control_point_history`, `plots/cl_cd_history`, `plots/history_separation_check.json`: per-iteration control-point paths with the surface box, C_L c / C_D c histories and drag polar, and the suction-side reversed-flow check of the baseline and final flows. Made by `scripts/plot_optimization_trajectory.py --root <this folder>` (added on branch `opt1031`).
- `baseline/`, `final/`, `best_feasible/`: surface coordinates, physical design variables, CGNS flow/mesh solutions, extracted fields, and result JSON.
- `ffd_geometry/`: actual pyGeo lattice coordinates and reconstruction verification.
- `validation.json`: numerical/geometry/mesh checks. `full_connectivity_validation/` contains the independent final CFD evaluation, log, settings, and mesh checks.
- `run.log`, `PSQP.out`, `history.jsonl`, `gradients.jsonl`, `optimization.hst`, and `optimizer_exit.json`: full run and optimization history.
- `source/`: executed scripts, input mesh, initial flow/geometry, software versions, derivative checks, and source checksums.

Optimization uses frozen overset interpolation connectivity while IDWarp deforms the whole mesh. The independent final check rebuilds connectivity. Grid independence has not been established. FFD and incidence remain unbounded explicitly; PSQP's step limit is numerical, and the actual surface box/thickness are the physical geometric constraints. The thickness maximum is piecewise differentiable when its active location changes.

The original October 27 optimization and the October 29 prepared setup remain in their separate folders.

The optimization itself took 44.27 minutes; baseline setup/checks and the final re-solve bring the driver wall time to 45.96 minutes. The independent full-connectivity evaluation took 14.21 seconds. `source/execution.json` records the exact image digest, command, and runtime environment. `final_checks.json` records zero overset orphans, positive saved grid-corner orientations, and verified source checksums. All 57 recorded CFD evaluations converged, and all 52 optimizer/final evaluations were feasible. One preoptimization finite-difference perturbation slightly crosses the design clearance constraint; it is not an accepted optimization iterate (see `final_checks.json`).
