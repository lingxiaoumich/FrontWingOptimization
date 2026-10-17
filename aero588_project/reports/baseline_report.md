# NACA 0012 ground-effect baseline

Status: CONVERGED. Zero-incidence baseline, not an optimized design.

- Chord: 0.70 m; speed: 100 km/h; lowest-point clearance: 0.20 m.
- Ground: ideal slip/symmetry plane. Road shear is omitted; this is not a moving no-slip road simulation.
- Air: 288.15 K and 101325 Pa; Reynolds number 1331194; Mach 0.08163.
- Fully turbulent Spalart-Allmaras RANS. One spanwise cell with symmetry boundaries models a two-dimensional section.
- Standard NACA 0012 finite trailing edge, 1.764 mm thick. The earlier sharp trailing-edge coordinate file is not used for this CFD mesh.
- 77093 cells before overset blanking; no nonpositive corner Jacobians or orphan overset cells found.
- Maximum wall y+: 0.512; median wall y+: 0.223.

## Computed baseline quantities

Use these only within the stated model and convergence status.

| Quantity | Value |
|---|---:|
| Lift coefficient, upward positive | -0.092717 |
| Drag coefficient | 0.011744 |
| Downforce per unit span | 30.673 N/m |
| Drag per unit span | 3.885 N/m |
| Downforce / drag | 7.895 |
| Residual / original initial residual | 4.487e-09 |

## Limitations and next checks

Mesh independence is not established. The ground has no viscous road boundary layer. Transition, finite-span effects, wheels, and the vehicle body are excluded. These are preliminary isolated-airfoil results. Verify grid sensitivity and ground treatment before using them for racing-car design decisions or optimization claims.

Angle-of-attack bounds are -4 to +4 degrees. For a ground-effect angle study, rotate the airfoil geometry while keeping the freestream horizontal and enforcing the 0.20 m minimum clearance; changing freestream angle alone would create a physically inconsistent road setup.

## Reproduction

Use mdolab/public:u24-gcc-ompi-stable with /home/mdolabuser/.bashrc_mdolab sourced. Mount a local working project at /project. Run scripts/build_ground_mesh.py, then run scripts/run_baseline.py with four MPI ranks. Add --restart to continue from runs/baseline_a0/baseline_a0_000_vol.cgns. Run scripts/plot_mesh.py and scripts/postprocess_baseline.py for the plots. Copy deliverables to G:\My Drive\aero588_project.
