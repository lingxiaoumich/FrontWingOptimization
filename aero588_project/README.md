# Racing-car airfoil optimization

Starting geometry: NACA 0012. Scope: two-dimensional airfoil only; no finite-span wing. Chord: 0.70 m. Vehicle speed: 100 km/h (27.7778 m/s). Ground clearance: 0.20 m measured from the lowest point of the airfoil to the ground. Clearance-to-chord ratio: 0.285714. Optimization objective: maximize positive downforce divided by drag; require positive downforce and positive drag.
Framework: MACH-Aero, with ADflow for CFD and adjoint derivatives, pyGeo for geometry parameterization, IDWarp for mesh deformation, and pyOptSparse for optimization.

## Project status
- **Solver, mesh and derivative studies** of the design case (0.70 m, 0.20 m, 100 km/h, α = 0 and −3°) with the repo's BCs: ANK on/off test, a three-level mesh study (77k–307k cells) and an adjoint vs finite-difference check of the optimization variables. See [study/README.md](study/README.md). The shape parameterization and optimization plan are in [docs/shape_control_and_optimization.md](docs/shape_control_and_optimization.md), with an interactive [FFD explorer](docs/ffd_explorer.html). The MACH-Aero container for this folder is `env/fwo.sh`.
- An additional CFD case has converged at 0.50 m chord, 0.05 m minimum clearance, 200 km/h, and −3° nose-down incidence with horizontal freestream. See [case report and plots](cases/naca0012_c050_h005_u200_am3/README.md). It uses the same SA RANS/slip-ground model, with a new 88,548-cell mesh; mesh independence is not established.
- Project folder is accessible; it was empty at initial inspection.
- Ubuntu 24.04 starts successfully after the earlier transient WSL service crash. Docker and mdolab/public:u24-gcc-ompi-stable are installed. Package imports, SLSQP initialization, and two-rank MPI checks passed; see reports/environment_check.json. Container image digest: sha256:f29af49b0877de4ca040c84cc9b838cfca3690a9b4bcdf21d61bef118a9535af.
- A zero-incidence SA RANS baseline has converged. See reports/baseline_report.md and results/baseline_a0.json. No optimization has been run; mesh independence is not established. The current ground is an ideal slip plane, not a moving no-slip road.

## Required design inputs
- Confirmed geometry: single NACA 0012 section with 0.70 m chord and 0.20 m clearance to its lowest point.
- Atmospheric conditions; design speed is confirmed as 100 km/h.
- Confirmed objective: maximize downforce-to-drag ratio. Minimum useful downforce, if required, remains to be specified. Compute airfoil forces only, excluding the ground.
- Incidence bounds: -4 to +4 degrees. Airfoil section-shape changes and their bounds remain to be confirmed.
- Thickness, packaging, and assignment or competition constraints.
- Available Linux/cluster environment and compute budget.

## Planned workflow
1. Confirm requirements and a working MACH-Aero environment.
2. Generate baseline NACA 0012 coordinates and a quasi-two-dimensional CFD setup.
3. A ground-effect overset mesh with an ideal slip ground is built. Verify moving no-slip road support and adjoint compatibility before extending this model.
4. Establish a converged baseline and check mesh sensitivity.
5. Parameterize geometry, add constraints, and verify selected derivatives.
6. Run optimization and re-evaluate the final design.
7. Save geometry, meshes, scripts, histories, plots, and a comparison report here.

NACA 0012 is symmetric; downforce requires an appropriate incidence or optimized camber. This isolated airfoil study includes ground effect. It excludes finite-span effects, endplates, the car body, and wheels. Maintain at least 0.20 m clearance from the lowest surface point to the ground during geometry optimization. Whether the airfoil may move vertically to maintain exactly 0.20 m clearance is not yet specified.

## Documentation
https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/
https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machAeroTutorials/airfoilanalysis_overview.html
https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machAeroTutorials/opt_overview.html





## Storage
The solver uses a local working directory because direct Google Drive mounts are unreliable inside Docker. Save all reports and aerodynamic deliverables to this project folder. No Ubuntu reinstall, distribution deletion, or Windows reboot was required.


