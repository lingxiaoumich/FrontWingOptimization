# Moving no-slip ground validation case

NACA 0012, chord 0.50 m, minimum clearance 0.05 m, incidence -3 degrees, speed 200 km/h. This is separate from the earlier slip-ground run.

The patched image is `aero588/adflow-moving-wall:2.12.2`; see ../../solver_extensions/moving_wall/README.md. The ground is a viscous adiabatic wall at (55.5555555556,0,0) m/s. The airfoil remains stationary. Ground first-cell height is 100 micrometres; airfoil first spacing is 2 micrometres.

The converged boundary-validation run used four MPI ranks and took 189.1 seconds. Residual reduction: 9.9000e-09. Airfoil-only CL: -1.14820011; CD: 0.02360433. Lift/drag iteration monitors include the road; the result JSON uses the separately registered airfoil-only functions.

Maximum y-plus: airfoil 0.470, ground 5.454. The boundary condition is validated, but the road grid needs a refinement/convergence study before optimization. The supplemental refinement_check_15um directory is explicitly **unconverged**, despite improved wall resolution; its forces must not be treated as final CFD results.

No optimization or grid-independence claim is made. Original cases and the proposal are unchanged.
