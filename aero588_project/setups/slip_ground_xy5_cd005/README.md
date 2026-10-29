# Five-pair FFD setup with x/y motion and CD limit 0.05

Prepared October 29, 2026. Geometry and derivative checks passed. The subsequent [50-iteration optimization and independent CFD validation](../../optimizations/slip_ground_xy5_cd005_20261029/README.md) are saved separately. That run reached its iteration limit before stationarity, giving 29.28% more downforce in the optimization model and 29.05% in the independent check. The completed October 27 optimization retains its original parameterization and drag cap.

## Design variables and constraints

- Five chordwise stations, with one upper and one lower FFD point at each station in the section view.
- Two tied span planes: 5 × 2 × 2 = 20 physical control points, projecting into 10 visible points (five upper/lower pairs).
- Each projected point has independent global x and y displacement: 10 x variables + 10 y variables + one unbounded incidence variable = 21 design variables. Displacements are in meters; incidence is in degrees. Both span planes move together. No z design variables or explicit x/y displacement bounds are added.
- CD ≤ 0.05, normalized with the fixed initial reference chord 0.5 m and reference area 0.025 m². This prescribed cap takes precedence over the original baseline CD and any old replay settings.
- The actual airfoil surface stays inside −0.125 ≤ x ≤ 0.625 m and 0.05 ≤ y ≤ 0.35 m. The FFD lattice itself is not constrained to that box.
- Maximum section thickness ≥ 0.06 × initial chord = 0.03 m. This is a lower limit on the thickest section, not a local minimum everywhere, and it remains normalized by the initial 0.5 m chord even if x motion changes actual chord.
- Slip ground, 200 km/h, initial incidence −3°, and the original 88,548-cell mesh/CFD controls are retained. The default optimizer remains pyOptSparse PSQP.

## Thickness after horizontal movement

The previous paired-material-point calculation can overestimate or underestimate thickness when corresponding upper/lower points move to different chordwise locations. The new calculation intersects lines perpendicular to the current chord with the **actual piecewise-linear CFD surface**, including the trailing-edge closure. It measures continuous interior line segments and takes their maximum length. All projected surface vertices are considered, so the maximum is exact for this surface representation up to floating-point tolerance; the intersection is not limited to the five FFD stations.

The chord direction follows the deformed original leading-edge point and the midpoint of the two original trailing-edge points. A concave section can have multiple interior segments on a normal; exterior gaps are excluded. Self-intersecting/collapsed surface boundaries are rejected before CFD. Existing mesh/solver failure checks still apply; box and thickness alone do not ensure a valid CFD grid.

The maximum is piecewise differentiable. At a switch in the active intersection or maximum location its derivative can have a kink. The implemented branch derivative includes changes in the chord direction and is chained through pyGeo to both x and y variables. It does not require material-point correspondence between the two surfaces.

## Checks and files

- `source/geometry_setup.json`: all requested geometry/drag settings and variable counts.
- `source/geometry_derivative_check.json`: automatic pyGeo/constraint derivative check.
- `geometry_derivative_check_extended.json`: independent x, y, and incidence directions at a nonzero deformed design. Maximum absolute derivative error was below 2 × 10⁻⁹.
- `section_thickness_checks.json`: analytical rectangle and offset-vertex examples, rotation/translation invariance, concave-section intervals, invalid boundaries, and surface-coordinate derivatives.
- `ffd_geometry/`: exported baseline FFD coordinates; reconstruction matches the baseline surface exactly.
- `baseline/`: geometry only, not a newly computed CFD solution.
- `source/ground.cgns` and `source/restart_vol.cgns`: original baseline mesh and converged baseline flow for a future run.
- `source/`: saved scripts, baseline coordinates, software provenance, and SHA-256 manifest.

## Future run

Use a **fresh run directory** and the existing MDO Lab Docker environment. Copy this folder's `source/` into that run's `source/`. All Python helpers needed for this configuration are included. Inside the container, with that folder mounted as `/project`:

```sh
source /home/mdolabuser/.bashrc_mdolab
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mpirun -np 4 python /project/source/run_slip_optimization.py --root /project --ffd-stations 5 --shape-directions xy --drag-cap 0.05 --optimizer PSQP
```

The driver recomputes the CFD baseline and its derivative check before optimization. Do not hot-start from the October 27 optimizer history or pass its design JSON directly: that run used ten stations and a single y-only variable group. Even though both configurations have 21 total variables, they represent different design spaces. Any warm start requires an explicit geometry fit into the new lattice.
