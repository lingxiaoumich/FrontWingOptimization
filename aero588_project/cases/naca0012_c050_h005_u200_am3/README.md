# NACA 0012 ground-effect CFD: 0.50 m chord, 5 cm clearance, 200 km/h, −3° nose down

Status: **CONVERGED**. This is a CFD analysis, not an optimized design.

## Conditions and geometry

| Parameter | Value |
|---|---:|
| Chord | 0.50 m |
| Minimum surface-to-ground clearance after rotation | 0.05 m |
| Clearance/chord | 0.10 |
| Freestream speed | 200 km/h (55.5556 m/s) |
| Aerodynamic incidence | −3° nose down |
| LE-to-TE chord-line angle relative to +x | +3° |
| Freestream angle | 0°, horizontal |
| Temperature | 288.15 K |
| Nominal atmospheric pressure | 101325 Pa |
| Mach number | 0.163256 |
| Reynolds number based on chord | 1901705 |
| Extruded span / spanwise cells | 0.05 m / 1 |

Negative aerodynamic incidence was implemented by tilting the leading edge below
its trailing edge. The airfoil was then shifted vertically until its lowest
surface point was exactly 0.05 m above the ground. The freestream was kept
horizontal. The standard finite NACA trailing-edge thickness is 1.26 mm.

## CFD and mesh

- ADflow 2.12.2 in the installed MDOLab arm64 Docker image, using four MPI ranks.
  The earlier Windows/WSL baseline used ADflow 2.13.1, so this is not a comparison
  with identical solver builds. Exact local versions are in reports/environment.json.
- Steady compressible RANS; fully turbulent Spalart–Allmaras; wall functions off.
- Second-order central finite volume mean-flow scheme with JST scalar artificial
  dissipation; first-order turbulence-equation discretization.
- ANK followed by NK; L2 convergence settings 1e-8; single-grid cycle.
- Ideal slip/symmetry ground, viscous no-slip airfoil, spanwise symmetry, farfield
  upstream/downstream/top boundaries.
- Near-field overset grid: 23,040 cells. Background grid: 65,508 cells.
  Total before blanking: **88,548 cells**.
- The ground-gap background spacing is 2.5 mm. First near-wall spacing is 2 µm.
- Near-field march distance was reduced to 0.025 m to avoid crossing the ground.
  Lowest near-field mesh point is 0.025314 m.
- No nonpositive checked corner Jacobians and no overset orphan cells.
- Wall y+: maximum **0.433**; median 0.168.

## Results

Forces are for the airfoil only. Lift is positive upwards. Dimensional forces
below are per metre of span, not the force on a full finite-span front wing.

| Quantity | Value |
|---|---:|
| CL | -1.151797 |
| CD | 0.022961 |
| Downforce / drag | 50.162 |
| Downforce | 1088.678 N/m |
| Drag | 21.703 N/m |
| Final / initial residual | 6.853e-09 |
| CFD including setup and output | 150.42 s |
| Nonlinear iterations | 148 |
| Solver work iterations | 5722 |

The final two logged force samples differ by 1.724e-04
in CL and 5.478e-07 in CD. Solver failure flags are false.
The solution started from a fresh initial state; no restart was used.

## Figures

- [Mesh](reports/plots/mesh.png), including a leading-edge detail.
- [Total pressure and wake loss](reports/plots/total_pressure.png).
- [Surface Cp](reports/plots/surface_cp.png), including a leading-edge detail.
- [Static-pressure Cp contour](reports/plots/static_pressure_cp.png).

PNG and SVG versions are saved in reports/plots. Surface Cp is plotted using
chord-aligned coordinates, not world-x divided by chord; the finite trailing-edge
face is excluded from the upper/lower curves. Data are not smoothed.

Total pressure is calculated from the saved normalized static pressure and Mach
using p0 = p [1 + (gamma−1) M²/2]^(gamma/(gamma−1)), with gamma = 1.4.
The exact CGNS pressure reference is 101325.948212 Pa,
freestream total pressure is 103228.993939 Pa, and
freestream dynamic pressure is 1890.416033 Pa. Pressure normalization
is checked against saved Cp and Mach is checked against pressure/density/velocity.
Rind cells and blanked/flooded cells are excluded. The near-grid periodic seam is
connected when contouring. The wake-loss color scale saturates above 1.10 and
below −0.08; the raw arrays are unchanged. Numerical total-pressure overshoots
are retained. Source hashes and field checks are saved alongside the plots.

## Reproduction

The exact run script is run_case.py and its inputs are case_config.json. Mount
this case directory at /project in mdolab/public:u24-gcc-ompi-stable and source
/home/mdolabuser/.bashrc_mdolab. Run:

```sh
mpirun -np 4 python /project/run_case.py build --case-root /project
mpirun -np 4 python /project/run_case.py solve --case-root /project
```

The shared project script scripts/plot_current_run.py accepts --case-directory
pointing here and uses the extracted reports/plots/baseline_fields.npz cache.
The shared extraction script scripts/extract_cfd_fields.py accepts this case
and a CGNS 4.x ADF shared library built with 64-bit cgsize_t.

Mesh independence is not established. This model omits viscous road shear,
transition, finite-span effects, endplates, wheels, and the car body.
