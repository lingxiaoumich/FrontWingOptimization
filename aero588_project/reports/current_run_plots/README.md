# Final NACA 0012 baseline plots

These figures use the saved final CFD volume solution and surface coefficients,
not a new CFD run. The model is fully turbulent SA RANS at zero incidence with
an ideal slip/symmetry ground. Mesh independence is not established.

## Figures

- `mesh.png` / `mesh.svg`: raw overlapping structured grids, with selected lines
  in the overview and every grid line in the leading-edge detail. There are
  77,093 cells before blanking. The near-airfoil grid is blue; the background
  grid is grey.
- `total_pressure.png` / `total_pressure.svg`: absolute local total pressure and
  a wake close-up colored by `(p0_inf - p0) / q_inf`.
- `surface_cp.png` / `surface_cp.svg`: unsmoothed upper/lower surface Cp and a
  leading-edge close-up. The finite trailing-edge face is excluded from these
  two curves. The Cp axis is inverted.
- `static_pressure_cp.png` / `static_pressure_cp.svg`: an additional static-Cp
  contour showing pressure changes around the airfoil and in the ground gap.

## Data and calculation

Volume source: `runs/baseline_a0/baseline_a0_000_vol.cgns`.
Surface source: `results/baseline_surface.csv`.
Source hashes are saved in `source_provenance.json`.

The CGNS file is ADF format. It was read with the official CGNS v4.4.0 ADF C
implementation in a temporary library. The normalized pressure was converted
using the CGNS reference pressure (101325.9482117028 Pa), rather than the rounded
101325 Pa input in the project report. CGNS reference Mach is 0.08162803501206718
and gamma is 1.4.

Local total pressure is reconstructed from static pressure and local Mach:

`p0 = p * (1 + (gamma - 1) * M**2 / 2)**(gamma / (gamma - 1))`

This is the ideal-gas local isentropic stagnation relation; it does not assume
that the full viscous flow field is isentropic. See
[NASA's isentropic flow equations](https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html).
The freestream total pressure is 101799.34000092941 Pa and freestream dynamic
pressure is 472.6040083724139 Pa.

Rind/ghost cells are removed. Compute cells (`Iblank=1`) and interpolation cells
(`Iblank=-1`) are plotted. Blanked, flooded, and flood-seed cells are excluded.
Contour interpolation uses triangles between adjacent cell centers within each
structured block; triangles involving excluded cells are masked. The airfoil
interior is covered with its saved surface outline. There is no interpolation
across the airfoil or a global triangulation across overset holes.

The wake-loss scale is limited to -0.05 through 0.50 to show wake detail; values
outside that range use the end colors. The underlying field is unchanged.
Near-wall loss coefficients reach 1.537. Small numerical total-pressure
overshoots above the freestream value (up to about 21 Pa in the near grid) are
retained, including negative loss coefficients; they have not been clipped to
zero. The absolute-pressure contour covers the full active-field range.

Pressure normalization was checked against the saved CoefPressure field; maximum
Cp reconstruction disagreement is about 1.28e-5, consistent with single-precision
saved fields. Mach reconstructed from saved pressure, density, and velocity
agrees with the saved Mach field within 1.08e-8. Details are in
`field_validation.json`.

## Reproduce the figures

From the project directory, run `python3 scripts/plot_current_run.py` with NumPy
and Matplotlib installed. The script uses the portable cached arrays in
`baseline_fields.npz`, so the CGNS library is not required for replotting.
