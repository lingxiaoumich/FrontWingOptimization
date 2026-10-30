FRONT WING OPTIMIZATION - MAIN FILE GUIDE
Repository: https://github.com/lingxiaoumich/FrontWingOptimization

HOW TO READ THIS GUIDE
Paths below are relative to the GitHub repository root. In the separate local
aero588_project working folder, remove the leading "aero588_project/" prefix.
The repository-root and project-root main_readme.txt files contain this same text.
The inventory covers each version-controlled file, including archived copies.
Git internals and new ignored caches are excluded. Archived source files should
stay with their original run; scripts/ contains the current working tools.

START HERE
1. aero588_project/README.md: current model, constraints, results, and reproduction links.
2. aero588_project/scripts/run_slip_optimization.py: main optimization driver.
3. aero588_project/scripts/airfoil_geometry.py: actual-surface thickness and topology.
4. aero588_project/setups/slip_ground_xy5_cd005/README.md: latest input setup and launch instructions.
5. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/comparison_results.txt: latest before/after results.
6. That run's settings.txt: CFD/optimizer configuration; plots/: PNG/SVG figures.
7. That run's final/design.json and final/surface_geometry.csv: physical final design.

CURRENT MODEL AND RESULTS
NACA 0012 initial section; initial chord 0.5 m; freestream 200 km/h; slip ground.
Five upper/lower FFD pairs are tied across two span planes: 20 physical points,
10 independent x variables, 10 y variables, and one free incidence variable.
Objective: minimize upward-positive CL (maximize downforce). CD <= 0.05.
Actual surface box: x = [-0.125, 0.625] m; y = [0.05, 0.35] m.
Maximum section thickness >= 0.03 m, measured by actual-surface chord-normal
intersections. This is not a 3 cm local minimum everywhere. Reference dimensions
remain based on the initial chord even if the actual chord changes.

October 29: CL -1.15181 -> -1.48912; CD 0.0229614 -> 0.0269161; downforce
+29.28%. Independent full-connectivity CFD check: +29.05% downforce.
October 27: older y-only, ten-pair parameterization and baseline-CD cap;
downforce +15.26%, independent check +15.13%. Different caps/parameterizations
mean these two runs do not isolate the benefit of x motion.
Both optimizations reached their 50-iteration limit before stationarity.
They are feasible single-grid improvements, not converged or grid-independent optima.

FOLDER MAP
scripts/           Current generators, CFD/optimizer drivers, checks, and plotting tools.
setups/            Prepared geometry/input snapshots and their pre-run checks.
optimizations/     Full run settings, histories, sources, mesh/flow files, plots, and checks.
cases/             Baseline slip-ground and separate moving-ground CFD cases.
solver_extensions/ Optional moving no-slip ground patch, Docker build, and validation.
geometry/, meshes/, results/, runs/  Original 0.70 m / 100 km/h legacy baseline.
reports/           Baseline summaries, environment checks, and proposal figure options.
output/pdf/ and proposal latex/      Proposal sources, graphics, and PDF deliverables.
tmp/pdfs/          Rendered proposal/presentation images retained from document preparation.

RUN/RESULT CONVENTIONS
baseline/        Initial CFD solution, or geometry-only data when under setups/.
final/           Optimizer-returned design; authoritative status is in run_complete.json.
best_feasible/   Lowest-CL recorded feasible candidate; can differ from final/.
source/          Input/driver snapshot and software/hash provenance for that exact run.
full_connectivity_validation/  Independent check of the selected final geometry.
restarts/ and feasible_recovery/  Earlier attempts/recovery evidence, not final results.

CGNS files retain raw meshes/flow fields; NPZ caches and CSV surface data feed
the plots. solution.txt may print scaled variables; design.json stores physical
metres/degrees. history.jsonl includes baseline and finite-difference checks as
well as optimizer trials, so CFD evaluation counts differ from major iterations.
A derivative-check perturbation can cross a design bound without being an
accepted optimization iterate. Read result/status JSON before using force data.

REPRODUCTION NOTES
CFD/adjoints and pyGeo/IDWarp require the recorded MDO Lab Linux/Docker environment;
they are not standalone macOS Python scripts. Most main optimization helpers take
--root RUN_DIRECTORY; the single-case driver takes --case-root CASE_DIRECTORY.
Use a fresh writable run directory with source/ inputs copied from the chosen
setup/run. Never overwrite a completed run or hot-start between incompatible FFD
parameterizations just because they happen to have the same total variable count.
The recorded source/execution.json and software_versions.json pin the actual run.
Plotting uses NumPy/Matplotlib; CGNS extraction uses the matching container library.
The moving-wall extension uses a separate patched image. The published optimization
runs use slip ground, not that extension. The 15 micrometre moving-road refinement
and stiff 2 micrometre attempts are explicitly not accepted converged results.

To refresh this guide after repository changes:
  python3 aero588_project/scripts/generate_main_readme.py --repo-root .
Stage any new project files first so git ls-files includes them.

COMPLETE FILE-BY-FILE INVENTORY

Total version-controlled file entries: 635

0001. .gitignore
    Ignore new operating-system caches, Python caches, and local credential files; CFD results/logs remain version controlled.

0002. README.md
    GitHub landing page: current model, two run summaries, status limitations, and navigation links.

0003. aero588_project/.DS_Store
    Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.

0004. aero588_project/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0005. aero588_project/case_config.json
    Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline.

0006. aero588_project/cases/naca0012_c050_h005_u200_am3/.DS_Store
    Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.

0007. aero588_project/cases/naca0012_c050_h005_u200_am3/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0008. aero588_project/cases/naca0012_c050_h005_u200_am3/case_config.json
    Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline.

0009. aero588_project/cases/naca0012_c050_h005_u200_am3/geometry/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.

0010. aero588_project/cases/naca0012_c050_h005_u200_am3/geometry/airfoil_surface.xyz
    Extruded physical airfoil surface coordinates for mesh generation; not an FFD control-point lattice.

0011. aero588_project/cases/naca0012_c050_h005_u200_am3/meshes/airfoil_near.cgns
    Body-fitted near-airfoil mesh block generated by pyHyp.

0012. aero588_project/cases/naca0012_c050_h005_u200_am3/meshes/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0013. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/.DS_Store
    Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.

0014. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/cfd.log
    Recorded execution output for the case CFD solve. Inspect the associated result/status before treating a run as converged.

0015. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/environment.json
    Recorded solver/library environment for the enclosing CFD case.

0016. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/mesh_build.log
    Recorded execution output for mesh generation. Inspect the associated result/status before treating a run as converged.

0017. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0018. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/naca0012_c050_h005_u200_am3_plots.zip
    Bundled CFD comparison figures for sharing.

0019. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/baseline_fields.npz
    Cached NumPy arrays extracted from a baseline CGNS solution for plotting, with ghost/rind data removed.

0020. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/field_validation.json
    Validation information for extracted plotting fields and their source solution.

0021. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/mesh.png
    Selected mesh lines around the airfoil and ground. PNG raster export.

0022. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/mesh.svg
    Selected mesh lines around the airfoil and ground. SVG vector export.

0023. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/source_provenance.json
    Origin paths/hashes and reference information for the plotting cache and source solution.

0024. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/static_pressure_cp.png
    Static-pressure coefficient contour of the saved baseline/single CFD case. PNG raster export.

0025. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/static_pressure_cp.svg
    Static-pressure coefficient contour of the saved baseline/single CFD case. SVG vector export.

0026. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/surface_cp.png
    Surface Cp diagram of the saved baseline/single CFD case. PNG raster export.

0027. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/surface_cp.svg
    Surface Cp diagram of the saved baseline/single CFD case. SVG vector export.

0028. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/total_pressure.png
    Total-pressure contour of the saved baseline/single CFD case. PNG raster export.

0029. aero588_project/cases/naca0012_c050_h005_u200_am3/reports/plots/total_pressure.svg
    Total-pressure contour of the saved baseline/single CFD case. SVG vector export.

0030. aero588_project/cases/naca0012_c050_h005_u200_am3/results/cfd_result.json
    Case CFD results and convergence flags; use the enclosing README to distinguish accepted and explicitly unconverged runs.

0031. aero588_project/cases/naca0012_c050_h005_u200_am3/results/surface.csv
    Case surface positions, normalized coordinates, Cp, skin-friction magnitude, and y-plus extracted from the CFD solution.

0032. aero588_project/cases/naca0012_c050_h005_u200_am3/run_case.py
    Case-specific build/solve driver saved with this case; use that case configuration and its recorded solver environment.

0033. aero588_project/cases/naca0012_c050_h005_u200_am3/runs/naca0012_c050_h005_u200_am3_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0034. aero588_project/cases/naca0012_c050_h005_u200_am3/runs/naca0012_c050_h005_u200_am3_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0035. aero588_project/cases/naca0012_c050_h005_u200_am3/runs/naca0012_c050_h005_u200_am3_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0036. aero588_project/cases/naca0012_c050_h005_u200_am3/runs/naca0012_c050_h005_u200_am3_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0037. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0038. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/case_config.json
    Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline.

0039. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/geometry/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.

0040. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/geometry/airfoil_surface.xyz
    Extruded physical airfoil surface coordinates for mesh generation; not an FFD control-point lattice.

0041. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/meshes/airfoil_near.cgns
    Body-fitted near-airfoil mesh block generated by pyHyp.

0042. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/meshes/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0043. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0044. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/case_config.json
    Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0045. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/geometry/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0046. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/geometry/airfoil_surface.xyz
    Extruded physical airfoil surface coordinates for mesh generation; not an FFD control-point lattice. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0047. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/meshes/airfoil_near.cgns
    Body-fitted near-airfoil mesh block generated by pyHyp. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0048. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/meshes/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0049. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/cfd.log
    Recorded execution output for the case CFD solve. Inspect the associated result/status before treating a run as converged. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0050. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/incomplete_restart.log
    Recorded execution output for the incomplete fine-road-grid restart. Inspect the associated result/status before treating a run as converged. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0051. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/mesh_build.log
    Recorded execution output for mesh generation. Inspect the associated result/status before treating a run as converged. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0052. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0053. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/pre_restart_result.json
    Preserved provisional result from before the fine-road-grid restart; not an accepted final aerodynamic result. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0054. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/startup_cfd.log
    Recorded execution output for the fine-road-grid startup attempt. Inspect the associated result/status before treating a run as converged. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0055. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/reports/wall_resolution.json
    Airfoil/road wall-spacing and y-plus assessment; a small y-plus alone does not establish CFD convergence. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0056. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/results/cfd_result.json
    Case CFD results and convergence flags; use the enclosing README to distinguish accepted and explicitly unconverged runs. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0057. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/run_case.py
    Case-specific build/solve driver saved with this case; use that case configuration and its recorded solver environment. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0058. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/runs/naca0012_c050_h005_u200_am3_moving_ground_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0059. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/runs/naca0012_c050_h005_u200_am3_moving_ground_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0060. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/runs/naca0012_c050_h005_u200_am3_moving_ground_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0061. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/refinement_check_15um/runs/naca0012_c050_h005_u200_am3_moving_ground_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.

0062. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/reports/cfd.log
    Recorded execution output for the case CFD solve. Inspect the associated result/status before treating a run as converged.

0063. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/reports/mesh_build.log
    Recorded execution output for mesh generation. Inspect the associated result/status before treating a run as converged.

0064. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/reports/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0065. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/reports/wall_resolution.json
    Airfoil/road wall-spacing and y-plus assessment; a small y-plus alone does not establish CFD convergence.

0066. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/results/cfd_result.json
    Case CFD results and convergence flags; use the enclosing README to distinguish accepted and explicitly unconverged runs.

0067. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/run_case.py
    Case-specific build/solve driver saved with this case; use that case configuration and its recorded solver environment.

0068. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/runs/naca0012_c050_h005_u200_am3_moving_ground_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0069. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/runs/naca0012_c050_h005_u200_am3_moving_ground_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0070. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/runs/naca0012_c050_h005_u200_am3_moving_ground_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0071. aero588_project/cases/naca0012_c050_h005_u200_am3_moving_ground/runs/naca0012_c050_h005_u200_am3_moving_ground_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0072. aero588_project/geometry/airfoil_surface.xyz
    Extruded physical airfoil surface coordinates for mesh generation; not an FFD control-point lattice.

0073. aero588_project/geometry/naca0012_baseline.dat
    Legacy NACA 0012 profile coordinates produced by generate_baseline.py.

0074. aero588_project/main_readme.txt
    This plain-text navigation guide and complete file-purpose inventory; identical copies are provided at repository and project roots.

0075. aero588_project/meshes/airfoil_near.cgns
    Body-fitted near-airfoil mesh block generated by pyHyp.

0076. aero588_project/meshes/ground_baseline.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0077. aero588_project/optimizations/slip_ground_box_20261027/PSQP.out
    Native optimizer iteration/stopping log; PSQP.out can also contain redirected Fortran CFD/adjoint solver output.

0078. aero588_project/optimizations/slip_ground_box_20261027/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0079. aero588_project/optimizations/slip_ground_box_20261027/RUN_STATUS.txt
    Human-readable lifecycle/stopping status; distinguish prepared geometry, completed CFD, and optimizer termination.

0080. aero588_project/optimizations/slip_ground_box_20261027/baseline/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0081. aero588_project/optimizations/slip_ground_box_20261027/baseline/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0082. aero588_project/optimizations/slip_ground_box_20261027/baseline/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0083. aero588_project/optimizations/slip_ground_box_20261027/baseline/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0084. aero588_project/optimizations/slip_ground_box_20261027/baseline/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0085. aero588_project/optimizations/slip_ground_box_20261027/baseline/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0086. aero588_project/optimizations/slip_ground_box_20261027/baseline/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0087. aero588_project/optimizations/slip_ground_box_20261027/baseline/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0088. aero588_project/optimizations/slip_ground_box_20261027/baseline/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0089. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0090. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0091. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0092. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0093. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0094. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0095. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0096. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0097. aero588_project/optimizations/slip_ground_box_20261027/best_feasible/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0098. aero588_project/optimizations/slip_ground_box_20261027/comparison_results.json
    Machine-readable before/after forces, geometry, and selected independent-validation results.

0099. aero588_project/optimizations/slip_ground_box_20261027/comparison_results.txt
    Human-readable before/after table, reference definitions, optimizer status, constraints, and interpretation limits.

0100. aero588_project/optimizations/slip_ground_box_20261027/export_repair.txt
    Explain the post-optimization JSON-export repair and how the recovered final design was verified.

0101. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0102. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0103. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0104. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0105. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0106. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0107. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0108. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0109. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0110. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/baseline/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0111. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0112. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0113. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0114. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0115. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0116. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0117. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0118. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0119. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/best_feasible/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0120. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/comparison_results.json
    Machine-readable before/after forces, geometry, and selected independent-validation results. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0121. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/comparison_results.txt
    Human-readable before/after table, reference definitions, optimizer status, constraints, and interpretation limits. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0122. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/connectivity_comparison.txt
    Compare the saved frozen-connectivity prediction with an independent rebuilt-connectivity CFD evaluation. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0123. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0124. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0125. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0126. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0127. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0128. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0129. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0130. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0131. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0132. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0133. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0134. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0135. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0136. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0137. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0138. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0139. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0140. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/baseline/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0141. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0142. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0143. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0144. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0145. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0146. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0147. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0148. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0149. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/best_feasible/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0150. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0151. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0152. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0153. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0154. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0155. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0156. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0157. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0158. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0159. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0160. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0161. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/latest.json
    Most recently completed CFD evaluation record; use final/result.json and run_complete.json for final conclusions. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0162. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0163. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/run_complete.json
    Driver completion record with returned optimizer status, final evaluation, best-feasible design, and wall times; file existence does not imply stationarity. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0164. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/seed_design.json
    Physical design used to initialize an earlier feasibility/recovery check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0165. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0166. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0167. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/aerodynamic_derivative_check.json
    Finite-difference versus adjoint aerodynamic derivative comparison and primal-failure flags. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0168. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0169. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0170. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0171. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0172. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/mesh_case_config.json
    Mesh-generation configuration preserved with the run input snapshot. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0173. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0174. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/replay_settings.json
    Original drag cap retained for the older continuation/recovery workflow; the current driver uses its explicit --drag-cap argument. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0175. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0176. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/full_connectivity_check/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0177. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0178. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/history_summary.csv
    Compact CSV summary of recorded CFD evaluation history; evaluations are not the same as optimizer major iterations. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0179. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0180. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/latest.json
    Most recently completed CFD evaluation record; use final/result.json and run_complete.json for final conclusions. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0181. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/cp_comparison.png
    Before/after surface pressure-coefficient curves in chord and global coordinates. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0182. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/cp_comparison.svg
    Before/after surface pressure-coefficient curves in chord and global coordinates. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0183. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/cp_contour_comparison.png
    Before/after static-pressure coefficient contours from saved CFD fields. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0184. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/cp_contour_comparison.svg
    Before/after static-pressure coefficient contours from saved CFD fields. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0185. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/geometry_comparison.png
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0186. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/geometry_comparison.svg
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0187. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/mesh_comparison.png
    Baseline and selected optimized meshes, plotted side by side. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0188. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/mesh_comparison.svg
    Baseline and selected optimized meshes, plotted side by side. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0189. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/optimization_history.png
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0190. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/optimization_history.svg
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0191. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/performance_comparison.png
    Before/after downforce, drag, and downforce/drag bars. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0192. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/performance_comparison.svg
    Before/after downforce, drag, and downforce/drag bars. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0193. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/pressure_loss_comparison.png
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0194. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/pressure_loss_comparison.svg
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0195. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/total_pressure_comparison.png
    Before/after total pressure normalized by the freestream stagnation pressure. PNG raster export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0196. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/plots/total_pressure_comparison.svg
    Before/after total pressure normalized by the freestream stagnation pressure. SVG vector export. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0197. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0198. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/run_complete.json
    Driver completion record with returned optimizer status, final evaluation, best-feasible design, and wall times; file existence does not imply stationarity. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0199. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/seed_design.json
    Physical design used to initialize an earlier feasibility/recovery check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0200. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0201. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0202. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/aerodynamic_derivative_check.json
    Finite-difference versus adjoint aerodynamic derivative comparison and primal-failure flags. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0203. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0204. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0205. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0206. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0207. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/mesh_case_config.json
    Mesh-generation configuration preserved with the run input snapshot. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0208. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0209. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/replay_settings.json
    Original drag cap retained for the older continuation/recovery workflow; the current driver uses its explicit --drag-cap argument. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0210. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0211. aero588_project/optimizations/slip_ground_box_20261027/feasible_recovery/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information. Earlier recovery design/check retained for provenance; not the final October 27 production design.

0212. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/baseline_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0213. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/baseline_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0214. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/control_points.npz
    Actual pyGeo baseline/final FFD coefficient arrays used by the FFD comparison plots.

0215. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/final_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0216. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/final_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0217. aero588_project/optimizations/slip_ground_box_20261027/ffd_geometry/verification.json
    FFD reconstruction errors and point counts; checks the exported lattice against the saved actual surface.

0218. aero588_project/optimizations/slip_ground_box_20261027/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0219. aero588_project/optimizations/slip_ground_box_20261027/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0220. aero588_project/optimizations/slip_ground_box_20261027/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0221. aero588_project/optimizations/slip_ground_box_20261027/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0222. aero588_project/optimizations/slip_ground_box_20261027/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0223. aero588_project/optimizations/slip_ground_box_20261027/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0224. aero588_project/optimizations/slip_ground_box_20261027/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0225. aero588_project/optimizations/slip_ground_box_20261027/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0226. aero588_project/optimizations/slip_ground_box_20261027/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0227. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0228. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0229. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0230. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0231. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0232. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0233. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0234. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0235. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0236. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0237. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0238. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0239. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0240. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0241. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0242. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/design_to_check.json
    Physical design variables selected as input to a geometry-derivative or independent-connectivity check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0243. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0244. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0245. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0246. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/mesh_case_config.json
    Mesh-generation configuration preserved with the run input snapshot. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0247. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0248. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/replay_settings.json
    Original drag cap retained for the older continuation/recovery workflow; the current driver uses its explicit --drag-cap argument. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0249. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0250. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0251. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0252. aero588_project/optimizations/slip_ground_box_20261027/full_connectivity_validation/validation_complete.json
    Completion/results record for the independent full-connectivity CFD evaluation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0253. aero588_project/optimizations/slip_ground_box_20261027/gradients.jsonl
    One JSON record per recorded adjoint sensitivity evaluation, with design-group gradients and failure flags.

0254. aero588_project/optimizations/slip_ground_box_20261027/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings.

0255. aero588_project/optimizations/slip_ground_box_20261027/history_summary.csv
    Compact CSV summary of recorded CFD evaluation history; evaluations are not the same as optimizer major iterations.

0256. aero588_project/optimizations/slip_ground_box_20261027/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result.

0257. aero588_project/optimizations/slip_ground_box_20261027/latest.json
    Most recently completed CFD evaluation record; use final/result.json and run_complete.json for final conclusions.

0258. aero588_project/optimizations/slip_ground_box_20261027/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem.

0259. aero588_project/optimizations/slip_ground_box_20261027/optimizer_exit.json
    Optimizer termination code/message, returned physical variables, timings, and evaluation count.

0260. aero588_project/optimizations/slip_ground_box_20261027/optimizer_iterations.csv
    PSQP major-iteration log data: objective, evaluation counts, constraint violation, and printed gradient norm.

0261. aero588_project/optimizations/slip_ground_box_20261027/plots/cp_comparison.png
    Before/after surface pressure-coefficient curves in chord and global coordinates. PNG raster export.

0262. aero588_project/optimizations/slip_ground_box_20261027/plots/cp_comparison.svg
    Before/after surface pressure-coefficient curves in chord and global coordinates. SVG vector export.

0263. aero588_project/optimizations/slip_ground_box_20261027/plots/cp_contour_comparison.png
    Before/after static-pressure coefficient contours from saved CFD fields. PNG raster export.

0264. aero588_project/optimizations/slip_ground_box_20261027/plots/cp_contour_comparison.svg
    Before/after static-pressure coefficient contours from saved CFD fields. SVG vector export.

0265. aero588_project/optimizations/slip_ground_box_20261027/plots/ffd_before_after_comparison.png
    Before/after airfoil and FFD overlays with control-point displacement plots. PNG raster export.

0266. aero588_project/optimizations/slip_ground_box_20261027/plots/ffd_before_after_comparison.svg
    Before/after airfoil and FFD overlays with control-point displacement plots. SVG vector export.

0267. aero588_project/optimizations/slip_ground_box_20261027/plots/ffd_mesh_comparison.png
    Actual FFD control points/lattice over the baseline and final computational meshes. PNG raster export.

0268. aero588_project/optimizations/slip_ground_box_20261027/plots/ffd_mesh_comparison.svg
    Actual FFD control points/lattice over the baseline and final computational meshes. SVG vector export.

0269. aero588_project/optimizations/slip_ground_box_20261027/plots/geometry_comparison.png
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. PNG raster export.

0270. aero588_project/optimizations/slip_ground_box_20261027/plots/geometry_comparison.svg
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. SVG vector export.

0271. aero588_project/optimizations/slip_ground_box_20261027/plots/mesh_comparison.png
    Baseline and selected optimized meshes, plotted side by side. PNG raster export.

0272. aero588_project/optimizations/slip_ground_box_20261027/plots/mesh_comparison.svg
    Baseline and selected optimized meshes, plotted side by side. SVG vector export.

0273. aero588_project/optimizations/slip_ground_box_20261027/plots/optimization_history.png
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. PNG raster export.

0274. aero588_project/optimizations/slip_ground_box_20261027/plots/optimization_history.svg
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. SVG vector export.

0275. aero588_project/optimizations/slip_ground_box_20261027/plots/optimizer_convergence.png
    PSQP objective, printed gradient norm, and constraint violation versus major iteration. PNG raster export.

0276. aero588_project/optimizations/slip_ground_box_20261027/plots/optimizer_convergence.svg
    PSQP objective, printed gradient norm, and constraint violation versus major iteration. SVG vector export.

0277. aero588_project/optimizations/slip_ground_box_20261027/plots/performance_comparison.png
    Before/after downforce, drag, and downforce/drag bars. PNG raster export.

0278. aero588_project/optimizations/slip_ground_box_20261027/plots/performance_comparison.svg
    Before/after downforce, drag, and downforce/drag bars. SVG vector export.

0279. aero588_project/optimizations/slip_ground_box_20261027/plots/pressure_loss_comparison.png
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. PNG raster export.

0280. aero588_project/optimizations/slip_ground_box_20261027/plots/pressure_loss_comparison.svg
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. SVG vector export.

0281. aero588_project/optimizations/slip_ground_box_20261027/plots/shape_change.png
    Close-up airfoil outline and changes in surface height. PNG raster export.

0282. aero588_project/optimizations/slip_ground_box_20261027/plots/shape_change.svg
    Close-up airfoil outline and changes in surface height. SVG vector export.

0283. aero588_project/optimizations/slip_ground_box_20261027/plots/total_pressure_comparison.png
    Before/after total pressure normalized by the freestream stagnation pressure. PNG raster export.

0284. aero588_project/optimizations/slip_ground_box_20261027/plots/total_pressure_comparison.svg
    Before/after total pressure normalized by the freestream stagnation pressure. SVG vector export.

0285. aero588_project/optimizations/slip_ground_box_20261027/plotting_summary.log
    Recorded execution output for figure generation. Inspect the associated result/status before treating a run as converged.

0286. aero588_project/optimizations/slip_ground_box_20261027/problem.txt
    Printed pyOptSparse objective, variables, scaling, and constraint definition for the run.

0287. aero588_project/optimizations/slip_ground_box_20261027/production_checks.json
    Audit counts for completed primal/adjoint solves and the feasibility/failure flags of the production optimization.

0288. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Archived setup/search attempt; not the selected final optimization result.

0289. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/gradients.jsonl
    One JSON record per recorded adjoint sensitivity evaluation, with design-group gradients and failure flags. Archived setup/search attempt; not the selected final optimization result.

0290. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Archived setup/search attempt; not the selected final optimization result.

0291. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem. Archived setup/search attempt; not the selected final optimization result.

0292. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Archived setup/search attempt; not the selected final optimization result.

0293. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_1/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Archived setup/search attempt; not the selected final optimization result.

0294. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_2/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Archived setup/search attempt; not the selected final optimization result.

0295. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_2/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Archived setup/search attempt; not the selected final optimization result.

0296. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_2/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem. Archived setup/search attempt; not the selected final optimization result.

0297. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_2/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Archived setup/search attempt; not the selected final optimization result.

0298. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_2/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Archived setup/search attempt; not the selected final optimization result.

0299. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Archived setup/search attempt; not the selected final optimization result.

0300. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Archived setup/search attempt; not the selected final optimization result.

0301. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result. Archived setup/search attempt; not the selected final optimization result.

0302. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem. Archived setup/search attempt; not the selected final optimization result.

0303. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Archived setup/search attempt; not the selected final optimization result.

0304. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_3/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Archived setup/search attempt; not the selected final optimization result.

0305. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package. Archived setup/search attempt; not the selected final optimization result.

0306. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/SLSQP.out
    Native optimizer iteration/stopping log; PSQP.out can also contain redirected Fortran CFD/adjoint solver output. Archived setup/search attempt; not the selected final optimization result.

0307. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/gradients.jsonl
    One JSON record per recorded adjoint sensitivity evaluation, with design-group gradients and failure flags. Archived setup/search attempt; not the selected final optimization result.

0308. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings. Archived setup/search attempt; not the selected final optimization result.

0309. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result. Archived setup/search attempt; not the selected final optimization result.

0310. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem. Archived setup/search attempt; not the selected final optimization result.

0311. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Archived setup/search attempt; not the selected final optimization result.

0312. aero588_project/optimizations/slip_ground_box_20261027/restarts/attempt_4_slsqp/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Archived setup/search attempt; not the selected final optimization result.

0313. aero588_project/optimizations/slip_ground_box_20261027/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged.

0314. aero588_project/optimizations/slip_ground_box_20261027/run_complete.json
    Driver completion record with returned optimizer status, final evaluation, best-feasible design, and wall times; file existence does not imply stationarity.

0315. aero588_project/optimizations/slip_ground_box_20261027/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check.

0316. aero588_project/optimizations/slip_ground_box_20261027/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded.

0317. aero588_project/optimizations/slip_ground_box_20261027/solution.txt
    Printed optimizer solution and stopping status; variables may be optimizer-scaled, so use design.json for physical values.

0318. aero588_project/optimizations/slip_ground_box_20261027/source/aerodynamic_derivative_check.json
    Finite-difference versus adjoint aerodynamic derivative comparison and primal-failure flags.

0319. aero588_project/optimizations/slip_ground_box_20261027/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.

0320. aero588_project/optimizations/slip_ground_box_20261027/source/check_export_serialization.py
    Archived regression check for serializing pyOptSparse termination data, NumPy arrays, and scalars to JSON. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0321. aero588_project/optimizations/slip_ground_box_20261027/source/check_optimization_geometry.py
    Check pyGeo derivatives of actual-surface box and thickness functions at a selected or perturbed design. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0322. aero588_project/optimizations/slip_ground_box_20261027/source/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0323. aero588_project/optimizations/slip_ground_box_20261027/source/executed_driver_before_export_fix.py
    Original executed optimization driver preserved before the JSON-export repair; provenance, not the current driver. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0324. aero588_project/optimizations/slip_ground_box_20261027/source/export_optimization_ffd.py
    Reconstruct saved pyGeo design variables, export actual FFD coefficients, and verify the reconstructed surfaces against saved geometry. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0325. aero588_project/optimizations/slip_ground_box_20261027/source/export_serialization_check.json
    Round-trip test result for the repaired optimizer-termination JSON serialization.

0326. aero588_project/optimizations/slip_ground_box_20261027/source/extract_optimization_fields.py
    Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0327. aero588_project/optimizations/slip_ground_box_20261027/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0328. aero588_project/optimizations/slip_ground_box_20261027/source/finalize_optimization_run.py
    Assemble resolved settings/provenance and check saved geometry, CFD convergence, mesh validity, and source checksums. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0329. aero588_project/optimizations/slip_ground_box_20261027/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation.

0330. aero588_project/optimizations/slip_ground_box_20261027/source/geometry_derivative_check_extended.json
    Additional independent geometry derivative directions and the design at which they were checked.

0331. aero588_project/optimizations/slip_ground_box_20261027/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0332. aero588_project/optimizations/slip_ground_box_20261027/source/mesh_case_config.json
    Mesh-generation configuration preserved with the run input snapshot.

0333. aero588_project/optimizations/slip_ground_box_20261027/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0334. aero588_project/optimizations/slip_ground_box_20261027/source/plot_ffd_mesh_comparison.py
    Plot saved baseline/final FFD lattices over meshes and compare actual x/y control-point motion; no CFD solve. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0335. aero588_project/optimizations/slip_ground_box_20261027/source/plot_optimization_comparison.py
    Generate before/after geometry, mesh, Cp, pressure/loss, load, and optimization-history figures plus numerical comparison reports. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0336. aero588_project/optimizations/slip_ground_box_20261027/source/plotting_versions.json
    Versions of host plotting libraries used for the saved figures.

0337. aero588_project/optimizations/slip_ground_box_20261027/source/repair_optimizer_export.py
    Archived recovery utility: verify the final history design and restore the result export after a JSON serialization failure. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0338. aero588_project/optimizations/slip_ground_box_20261027/source/replay_settings.json
    Original drag cap retained for the older continuation/recovery workflow; the current driver uses its explicit --drag-cap argument.

0339. aero588_project/optimizations/slip_ground_box_20261027/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check.

0340. aero588_project/optimizations/slip_ground_box_20261027/source/run_ground_case.py
    Build/solve a case from case_config.json; supports slip, stationary no-slip, and patched moving no-slip ground. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0341. aero588_project/optimizations/slip_ground_box_20261027/source/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0342. aero588_project/optimizations/slip_ground_box_20261027/source/setup_notes.txt
    Notes explaining the numerical setup or an archived attempted configuration.

0343. aero588_project/optimizations/slip_ground_box_20261027/source/sha256_manifest.json
    SHA-256 hashes of the associated input/source files; verifies the preserved run snapshot.

0344. aero588_project/optimizations/slip_ground_box_20261027/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information.

0345. aero588_project/optimizations/slip_ground_box_20261027/source/start_design.json
    Physical starting variables for an archived optimization continuation.

0346. aero588_project/optimizations/slip_ground_box_20261027/validation.json
    Saved solution/geometry/grid/convergence checks and recorded derivative-check information.

0347. aero588_project/optimizations/slip_ground_box_20261027/validation_summary.log
    Recorded execution output for the saved final-result checks. Inspect the associated result/status before treating a run as converged.

0348. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/PSQP.out
    Native optimizer iteration/stopping log; PSQP.out can also contain redirected Fortran CFD/adjoint solver output.

0349. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0350. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/RUN_STATUS.txt
    Human-readable lifecycle/stopping status; distinguish prepared geometry, completed CFD, and optimizer termination.

0351. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0352. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0353. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0354. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0355. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0356. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0357. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0358. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0359. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/baseline/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0360. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0361. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0362. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0363. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0364. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0365. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0366. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0367. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0368. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/best_feasible/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0369. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/comparison_results.json
    Machine-readable before/after forces, geometry, and selected independent-validation results.

0370. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/comparison_results.txt
    Human-readable before/after table, reference definitions, optimizer status, constraints, and interpretation limits.

0371. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/completion_summary.json
    Final completion-pipeline summary: selected result, optimizer exit, gains, independent validation, and finish time.

0372. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/baseline_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0373. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/baseline_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0374. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/best_feasible_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0375. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/best_feasible_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0376. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/control_points.npz
    Actual pyGeo baseline/final FFD coefficient arrays used by the FFD comparison plots.

0377. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/final_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0378. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/final_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0379. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/ffd_geometry/verification.json
    FFD reconstruction errors and point counts; checks the exported lattice against the saved actual surface.

0380. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0381. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.

0382. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0383. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0384. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0385. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0386. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.

0387. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.

0388. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0389. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/final_checks.json
    Final numerical/data audit, including saved-grid validity, convergence, source hashes, orphan count, and derivative-check exceptions.

0390. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0391. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/fields.npz
    Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0392. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/flow_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0393. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/flow_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0394. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/flow_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0395. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/flow_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0396. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/result.json
    Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0397. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/surface_cp.csv
    Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0398. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/final/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0399. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0400. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0401. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0402. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0403. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/airfoil_geometry.py
    Compute actual-surface maximum chord-normal thickness and its active-feature derivative; detect crossed/collapsed section boundaries. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0404. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0405. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/design_to_check.json
    Physical design variables selected as input to a geometry-derivative or independent-connectivity check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0406. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/extract_optimization_fields.py
    Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0407. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0408. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0409. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/geometry_setup.json
    FFD dimensions, variable groups/counts/units, drag cap, box bounds, and thickness definition for the setup. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0410. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0411. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/initial_design.json
    Zero-shape-change initial FFD variables and initial incidence for this parameterization. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0412. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0413. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0414. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0415. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0416. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/full_connectivity_validation/validation_complete.json
    Completion/results record for the independent full-connectivity CFD evaluation. Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.

0417. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/geometry_derivative_check_extended.json
    Additional independent geometry derivative directions and the design at which they were checked.

0418. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/gradients.jsonl
    One JSON record per recorded adjoint sensitivity evaluation, with design-group gradients and failure flags.

0419. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/history.jsonl
    One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings.

0420. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/history_summary.csv
    Compact CSV summary of recorded CFD evaluation history; evaluations are not the same as optimizer major iterations.

0421. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/incoming_design.json
    Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result.

0422. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/latest.json
    Most recently completed CFD evaluation record; use final/result.json and run_complete.json for final conclusions.

0423. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/optimization.hst
    Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem.

0424. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/optimizer_exit.json
    Optimizer termination code/message, returned physical variables, timings, and evaluation count.

0425. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/optimizer_iterations.csv
    PSQP major-iteration log data: objective, evaluation counts, constraint violation, and printed gradient norm.

0426. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/cp_comparison.png
    Before/after surface pressure-coefficient curves in chord and global coordinates. PNG raster export.

0427. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/cp_comparison.svg
    Before/after surface pressure-coefficient curves in chord and global coordinates. SVG vector export.

0428. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/cp_contour_comparison.png
    Before/after static-pressure coefficient contours from saved CFD fields. PNG raster export.

0429. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/cp_contour_comparison.svg
    Before/after static-pressure coefficient contours from saved CFD fields. SVG vector export.

0430. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/ffd_before_after_comparison.png
    Before/after airfoil and FFD overlays with control-point displacement plots. PNG raster export.

0431. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/ffd_before_after_comparison.svg
    Before/after airfoil and FFD overlays with control-point displacement plots. SVG vector export.

0432. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/ffd_mesh_comparison.png
    Actual FFD control points/lattice over the baseline and final computational meshes. PNG raster export.

0433. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/ffd_mesh_comparison.svg
    Actual FFD control points/lattice over the baseline and final computational meshes. SVG vector export.

0434. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/geometry_comparison.png
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. PNG raster export.

0435. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/geometry_comparison.svg
    Baseline and selected airfoil outlines relative to the allowed surface box and ground. SVG vector export.

0436. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/mesh_comparison.png
    Baseline and selected optimized meshes, plotted side by side. PNG raster export.

0437. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/mesh_comparison.svg
    Baseline and selected optimized meshes, plotted side by side. SVG vector export.

0438. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/optimization_history.png
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. PNG raster export.

0439. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/optimization_history.svg
    Loads, incidence, and CFD residual history over evaluations, including trials/checks. SVG vector export.

0440. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/optimizer_convergence.png
    PSQP objective, printed gradient norm, and constraint violation versus major iteration. PNG raster export.

0441. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/optimizer_convergence.svg
    PSQP objective, printed gradient norm, and constraint violation versus major iteration. SVG vector export.

0442. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/performance_comparison.png
    Before/after downforce, drag, and downforce/drag bars. PNG raster export.

0443. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/performance_comparison.svg
    Before/after downforce, drag, and downforce/drag bars. SVG vector export.

0444. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/pressure_loss_comparison.png
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. PNG raster export.

0445. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/pressure_loss_comparison.svg
    Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained. SVG vector export.

0446. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/shape_change.png
    Close-up airfoil outline and changes in surface height. PNG raster export.

0447. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/shape_change.svg
    Close-up airfoil outline and changes in surface height. SVG vector export.

0448. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/total_pressure_comparison.png
    Before/after total pressure normalized by the freestream stagnation pressure. PNG raster export.

0449. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/plots/total_pressure_comparison.svg
    Before/after total pressure normalized by the freestream stagnation pressure. SVG vector export.

0450. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/postprocess.log
    Recorded execution output for the extraction, independent CFD check, plotting, and completion pipeline. Inspect the associated result/status before treating a run as converged.

0451. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/problem.txt
    Printed pyOptSparse objective, variables, scaling, and constraint definition for the run.

0452. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/progress.json
    Mirrored optimization progress snapshot; retained for provenance rather than used as the authoritative final result.

0453. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/run.log
    Recorded execution output for the enclosing CFD/optimization driver. Inspect the associated result/status before treating a run as converged.

0454. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/run_complete.json
    Driver completion record with returned optimizer status, final evaluation, best-feasible design, and wall times; file existence does not imply stationarity.

0455. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/section_thickness_checks.json
    Analytical and derivative tests of the actual-surface chord-normal thickness calculation.

0456. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/settings.json
    Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check.

0457. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/settings.txt
    Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded.

0458. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/solution.txt
    Printed optimizer solution and stopping status; variables may be optimizer-scaled, so use design.json for physical values.

0459. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/aerodynamic_derivative_check.json
    Finite-difference versus adjoint aerodynamic derivative comparison and primal-failure flags.

0460. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.

0461. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/airfoil_geometry.py
    Compute actual-surface maximum chord-normal thickness and its active-feature derivative; detect crossed/collapsed section boundaries. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0462. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/check_optimization_geometry.py
    Check pyGeo derivatives of actual-surface box and thickness functions at a selected or perturbed design. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0463. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0464. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/check_section_thickness.py
    Test analytical thickness examples, shifted upper/lower vertices, concave sections, rotation invariance, invalid topology, and derivatives. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0465. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/execution.json
    Exact Docker image digest, launch command, thread settings, and staging/result locations used for the run.

0466. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/export_optimization_ffd.py
    Reconstruct saved pyGeo design variables, export actual FFD coefficients, and verify the reconstructed surfaces against saved geometry. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0467. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/extract_optimization_fields.py
    Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0468. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0469. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/finalize_optimization_run.py
    Assemble resolved settings/provenance and check saved geometry, CFD convergence, mesh validity, and source checksums. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0470. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/finish_xy_optimization.py
    Host completion worker for the October 29 run: wait for Docker, validate/extract final results, generate plots/reports, and mirror the results. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0471. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation.

0472. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/geometry_setup.json
    FFD dimensions, variable groups/counts/units, drag cap, box bounds, and thickness definition for the setup.

0473. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0474. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/initial_design.json
    Zero-shape-change initial FFD variables and initial incidence for this parameterization.

0475. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0476. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/plot_ffd_mesh_comparison.py
    Plot saved baseline/final FFD lattices over meshes and compare actual x/y control-point motion; no CFD solve. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0477. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/plot_optimization_comparison.py
    Generate before/after geometry, mesh, Cp, pressure/loss, load, and optimization-history figures plus numerical comparison reports. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0478. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check.

0479. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0480. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/sha256_manifest.json
    SHA-256 hashes of the associated input/source files; verifies the preserved run snapshot.

0481. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information.

0482. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/source/sync_optimization_progress.py
    Mirror stable files from local Docker staging into the project and write progress summaries; does not compute CFD. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0483. aero588_project/optimizations/slip_ground_xy5_cd005_20261029/validation.json
    Saved solution/geometry/grid/convergence checks and recorded derivative-check information.

0484. aero588_project/output/pdf/aer588_project_proposal.log
    LaTeX compiler diagnostics for the proposal PDF; not a CFD log.

0485. aero588_project/output/pdf/aer588_project_proposal.pdf
    Compiled/exported project-proposal PDF, not an aerodynamic result.

0486. aero588_project/output/pdf/aer588_project_proposal.tex
    Editable LaTeX project-proposal source; main.tex is the packaged Overleaf entry point.

0487. aero588_project/proposal latex.zip
    Packaged proposal source/graphics for Overleaf.

0488. aero588_project/proposal latex/.DS_Store
    Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.

0489. aero588_project/proposal latex/AERO588_Project_Proposal (1).pdf
    Compiled/exported project-proposal PDF, not an aerodynamic result.

0490. aero588_project/proposal latex/README.txt
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0491. aero588_project/proposal latex/graphics/mesh.png
    Selected mesh lines around the airfoil and ground. PNG raster export.

0492. aero588_project/proposal latex/graphics/mesh.svg
    Selected mesh lines around the airfoil and ground. SVG vector export.

0493. aero588_project/proposal latex/graphics/option_a_actual_mesh.png
    Proposal figure option showing selected lines of the actual CFD grid. PNG raster export.

0494. aero588_project/proposal latex/graphics/option_a_actual_mesh.svg
    Proposal figure option showing selected lines of the actual CFD grid. SVG vector export.

0495. aero588_project/proposal latex/graphics/option_b_compact.png
    Compact proposal schematic of the design box, airfoil shape changes, and moving ground. PNG raster export.

0496. aero588_project/proposal latex/graphics/option_b_compact.svg
    Compact proposal schematic of the design box, airfoil shape changes, and moving ground. SVG vector export.

0497. aero588_project/proposal latex/graphics/option_b_shape_schematic.png
    Proposal schematic of the FFD/design-box optimization concept; not a CFD solution. PNG raster export.

0498. aero588_project/proposal latex/graphics/option_b_shape_schematic.svg
    Proposal schematic of the FFD/design-box optimization concept; not a CFD solution. SVG vector export.

0499. aero588_project/proposal latex/graphics/static_pressure_cp.png
    Static-pressure coefficient contour of the saved baseline/single CFD case. PNG raster export.

0500. aero588_project/proposal latex/graphics/static_pressure_cp.svg
    Static-pressure coefficient contour of the saved baseline/single CFD case. SVG vector export.

0501. aero588_project/proposal latex/graphics/surface_cp.png
    Surface Cp diagram of the saved baseline/single CFD case. PNG raster export.

0502. aero588_project/proposal latex/graphics/surface_cp.svg
    Surface Cp diagram of the saved baseline/single CFD case. SVG vector export.

0503. aero588_project/proposal latex/graphics/total_pressure.png
    Total-pressure contour of the saved baseline/single CFD case. PNG raster export.

0504. aero588_project/proposal latex/graphics/total_pressure.svg
    Total-pressure contour of the saved baseline/single CFD case. SVG vector export.

0505. aero588_project/proposal latex/main.tex
    Editable LaTeX project-proposal source; main.tex is the packaged Overleaf entry point.

0506. aero588_project/reports/.DS_Store
    Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.

0507. aero588_project/reports/baseline_a0.json
    Legacy baseline aerodynamic coefficients, dimensional/reference conditions, runtime, and convergence results.

0508. aero588_project/reports/baseline_a0.log
    Recorded execution output for the legacy baseline solve. Inspect the associated result/status before treating a run as converged.

0509. aero588_project/reports/baseline_a0_attempt1.log
    Recorded execution output for the archived attempt named by this file. Inspect the associated result/status before treating a run as converged.

0510. aero588_project/reports/baseline_a0_attempt2.log
    Recorded execution output for the archived attempt named by this file. Inspect the associated result/status before treating a run as converged.

0511. aero588_project/reports/baseline_a0_converged_history.log
    Recorded execution output for the accepted legacy baseline residual/force history. Inspect the associated result/status before treating a run as converged.

0512. aero588_project/reports/baseline_convergence.png
    Residual and aerodynamic convergence history of the legacy baseline. PNG raster export.

0513. aero588_project/reports/baseline_pressure.png
    Surface pressure-coefficient plot of the legacy baseline. PNG raster export.

0514. aero588_project/reports/baseline_report.md
    Human-readable legacy baseline conditions, force/convergence results, and model/mesh limitations.

0515. aero588_project/reports/baseline_surface.csv
    Extracted legacy surface positions, normalized coordinates, Cp, skin-friction magnitude, and y-plus for plotting.

0516. aero588_project/reports/current_run_plots/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0517. aero588_project/reports/current_run_plots/baseline_fields.npz
    Cached NumPy arrays extracted from a baseline CGNS solution for plotting, with ghost/rind data removed.

0518. aero588_project/reports/current_run_plots/field_validation.json
    Validation information for extracted plotting fields and their source solution.

0519. aero588_project/reports/current_run_plots/mesh.png
    Selected mesh lines around the airfoil and ground. PNG raster export.

0520. aero588_project/reports/current_run_plots/mesh.svg
    Selected mesh lines around the airfoil and ground. SVG vector export.

0521. aero588_project/reports/current_run_plots/source_provenance.json
    Origin paths/hashes and reference information for the plotting cache and source solution.

0522. aero588_project/reports/current_run_plots/static_pressure_cp.png
    Static-pressure coefficient contour of the saved baseline/single CFD case. PNG raster export.

0523. aero588_project/reports/current_run_plots/static_pressure_cp.svg
    Static-pressure coefficient contour of the saved baseline/single CFD case. SVG vector export.

0524. aero588_project/reports/current_run_plots/surface_cp.png
    Surface Cp diagram of the saved baseline/single CFD case. PNG raster export.

0525. aero588_project/reports/current_run_plots/surface_cp.svg
    Surface Cp diagram of the saved baseline/single CFD case. SVG vector export.

0526. aero588_project/reports/current_run_plots/total_pressure.png
    Total-pressure contour of the saved baseline/single CFD case. PNG raster export.

0527. aero588_project/reports/current_run_plots/total_pressure.svg
    Total-pressure contour of the saved baseline/single CFD case. SVG vector export.

0528. aero588_project/reports/environment_check.json
    Results of MACH-Aero package, optimizer, and MPI smoke checks.

0529. aero588_project/reports/mesh_build.log
    Recorded execution output for mesh generation. Inspect the associated result/status before treating a run as converged.

0530. aero588_project/reports/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0531. aero588_project/reports/mesh_preview.png
    Preview of the original baseline overset mesh. PNG raster export.

0532. aero588_project/reports/naca0012_current_run_plots.zip
    Bundled CFD comparison figures for sharing.

0533. aero588_project/reports/proposal_figure_options/__pycache__/generate_figures.cpython-312.pyc
    Legacy tracked Python bytecode cache; regenerated by Python and not needed to reproduce results.

0534. aero588_project/reports/proposal_figure_options/generate_figures.py
    Generate proposal mesh/schematic figure alternatives from existing data; a presentation helper, not a CFD solver.

0535. aero588_project/reports/proposal_figure_options/option_a_actual_mesh.png
    Proposal figure option showing selected lines of the actual CFD grid. PNG raster export.

0536. aero588_project/reports/proposal_figure_options/option_a_actual_mesh.svg
    Proposal figure option showing selected lines of the actual CFD grid. SVG vector export.

0537. aero588_project/reports/proposal_figure_options/option_b_compact.png
    Compact proposal schematic of the design box, airfoil shape changes, and moving ground. PNG raster export.

0538. aero588_project/reports/proposal_figure_options/option_b_compact.svg
    Compact proposal schematic of the design box, airfoil shape changes, and moving ground. SVG vector export.

0539. aero588_project/reports/proposal_figure_options/option_b_shape_schematic.png
    Proposal schematic of the FFD/design-box optimization concept; not a CFD solution. PNG raster export.

0540. aero588_project/reports/proposal_figure_options/option_b_shape_schematic.svg
    Proposal schematic of the FFD/design-box optimization concept; not a CFD solution. SVG vector export.

0541. aero588_project/results/baseline_a0.json
    Legacy baseline aerodynamic coefficients, dimensional/reference conditions, runtime, and convergence results.

0542. aero588_project/results/baseline_surface.csv
    Extracted legacy surface positions, normalized coordinates, Cp, skin-friction magnitude, and y-plus for plotting.

0543. aero588_project/runs/baseline_a0/baseline_a0_000_slices.dat
    Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.

0544. aero588_project/runs/baseline_a0/baseline_a0_000_surf.cgns
    Saved surface mesh and CFD solution in CGNS format.

0545. aero588_project/runs/baseline_a0/baseline_a0_000_surf.plt
    Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.

0546. aero588_project/runs/baseline_a0/baseline_a0_000_vol.cgns
    Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.

0547. aero588_project/scripts/__pycache__/adf_reader.cpython-312.pyc
    Legacy tracked Python bytecode cache; regenerated by Python and not needed to reproduce results.

0548. aero588_project/scripts/adf_reader.py
    Read CGNS ADF trees through the official shared library; does not modify the database.

0549. aero588_project/scripts/airfoil_geometry.py
    Compute actual-surface maximum chord-normal thickness and its active-feature derivative; detect crossed/collapsed section boundaries.

0550. aero588_project/scripts/build_ground_mesh.py
    Build the original NACA 0012 body-fitted/background overset mesh with slip ground; legacy /project container workflow.

0551. aero588_project/scripts/check_environment.ps1
    Windows/WSL PowerShell launcher for the MACH-Aero import, MPI, and optimizer environment checks.

0552. aero588_project/scripts/check_optimization_geometry.py
    Check pyGeo derivatives of actual-surface box and thickness functions at a selected or perturbed design.

0553. aero588_project/scripts/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer.

0554. aero588_project/scripts/check_section_thickness.py
    Test analytical thickness examples, shifted upper/lower vertices, concave sections, rotation invariance, invalid topology, and derivatives.

0555. aero588_project/scripts/export_optimization_ffd.py
    Reconstruct saved pyGeo design variables, export actual FFD coefficients, and verify the reconstructed surfaces against saved geometry.

0556. aero588_project/scripts/extract_cfd_fields.py
    Extract plotting fields and surface Cp/Cf/y-plus from a completed single CFD case using the ADF reader.

0557. aero588_project/scripts/extract_optimization_fields.py
    Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations.

0558. aero588_project/scripts/finalize_optimization_run.py
    Assemble resolved settings/provenance and check saved geometry, CFD convergence, mesh validity, and source checksums.

0559. aero588_project/scripts/finish_xy_optimization.py
    Host completion worker for the October 29 run: wait for Docker, validate/extract final results, generate plots/reports, and mirror the results.

0560. aero588_project/scripts/generate_baseline.py
    Generate the original unrotated 0.70 m NACA 0012 geometry at 0.20 m clearance; does not run CFD.

0561. aero588_project/scripts/generate_main_readme.py
    Generate this file-by-file plain-text guide from the repository inventory and save identical local/GitHub copies.

0562. aero588_project/scripts/plot_current_run.py
    Plot an existing baseline/single-case mesh, surface Cp, static-pressure contour, and total-pressure contour from extracted fields.

0563. aero588_project/scripts/plot_ffd_mesh_comparison.py
    Plot saved baseline/final FFD lattices over meshes and compare actual x/y control-point motion; no CFD solve.

0564. aero588_project/scripts/plot_mesh.py
    Plot selected lines from the original combined baseline CGNS mesh using the container libraries.

0565. aero588_project/scripts/plot_optimization_comparison.py
    Generate before/after geometry, mesh, Cp, pressure/loss, load, and optimization-history figures plus numerical comparison reports.

0566. aero588_project/scripts/postprocess_baseline.py
    Read the legacy baseline CFD history and create baseline convergence/pressure summaries and plots.

0567. aero588_project/scripts/run_baseline.py
    Solve the original zero-incidence slip-ground RANS baseline; legacy /project container workflow.

0568. aero588_project/scripts/run_ground_case.py
    Build/solve a case from case_config.json; supports slip, stationary no-slip, and patched moving no-slip ground.

0569. aero588_project/scripts/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories.

0570. aero588_project/scripts/sync_optimization_progress.py
    Mirror stable files from local Docker staging into the project and write progress summaries; does not compute CFD.

0571. aero588_project/scripts/verify_environment.py
    Smoke-check MACH-Aero imports, MPI behavior, and optimizer initialization without running CFD.

0572. aero588_project/setups/slip_ground_xy5_cd005/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0573. aero588_project/setups/slip_ground_xy5_cd005/RUN_STATUS.txt
    Human-readable lifecycle/stopping status; distinguish prepared geometry, completed CFD, and optimizer termination.

0574. aero588_project/setups/slip_ground_xy5_cd005/baseline/design.json
    Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.

0575. aero588_project/setups/slip_ground_xy5_cd005/baseline/surface_geometry.csv
    Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.

0576. aero588_project/setups/slip_ground_xy5_cd005/ffd_geometry/baseline_control_points.csv
    Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.

0577. aero588_project/setups/slip_ground_xy5_cd005/ffd_geometry/baseline_ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0578. aero588_project/setups/slip_ground_xy5_cd005/ffd_geometry/control_points.npz
    Actual pyGeo baseline/final FFD coefficient arrays used by the FFD comparison plots.

0579. aero588_project/setups/slip_ground_xy5_cd005/ffd_geometry/verification.json
    FFD reconstruction errors and point counts; checks the exported lattice against the saved actual surface.

0580. aero588_project/setups/slip_ground_xy5_cd005/geometry_derivative_check_extended.json
    Additional independent geometry derivative directions and the design at which they were checked.

0581. aero588_project/setups/slip_ground_xy5_cd005/geometry_validation.json
    Geometry-only baseline metrics and derivative checks performed before the optimization was launched.

0582. aero588_project/setups/slip_ground_xy5_cd005/section_thickness_checks.json
    Analytical and derivative tests of the actual-surface chord-normal thickness calculation.

0583. aero588_project/setups/slip_ground_xy5_cd005/source/airfoil_coordinates.csv
    Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.

0584. aero588_project/setups/slip_ground_xy5_cd005/source/airfoil_geometry.py
    Compute actual-surface maximum chord-normal thickness and its active-feature derivative; detect crossed/collapsed section boundaries. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0585. aero588_project/setups/slip_ground_xy5_cd005/source/check_optimization_geometry.py
    Check pyGeo derivatives of actual-surface box and thickness functions at a selected or perturbed design. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0586. aero588_project/setups/slip_ground_xy5_cd005/source/check_optimized_connectivity.py
    Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0587. aero588_project/setups/slip_ground_xy5_cd005/source/check_section_thickness.py
    Test analytical thickness examples, shifted upper/lower vertices, concave sections, rotation invariance, invalid topology, and derivatives. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0588. aero588_project/setups/slip_ground_xy5_cd005/source/export_optimization_ffd.py
    Reconstruct saved pyGeo design variables, export actual FFD coefficients, and verify the reconstructed surfaces against saved geometry. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0589. aero588_project/setups/slip_ground_xy5_cd005/source/extract_optimization_fields.py
    Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0590. aero588_project/setups/slip_ground_xy5_cd005/source/ffd.xyz
    Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design's actual coefficients.

0591. aero588_project/setups/slip_ground_xy5_cd005/source/finalize_optimization_run.py
    Assemble resolved settings/provenance and check saved geometry, CFD convergence, mesh validity, and source checksums. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0592. aero588_project/setups/slip_ground_xy5_cd005/source/geometry_derivative_check.json
    Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation.

0593. aero588_project/setups/slip_ground_xy5_cd005/source/geometry_setup.json
    FFD dimensions, variable groups/counts/units, drag cap, box bounds, and thickness definition for the setup.

0594. aero588_project/setups/slip_ground_xy5_cd005/source/ground.cgns
    Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.

0595. aero588_project/setups/slip_ground_xy5_cd005/source/initial_design.json
    Zero-shape-change initial FFD variables and initial incidence for this parameterization.

0596. aero588_project/setups/slip_ground_xy5_cd005/source/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0597. aero588_project/setups/slip_ground_xy5_cd005/source/plot_ffd_mesh_comparison.py
    Plot saved baseline/final FFD lattices over meshes and compare actual x/y control-point motion; no CFD solve. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0598. aero588_project/setups/slip_ground_xy5_cd005/source/plot_optimization_comparison.py
    Generate before/after geometry, mesh, Cp, pressure/loss, load, and optimization-history figures plus numerical comparison reports. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0599. aero588_project/setups/slip_ground_xy5_cd005/source/restart_vol.cgns
    Saved volume flow/mesh solution used to warm-start this run or independent final-design check.

0600. aero588_project/setups/slip_ground_xy5_cd005/source/run_slip_optimization.py
    Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories. Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.

0601. aero588_project/setups/slip_ground_xy5_cd005/source/sha256_manifest.json
    SHA-256 hashes of the associated input/source files; verifies the preserved run snapshot.

0602. aero588_project/setups/slip_ground_xy5_cd005/source/software_versions.json
    Recorded solver/library versions, source commits, MPI, Python, and platform information.

0603. aero588_project/solver_extensions/moving_wall/Dockerfile
    Build the separate patched ADflow moving-wall image; does not replace the stock slip-ground image.

0604. aero588_project/solver_extensions/moving_wall/IMPLEMENTATION_REPORT.md
    Detailed moving-wall implementation, validation findings, limitations, and airfoil boundary-test results.

0605. aero588_project/solver_extensions/moving_wall/README.md
    Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.

0606. aero588_project/solver_extensions/moving_wall/adflow-moving-wall.patch
    Exact upstream Python/Fortran source differences implementing the moving-wall API and associated spatial-derivative changes.

0607. aero588_project/solver_extensions/moving_wall/inspect_wall_resolution.py
    Inspect moving-ground mesh/wall-resolution data and report airfoil/road y-plus diagnostics.

0608. aero588_project/solver_extensions/moving_wall/patch_adflow.py
    Apply the version-checked moving-wall API and primal/adjoint source patch to ADflow 2.12.2.

0609. aero588_project/solver_extensions/moving_wall/prescribedWallMotion.F90
    Fortran module storing prescribed family wall velocities and exposing dimensional diagnostic statistics.

0610. aero588_project/solver_extensions/moving_wall/validate_moving_wall.py
    Run moving/stationary floor and wall-velocity/spatial-derivative validation cases, including MPI diagnostics.

0611. aero588_project/solver_extensions/moving_wall/validation/airfoil_100um/cfd.log
    Recorded execution output for the case CFD solve. Inspect the associated result/status before treating a run as converged.

0612. aero588_project/solver_extensions/moving_wall/validation/airfoil_100um/cfd_result.json
    Case CFD results and convergence flags; use the enclosing README to distinguish accepted and explicitly unconverged runs.

0613. aero588_project/solver_extensions/moving_wall/validation/airfoil_100um/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.

0614. aero588_project/solver_extensions/moving_wall/validation/airfoil_100um/wall_resolution.json
    Airfoil/road wall-spacing and y-plus assessment; a small y-plus alone does not establish CFD convergence.

0615. aero588_project/solver_extensions/moving_wall/validation/build.log
    Recorded execution output for the patched ADflow Docker build. Inspect the associated result/status before treating a run as converged.

0616. aero588_project/solver_extensions/moving_wall/validation/build_metadata.json
    Docker build/image/source information recorded for the moving-wall extension.

0617. aero588_project/solver_extensions/moving_wall/validation/flat_moving.log
    Recorded execution output for the prescribed moving-floor validation case. Inspect the associated result/status before treating a run as converged.

0618. aero588_project/solver_extensions/moving_wall/validation/flat_moving/validation.json
    Moving-wall/control-case validation metrics, wall-velocity checks, and spatial-derivative results.

0619. aero588_project/solver_extensions/moving_wall/validation/flat_stationary.log
    Recorded execution output for the stationary-floor control case. Inspect the associated result/status before treating a run as converged.

0620. aero588_project/solver_extensions/moving_wall/validation/flat_stationary/validation.json
    Moving-wall/control-case validation metrics, wall-velocity checks, and spatial-derivative results.

0621. aero588_project/solver_extensions/moving_wall/validation/stiff_2um_attempt/airfoil_cfd.log
    Recorded execution output for the airfoil boundary-validation attempt. Inspect the associated result/status before treating a run as converged. Archived unsuccessful fine-road-grid attempt; not an accepted final aerodynamic solution.

0622. aero588_project/solver_extensions/moving_wall/validation/stiff_2um_attempt/case_config.json
    Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline. Archived unsuccessful fine-road-grid attempt; not an accepted final aerodynamic solution.

0623. aero588_project/solver_extensions/moving_wall/validation/stiff_2um_attempt/mesh_check.json
    Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation. Archived unsuccessful fine-road-grid attempt; not an accepted final aerodynamic solution.

0624. aero588_project/tmp/pdfs/proposal_slides/contact.png
    Contact sheet of rendered presentation slides used during proposal preparation. PNG raster export.

0625. aero588_project/tmp/pdfs/proposal_slides/slide-1.png
    Rendered original presentation slide 1; proposal-reference material. PNG raster export.

0626. aero588_project/tmp/pdfs/proposal_slides/slide-2.png
    Rendered original presentation slide 2; proposal-reference material. PNG raster export.

0627. aero588_project/tmp/pdfs/proposal_slides/slide-3.png
    Rendered original presentation slide 3; proposal-reference material. PNG raster export.

0628. aero588_project/tmp/pdfs/proposal_slides/slide-4.png
    Rendered original presentation slide 4; proposal-reference material. PNG raster export.

0629. aero588_project/tmp/pdfs/proposal_slides/slide-5.png
    Rendered original presentation slide 5; proposal-reference material. PNG raster export.

0630. aero588_project/tmp/pdfs/proposal_slides/slide-6.png
    Rendered original presentation slide 6; proposal-reference material. PNG raster export.

0631. aero588_project/tmp/pdfs/proposal_slides/slide-7.png
    Rendered original presentation slide 7; proposal-reference material. PNG raster export.

0632. aero588_project/tmp/pdfs/proposal_slides/slide-8.png
    Rendered original presentation slide 8; proposal-reference material. PNG raster export.

0633. aero588_project/tmp/pdfs/proposal_slides/slide-9.png
    Rendered original presentation slide 9; proposal-reference material. PNG raster export.

0634. aero588_project/tmp/pdfs/proposal_validation.png
    Rendered proposal page used to inspect layout; not aerodynamic validation. PNG raster export.

0635. main_readme.txt
    This plain-text navigation guide and complete file-purpose inventory; identical copies are provided at repository and project roots.

LOCAL-ONLY BOOKKEEPING
reports/github_upload_plan_20261029.md, when present locally, records the
upload preparation/approval and commit bookkeeping. It is not a CFD input
and is not included in the public simulation-data snapshot.
