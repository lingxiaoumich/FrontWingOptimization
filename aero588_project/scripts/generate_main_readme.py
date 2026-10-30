"""Generate a plain-text guide with a purpose for every version-controlled file.

Use --repo-root for a checkout containing aero588_project/. Optionally use
--local-project to save the identical guide in the separate working project.
Only documentation is generated; CFD data and archived source are not modified.
"""
import argparse
from pathlib import Path
import subprocess


SCRIPT_PURPOSES = {
    'adf_reader.py': 'Read CGNS ADF trees through the official shared library; does not modify the database.',
    'airfoil_geometry.py': 'Compute actual-surface maximum chord-normal thickness and its active-feature derivative; detect crossed/collapsed section boundaries.',
    'build_ground_mesh.py': 'Build the original NACA 0012 body-fitted/background overset mesh with slip ground; legacy /project container workflow.',
    'check_environment.ps1': 'Windows/WSL PowerShell launcher for the MACH-Aero import, MPI, and optimizer environment checks.',
    'check_optimization_geometry.py': 'Check pyGeo derivatives of actual-surface box and thickness functions at a selected or perturbed design.',
    'check_optimized_connectivity.py': 'Independently re-evaluate a saved design with full overset connectivity updates; does not run an optimizer.',
    'check_section_thickness.py': 'Test analytical thickness examples, shifted upper/lower vertices, concave sections, rotation invariance, invalid topology, and derivatives.',
    'check_export_serialization.py': 'Archived regression check for serializing pyOptSparse termination data, NumPy arrays, and scalars to JSON.',
    'export_optimization_ffd.py': 'Reconstruct saved pyGeo design variables, export actual FFD coefficients, and verify the reconstructed surfaces against saved geometry.',
    'extract_cfd_fields.py': 'Extract plotting fields and surface Cp/Cf/y-plus from a completed single CFD case using the ADF reader.',
    'extract_optimization_fields.py': 'Extract CGNS fields for saved optimization designs, compute pressure-derived quantities, and check mesh corner orientations.',
    'finalize_optimization_run.py': 'Assemble resolved settings/provenance and check saved geometry, CFD convergence, mesh validity, and source checksums.',
    'finish_xy_optimization.py': 'Host completion worker for the October 29 run: wait for Docker, validate/extract final results, generate plots/reports, and mirror the results.',
    'generate_baseline.py': 'Generate the original unrotated 0.70 m NACA 0012 geometry at 0.20 m clearance; does not run CFD.',
    'generate_figures.py': 'Generate proposal mesh/schematic figure alternatives from existing data; a presentation helper, not a CFD solver.',
    'generate_main_readme.py': 'Generate this file-by-file plain-text guide from the repository inventory and save identical local/GitHub copies.',
    'plot_current_run.py': 'Plot an existing baseline/single-case mesh, surface Cp, static-pressure contour, and total-pressure contour from extracted fields.',
    'plot_ffd_mesh_comparison.py': 'Plot saved baseline/final FFD lattices over meshes and compare actual x/y control-point motion; no CFD solve.',
    'plot_mesh.py': 'Plot selected lines from the original combined baseline CGNS mesh using the container libraries.',
    'plot_optimization_comparison.py': 'Generate before/after geometry, mesh, Cp, pressure/loss, load, and optimization-history figures plus numerical comparison reports.',
    'postprocess_baseline.py': 'Read the legacy baseline CFD history and create baseline convergence/pressure summaries and plots.',
    'run_baseline.py': 'Solve the original zero-incidence slip-ground RANS baseline; legacy /project container workflow.',
    'run_ground_case.py': 'Build/solve a case from case_config.json; supports slip, stationary no-slip, and patched moving no-slip ground.',
    'run_case.py': 'Case-specific build/solve driver saved with this case; use that case configuration and its recorded solver environment.',
    'run_slip_optimization.py': 'Drive adjoint-based slip-ground downforce optimization, actual-surface constraints, saved designs, and CFD/gradient histories.',
    'sync_optimization_progress.py': 'Mirror stable files from local Docker staging into the project and write progress summaries; does not compute CFD.',
    'verify_environment.py': 'Smoke-check MACH-Aero imports, MPI behavior, and optimizer initialization without running CFD.',
    'inspect_wall_resolution.py': 'Inspect moving-ground mesh/wall-resolution data and report airfoil/road y-plus diagnostics.',
    'patch_adflow.py': 'Apply the version-checked moving-wall API and primal/adjoint source patch to ADflow 2.12.2.',
    'validate_moving_wall.py': 'Run moving/stationary floor and wall-velocity/spatial-derivative validation cases, including MPI diagnostics.',
    'repair_optimizer_export.py': 'Archived recovery utility: verify the final history design and restore the result export after a JSON serialization failure.',
    'executed_driver_before_export_fix.py': 'Original executed optimization driver preserved before the JSON-export repair; provenance, not the current driver.',
}

DATA_PURPOSES = {
    '.gitignore': 'Ignore new operating-system caches, Python caches, and local credential files; CFD results/logs remain version controlled.',
    'main_readme.txt': 'This plain-text navigation guide and complete file-purpose inventory; identical copies are provided at repository and project roots.',
    'case_config.json': 'Operating conditions, dimensions, mesh controls, and ground model for the enclosing case; the top-level version is the legacy 0.70 m baseline.',
    'mesh_case_config.json': 'Mesh-generation configuration preserved with the run input snapshot.',
    'aerodynamic_derivative_check.json': 'Finite-difference versus adjoint aerodynamic derivative comparison and primal-failure flags.',
    'airfoil_coordinates.csv': 'Initial physical airfoil section coordinates used for FFD embedding and mesh/geometry reconstruction.',
    'airfoil_surface.xyz': 'Extruded physical airfoil surface coordinates for mesh generation; not an FFD control-point lattice.',
    'naca0012_baseline.dat': 'Legacy NACA 0012 profile coordinates produced by generate_baseline.py.',
    'baseline_a0.json': 'Legacy baseline aerodynamic coefficients, dimensional/reference conditions, runtime, and convergence results.',
    'baseline_surface.csv': 'Extracted legacy surface positions, normalized coordinates, Cp, skin-friction magnitude, and y-plus for plotting.',
    'baseline_fields.npz': 'Cached NumPy arrays extracted from a baseline CGNS solution for plotting, with ghost/rind data removed.',
    'baseline_report.md': 'Human-readable legacy baseline conditions, force/convergence results, and model/mesh limitations.',
    'cfd_result.json': 'Case CFD results and convergence flags; use the enclosing README to distinguish accepted and explicitly unconverged runs.',
    'comparison_results.json': 'Machine-readable before/after forces, geometry, and selected independent-validation results.',
    'comparison_results.txt': 'Human-readable before/after table, reference definitions, optimizer status, constraints, and interpretation limits.',
    'completion_summary.json': 'Final completion-pipeline summary: selected result, optimizer exit, gains, independent validation, and finish time.',
    'connectivity_comparison.txt': 'Compare the saved frozen-connectivity prediction with an independent rebuilt-connectivity CFD evaluation.',
    'control_points.npz': 'Actual pyGeo baseline/final FFD coefficient arrays used by the FFD comparison plots.',
    'design.json': 'Physical design variables for this saved geometry: displacements in metres and incidence in degrees; use this instead of scaled solution.txt variables.',
    'design_to_check.json': 'Physical design variables selected as input to a geometry-derivative or independent-connectivity check.',
    'environment.json': 'Recorded solver/library environment for the enclosing CFD case.',
    'environment_check.json': 'Results of MACH-Aero package, optimizer, and MPI smoke checks.',
    'execution.json': 'Exact Docker image digest, launch command, thread settings, and staging/result locations used for the run.',
    'export_repair.txt': 'Explain the post-optimization JSON-export repair and how the recovered final design was verified.',
    'export_serialization_check.json': 'Round-trip test result for the repaired optimizer-termination JSON serialization.',
    'field_validation.json': 'Validation information for extracted plotting fields and their source solution.',
    'fields.npz': 'Extracted near-body/background mesh coordinates, flow fields, blanking, and pressure reference values for a saved design.',
    'final_checks.json': 'Final numerical/data audit, including saved-grid validity, convergence, source hashes, orphan count, and derivative-check exceptions.',
    'geometry_derivative_check.json': 'Automatic box/thickness sensitivity check errors at a small nonzero FFD perturbation.',
    'geometry_derivative_check_extended.json': 'Additional independent geometry derivative directions and the design at which they were checked.',
    'geometry_setup.json': 'FFD dimensions, variable groups/counts/units, drag cap, box bounds, and thickness definition for the setup.',
    'geometry_validation.json': 'Geometry-only baseline metrics and derivative checks performed before the optimization was launched.',
    'gradients.jsonl': 'One JSON record per recorded adjoint sensitivity evaluation, with design-group gradients and failure flags.',
    'history.jsonl': 'One JSON record per CFD evaluation, including baseline/check perturbations, optimizer trial designs, loads, feasibility, residuals, and timings.',
    'history_summary.csv': 'Compact CSV summary of recorded CFD evaluation history; evaluations are not the same as optimizer major iterations.',
    'incoming_design.json': 'Most recently submitted physical design variables; a debugging/checkpoint aid, not necessarily the selected final result.',
    'initial_design.json': 'Zero-shape-change initial FFD variables and initial incidence for this parameterization.',
    'latest.json': 'Most recently completed CFD evaluation record; use final/result.json and run_complete.json for final conclusions.',
    'mesh_check.json': 'Mesh/geometry dimensions, first spacings, cell counts, and validity information recorded during mesh generation.',
    'optimizer_exit.json': 'Optimizer termination code/message, returned physical variables, timings, and evaluation count.',
    'optimizer_iterations.csv': 'PSQP major-iteration log data: objective, evaluation counts, constraint violation, and printed gradient norm.',
    'plotting_versions.json': 'Versions of host plotting libraries used for the saved figures.',
    'pre_restart_result.json': 'Preserved provisional result from before the fine-road-grid restart; not an accepted final aerodynamic result.',
    'problem.txt': 'Printed pyOptSparse objective, variables, scaling, and constraint definition for the run.',
    'production_checks.json': 'Audit counts for completed primal/adjoint solves and the feasibility/failure flags of the production optimization.',
    'progress.json': 'Mirrored optimization progress snapshot; retained for provenance rather than used as the authoritative final result.',
    'replay_settings.json': 'Original drag cap retained for the older continuation/recovery workflow; the current driver uses its explicit --drag-cap argument.',
    'result.json': 'Aerodynamic loads, actual-geometry metrics, convergence/feasibility flags, and any extracted grid/field diagnostics for this saved design.',
    'run_complete.json': 'Driver completion record with returned optimizer status, final evaluation, best-feasible design, and wall times; file existence does not imply stationarity.',
    'section_thickness_checks.json': 'Analytical and derivative tests of the actual-surface chord-normal thickness calculation.',
    'seed_design.json': 'Physical design used to initialize an earlier feasibility/recovery check.',
    'settings.json': 'Machine-readable CFD, mesh, geometry, reference, and optimizer settings for this run/check.',
    'settings.txt': 'Human-readable CFD/optimizer settings, software provenance, and resolved solver options where recorded.',
    'setup_notes.txt': 'Notes explaining the numerical setup or an archived attempted configuration.',
    'sha256_manifest.json': 'SHA-256 hashes of the associated input/source files; verifies the preserved run snapshot.',
    'software_versions.json': 'Recorded solver/library versions, source commits, MPI, Python, and platform information.',
    'solution.txt': 'Printed optimizer solution and stopping status; variables may be optimizer-scaled, so use design.json for physical values.',
    'source_provenance.json': 'Origin paths/hashes and reference information for the plotting cache and source solution.',
    'start_design.json': 'Physical starting variables for an archived optimization continuation.',
    'surface.csv': 'Case surface positions, normalized coordinates, Cp, skin-friction magnitude, and y-plus extracted from the CFD solution.',
    'surface_cp.csv': 'Saved surface positions, reference-chord coordinates, Cp, Cf, and y-plus used for before/after plots.',
    'surface_geometry.csv': 'Actual deformed CFD section boundary coordinates, including the trailing-edge closure; this is geometry, not the FFD lattice.',
    'validation_complete.json': 'Completion/results record for the independent full-connectivity CFD evaluation.',
    'verification.json': 'FFD reconstruction errors and point counts; checks the exported lattice against the saved actual surface.',
    'wall_resolution.json': 'Airfoil/road wall-spacing and y-plus assessment; a small y-plus alone does not establish CFD convergence.',
    'Dockerfile': 'Build the separate patched ADflow moving-wall image; does not replace the stock slip-ground image.',
    'adflow-moving-wall.patch': 'Exact upstream Python/Fortran source differences implementing the moving-wall API and associated spatial-derivative changes.',
    'prescribedWallMotion.F90': 'Fortran module storing prescribed family wall velocities and exposing dimensional diagnostic statistics.',
    'build_metadata.json': 'Docker build/image/source information recorded for the moving-wall extension.',
    'IMPLEMENTATION_REPORT.md': 'Detailed moving-wall implementation, validation findings, limitations, and airfoil boundary-test results.',
}

PLOTS = {
    'mesh': 'Selected mesh lines around the airfoil and ground.',
    'mesh_preview': 'Preview of the original baseline overset mesh.',
    'mesh_comparison': 'Baseline and selected optimized meshes, plotted side by side.',
    'ffd_mesh_comparison': 'Actual FFD control points/lattice over the baseline and final computational meshes.',
    'ffd_before_after_comparison': 'Before/after airfoil and FFD overlays with control-point displacement plots.',
    'geometry_comparison': 'Baseline and selected airfoil outlines relative to the allowed surface box and ground.',
    'shape_change': 'Close-up airfoil outline and changes in surface height.',
    'cp_comparison': 'Before/after surface pressure-coefficient curves in chord and global coordinates.',
    'cp_contour_comparison': 'Before/after static-pressure coefficient contours from saved CFD fields.',
    'total_pressure_comparison': 'Before/after total pressure normalized by the freestream stagnation pressure.',
    'pressure_loss_comparison': 'Before/after total-pressure-loss coefficient; numerical pressure overshoots are retained.',
    'performance_comparison': 'Before/after downforce, drag, and downforce/drag bars.',
    'optimization_history': 'Loads, incidence, and CFD residual history over evaluations, including trials/checks.',
    'optimizer_convergence': 'PSQP objective, printed gradient norm, and constraint violation versus major iteration.',
    'static_pressure_cp': 'Static-pressure coefficient contour of the saved baseline/single CFD case.',
    'surface_cp': 'Surface Cp diagram of the saved baseline/single CFD case.',
    'total_pressure': 'Total-pressure contour of the saved baseline/single CFD case.',
    'baseline_convergence': 'Residual and aerodynamic convergence history of the legacy baseline.',
    'baseline_pressure': 'Surface pressure-coefficient plot of the legacy baseline.',
    'option_a_actual_mesh': 'Proposal figure option showing selected lines of the actual CFD grid.',
    'option_b_compact': 'Compact proposal schematic of the design box, airfoil shape changes, and moving ground.',
    'option_b_shape_schematic': 'Proposal schematic of the FFD/design-box optimization concept; not a CFD solution.',
    'contact': 'Contact sheet of rendered presentation slides used during proposal preparation.',
    'proposal_validation': 'Rendered proposal page used to inspect layout; not aerodynamic validation.',
}


def describe(path, repository):
    name=path.name
    if name in SCRIPT_PURPOSES:
        result=SCRIPT_PURPOSES[name]
        if 'source' in path.parts:
            result+=' Run/setup-specific archived copy; preserve it for reproduction, since it can differ from scripts/.'
        return result
    if name in DATA_PURPOSES:return DATA_PURPOSES[name]
    if name in {'README.md','README.txt'}:
        if path.as_posix()=='README.md':return 'GitHub landing page: current model, two run summaries, status limitations, and navigation links.'
        return 'Overview, interpretation, reproduction/use notes, and limitations for the enclosing project, case, run, or package.'
    if name=='RUN_STATUS.txt':return 'Human-readable lifecycle/stopping status; distinguish prepared geometry, completed CFD, and optimizer termination.'
    if name=='validation.json':
        if 'solver_extensions' in path.parts:return 'Moving-wall/control-case validation metrics, wall-velocity checks, and spatial-derivative results.'
        return 'Saved solution/geometry/grid/convergence checks and recorded derivative-check information.'
    if name.endswith('_control_points.csv'):return 'Actual exported FFD point indices and x/y/z coordinates for the design named in the filename.'
    if name=='ffd.xyz' or name.endswith('_ffd.xyz'):return 'Plot3D FFD lattice coordinates; source/ffd.xyz is the reference lattice, while named exports contain that design\'s actual coefficients.'
    if name.endswith('.cgns'):
        if name=='airfoil_near.cgns':return 'Body-fitted near-airfoil mesh block generated by pyHyp.'
        if name in {'ground.cgns','ground_baseline.cgns'}:return 'Combined near-body/background overset input mesh, including physical and symmetry boundary definitions.'
        if name=='restart_vol.cgns':return 'Saved volume flow/mesh solution used to warm-start this run or independent final-design check.'
        if '_surf.cgns' in name:return 'Saved surface mesh and CFD solution in CGNS format.'
        return 'Saved volume mesh and flow fields in CGNS format; usable as a matching-case restart and as plotting/verification data.'
    if name.endswith('_slices.dat'):return 'Tecplot-style section slice containing surface positions, Cp, Cf, and y-plus for the extruded section.'
    if name.endswith('.plt'):return 'Tecplot-format surface solution written by ADflow; an optional visualization export alongside the CGNS data.'
    if name=='optimization.hst':return 'Binary pyOptSparse history/checkpoint with objective, constraint, variable, and sensitivity records; reuse only with a compatible optimization problem.'
    if name in {'PSQP.out','SLSQP.out'}:return 'Native optimizer iteration/stopping log; PSQP.out can also contain redirected Fortran CFD/adjoint solver output.'
    if name.endswith('.log'):
        if name.startswith('aer588_project_proposal'):return 'LaTeX compiler diagnostics for the proposal PDF; not a CFD log.'
        return 'Recorded execution output for '+{
            'cfd.log':'the case CFD solve', 'airfoil_cfd.log':'the airfoil boundary-validation attempt',
            'run.log':'the enclosing CFD/optimization driver', 'postprocess.log':'the extraction, independent CFD check, plotting, and completion pipeline',
            'mesh_build.log':'mesh generation', 'plotting_summary.log':'figure generation',
            'validation_summary.log':'the saved final-result checks', 'build.log':'the patched ADflow Docker build',
            'flat_moving.log':'the prescribed moving-floor validation case', 'flat_stationary.log':'the stationary-floor control case',
            'startup_cfd.log':'the fine-road-grid startup attempt', 'incomplete_restart.log':'the incomplete fine-road-grid restart',
            'baseline_a0.log':'the legacy baseline solve', 'baseline_a0_converged_history.log':'the accepted legacy baseline residual/force history',
        }.get(name,'the archived attempt named by this file')+'. Inspect the associated result/status before treating a run as converged.'
    if path.suffix in {'.png','.svg'}:
        if path.stem.startswith('slide-'):purpose='Rendered original presentation slide '+path.stem.split('-')[-1]+'; proposal-reference material.'
        else:purpose=PLOTS.get(path.stem)
        if purpose:return purpose+(' PNG raster export.' if path.suffix=='.png' else ' SVG vector export.')
    if name in {'main.tex','aer588_project_proposal.tex'}:return 'Editable LaTeX project-proposal source; main.tex is the packaged Overleaf entry point.'
    if path.suffix=='.pdf':return 'Compiled/exported project-proposal PDF, not an aerodynamic result.'
    if path.suffix=='.zip':return 'Packaged proposal source/graphics for Overleaf.' if 'proposal' in name else 'Bundled CFD comparison figures for sharing.'
    if name=='.DS_Store':return 'Legacy tracked macOS Finder display metadata; has no role in geometry, CFD, or optimization.'
    if path.suffix=='.pyc':return 'Legacy tracked Python bytecode cache; regenerated by Python and not needed to reproduce results.'
    raise ValueError('Missing a specific file description: '+str(path))


def context(path):
    s=path.as_posix()
    if '/restarts/' in s:return 'Archived setup/search attempt; not the selected final optimization result.'
    if '/feasible_recovery/' in s:return 'Earlier recovery design/check retained for provenance; not the final October 27 production design.'
    if '/refinement_check_15um/' in s:return 'Unconverged moving-ground refinement attempt; its provisional forces are not accepted final results.'
    if '/stiff_2um_attempt/' in s:return 'Archived unsuccessful fine-road-grid attempt; not an accepted final aerodynamic solution.'
    if '/full_connectivity_validation/' in s:return 'Independent final-design CFD check with rebuilt connectivity on the same mesh; not grid independence.'
    return ''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--local-project',type=Path)
    args=parser.parse_args()
    repo=args.repo_root.resolve()
    output=subprocess.check_output(['git','-C',str(repo),'ls-files','-z'])
    files={Path(p.decode()) for p in output.split(b'\0') if p}
    files.update([Path('main_readme.txt'),Path('aero588_project/main_readme.txt'),
                  Path('aero588_project/scripts/generate_main_readme.py')])
    files=sorted(files,key=lambda p:p.as_posix())
    header='''FRONT WING OPTIMIZATION - MAIN FILE GUIDE
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
'''
    lines=[header,f'Total version-controlled file entries: {len(files)}\n']
    for number,path in enumerate(files,1):
        description=describe(path,repo)
        note=context(path)
        lines.append(f'{number:04d}. {path.as_posix()}\n    {description}'+(' '+note if note else '')+'\n')
    lines.append('LOCAL-ONLY BOOKKEEPING\n'
                 'reports/github_upload_plan_20261029.md, when present locally, records the\n'
                 'upload preparation/approval and commit bookkeeping. It is not a CFD input\n'
                 'and is not included in the public simulation-data snapshot.\n')
    text='\n'.join(lines)
    destinations=[repo/'main_readme.txt',repo/'aero588_project/main_readme.txt']
    if args.local_project:destinations.append(args.local_project.resolve()/'main_readme.txt')
    for path in destinations:path.write_text(text)
    print(f'Wrote {len(files)} file descriptions to '+', '.join(map(str,destinations)))


if __name__=='__main__':main()
