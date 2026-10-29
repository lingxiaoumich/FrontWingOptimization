"""Finish the authorized xy5 run after Docker returns; save checks and plots.

Runs on the host, using the recorded MDO Lab image for CGNS/pyGeo/CFD work.
This worker does not restart or extend the optimizer.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

p=argparse.ArgumentParser()
p.add_argument('--stage',type=Path,required=True)
p.add_argument('--destination',type=Path,required=True)
p.add_argument('--container',required=True)
a=p.parse_args()
stage=a.stage.resolve();run=stage/'run';destination=a.destination.resolve()
image='mdolab/public@sha256:bfb64541108e7f819b7abe0f2f119ca08f1bbd7a5ea9ebe7377e90ad8f0d9333'
env=os.environ.copy();env['MPLCONFIGDIR']='/private/tmp/aero588-matplotlib'


def status(message):
    print(message,flush=True)
    for path in [run/'RUN_STATUS.txt',destination/'RUN_STATUS.txt']:
        path.write_text(message+'\n')


def docker(command):
    with (run/'postprocess.log').open('a') as log:
        subprocess.run(['docker','run','--rm','-v',str(stage)+':/project',image,
                        'bash','-lc','set -e\nsource /home/mdolabuser/.bashrc_mdolab\n'
                        'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n'+command],
                       stdout=log,stderr=subprocess.STDOUT,check=True)


def host(script,*args):
    with (run/'postprocess.log').open('a') as log:
        subprocess.run([sys.executable,str(run/'source'/script),*map(str,args)],
                       env=env,stdout=log,stderr=subprocess.STDOUT,check=True)


try:
    status('RUNNING: PSQP optimization, five FFD pairs with x/y motion, CD <= 0.05; maximum 50 major iterations.')
    waited=subprocess.run(['docker','wait',a.container],capture_output=True,text=True,check=True)
    if waited.stdout.strip()!='0' or not (run/'run_complete.json').exists():
        raise RuntimeError('Optimization process did not return successfully: '+waited.stdout.strip())
    status('OPTIMIZER RETURNED: checking saved grids, exporting FFD coordinates, and validating full overset connectivity.')
    docker('python /project/run/source/extract_optimization_fields.py --root /project/run\n'
           'python /project/run/source/export_optimization_ffd.py --root /project/run --labels baseline final best_feasible')
    result=json.loads((run/'final/result.json').read_text())
    selected='final' if result['feasible'] and not result['cfd_failed'] else 'best_feasible'
    validation=run/'full_connectivity_validation'
    (validation/'source').mkdir(parents=True,exist_ok=True)
    for name in ['ground.cgns','airfoil_coordinates.csv','mesh_check.json','software_versions.json',
                 'airfoil_geometry.py','run_slip_optimization.py','check_optimized_connectivity.py',
                 'extract_optimization_fields.py']:
        shutil.copy2(run/'source'/name,validation/'source'/name)
    shutil.copy2(run/selected/'flow_000_vol.cgns',validation/'source/restart_vol.cgns')
    shutil.copy2(run/selected/'design.json',validation/'source/design_to_check.json')
    docker('mpirun -np 4 python /project/run/full_connectivity_validation/source/check_optimized_connectivity.py '
           '--root /project/run/full_connectivity_validation --ffd-stations 5 --shape-directions xy --drag-cap 0.05 '
           '> /project/run/full_connectivity_validation/run.log 2>&1\n'
           'python /project/run/full_connectivity_validation/source/extract_optimization_fields.py '
           '--root /project/run/full_connectivity_validation')
    validated=json.loads((validation/'final/result.json').read_text())
    if not validated['accepted_convergence']:
        raise RuntimeError('Full-connectivity CFD check did not converge')
    status('POSTPROCESSING: generating comparison figures, settings report, and saved-result checks.')
    host('plot_optimization_comparison.py','--root',run,'--after',selected)
    host('plot_ffd_mesh_comparison.py','--root',run,'--after',selected)
    host('finalize_optimization_run.py','--root',run)
    completed=json.loads((run/'run_complete.json').read_text())
    baseline=json.loads((run/'baseline/result.json').read_text())
    chosen=json.loads((run/selected/'result.json').read_text())
    gain=100*(chosen['cl']/baseline['cl']-1)
    checked_gain=100*(validated['cl']/baseline['cl']-1)
    inform=completed['inform']
    optimizer_status=('stopping criterion reached' if inform['value']==0 else
                      'iteration limit reached; stationarity not reached' if inform['value']==11 else str(inform))
    summary={'optimizer_inform':inform,'selected_folder':selected,
             'downforce_increase_percent':gain,'full_connectivity_downforce_increase_percent':checked_gain,
             'independent_full_connectivity_feasible':validated['feasible'],
             'single_grid_only':True,'finished_utc':time.strftime('%Y-%m-%d %H:%M:%S UTC',time.gmtime())}
    (run/'completion_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (run/'README.md').write_text(f'''# Five-pair x/y slip-ground optimization

The PSQP run returned: **{optimizer_status}**. This is a preliminary single-grid result. The selected design is `{selected}/`, with {gain:.3f}% more downforce than the baseline. An independent CFD run with full overset connectivity gives {checked_gain:.3f}% more downforce and reports feasibility = {validated['feasible']}.

The model uses five upper/lower FFD pairs, tied across two span planes: 10 independent x displacements, 10 y displacements, and free incidence (21 variables). CD ≤ 0.05. Initial chord is 0.5 m, freestream is 200 km/h, and the ground is an ideal slip wall. The actual airfoil surface is constrained to x = [−0.125, 0.625] m and y = [0.05, 0.35] m. Maximum section thickness must be at least 0.03 m, calculated from actual-surface chord-normal intersections.

## Results and provenance

- `comparison_results.txt` and `.json`: baseline/selected force coefficients, dimensional loads, clearance, chord, thickness, incidence, solver status, and independent validation.
- `settings.txt` and `.json`: geometry, CFD, mesh, optimizer settings, units, and resolved solver defaults.
- `plots/`: PNG/SVG comparisons of geometry, mesh, Cp, pressure contours/loss, loads, iteration history, and FFD movement in x and y.
- `baseline/`, `final/`, `best_feasible/`: surface coordinates, physical design variables, CGNS flow/mesh solutions, extracted fields, and result JSON.
- `ffd_geometry/`: actual pyGeo lattice coordinates and reconstruction verification.
- `validation.json`: numerical/geometry/mesh checks. `full_connectivity_validation/` contains the independent final CFD evaluation, log, settings, and mesh checks.
- `run.log`, `PSQP.out`, `history.jsonl`, `gradients.jsonl`, `optimization.hst`, and `optimizer_exit.json`: full run and optimization history.
- `source/`: executed scripts, input mesh, initial flow/geometry, software versions, derivative checks, and source checksums.

Optimization uses frozen overset interpolation connectivity while IDWarp deforms the whole mesh. The independent final check rebuilds connectivity. Grid independence has not been established. FFD and incidence remain unbounded explicitly; PSQP's step limit is numerical, and the actual surface box/thickness are the physical geometric constraints. The thickness maximum is piecewise differentiable when its active location changes.

The original October 27 optimization and the October 29 prepared setup remain in their separate folders.
''')
    status('COMPLETED: '+optimizer_status+f'. Selected downforce gain {gain:.3f}%; independent full-connectivity gain {checked_gain:.3f}%.')
    shutil.copytree(run,destination,dirs_exist_ok=True)
    print(json.dumps(summary),flush=True)
except Exception as exc:
    status('ACTION REQUIRED: '+str(exc)+'. See run.log and postprocess.log; completed CFD data are preserved.')
    shutil.copytree(run,destination,dirs_exist_ok=True)
    raise
