"""Recreate the saved flow of an iterate whose flow file was removed by prune_iterate_flows.py.

  python scripts/regenerate_iterate_flow.py optimizations/<run>/iterates/eval_NNN [more iterates ...]

A pruned iterate keeps design.json (its design variables), result.json and the surface. This
solves that design again with the run's own script, grid and options (<run>/source and
command.txt), starting from the nearest iterate of the run that still has its flow, puts
flow_000_vol.cgns back in the iterate folder and compares the forces with result.json.

Run it in the mdolab environment with the machine free: it uses the run's rankfile (64 ranks).
--check only compares the forces and leaves the iterate folder as it is.
"""
import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('iterates', type=Path, nargs='+')
p.add_argument('--start', type=Path, help='iterate folder to start from (default: the nearest one with a flow file)')
p.add_argument('--check', action='store_true', help='compare the forces only; do not put the flow file back')
a = p.parse_args()

number = lambda d: int(d.name.split('_')[1])
for it in a.iterates:
    it = it.resolve(); run = it.parents[1]
    want = json.loads((it/'result.json').read_text())
    if (it/'flow_000_vol.cgns').exists() and not a.check:
        print(it, 'already has its flow file; nothing done'); continue
    kept = [d for d in (run/'iterates').glob('eval_*') if d != it and (d/'flow_000_vol.cgns').exists() and (d/'design.json').exists()]
    start = a.start.resolve() if a.start else min(kept, key=lambda d: abs(number(d)-number(it)))
    root = run/f'regenerate_{it.name}'
    shutil.copytree(run/'source', root/'source', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy(run/'rankfile', root/'rankfile')
    cmd = (run/'command.txt').read_text().split('\n')[0].split()
    for flag, values in [('--restart-design', 1), ('--restart-flow', 1), ('--hotstart', 1), ('--reuse-checks', 0), ('--evaluate-design', 1),
                         ('--save-iterate-flows', 0), ('--iterate-volume-only', 0)]:
        if flag in cmd: del cmd[cmd.index(flag):cmd.index(flag)+1+values]
    cmd[cmd.index('--root')+1] = str(root)
    cmd[[i for i, w in enumerate(cmd) if w.endswith('/source/run_slip_optimization.py')][0]] = str(root/'source/run_slip_optimization.py')
    # the first solve is the start iterate, from its own converged flow; then the wanted design is solved from there
    cmd += ['--restart-design', str(start/'design.json'), '--restart-flow', str(start/'flow_000_vol.cgns'), '--reuse-checks',
            '--evaluate-design', str(it/'design.json')]
    t = time.time()
    with (root/'run.log').open('w') as log:
        subprocess.run(' '.join(cmd), shell=True, cwd=root, stdout=log, stderr=subprocess.STDOUT)
    out = json.loads((root/'run_complete.json').read_text()); got = out['final']
    report = {'iterate': it.name, 'started_from': start.name, 'seconds': round(time.time()-t), 'failed': bool(out['final_failed']),
              'cl_saved': want['cl'], 'cl_regenerated': got['cl'], 'cl_relative_difference': abs(got['cl']-want['cl'])/abs(want['cl']),
              'cd_saved': want['cd'], 'cd_regenerated': got['cd'], 'cd_relative_difference': abs(got['cd']-want['cd'])/abs(want['cd'])}
    print(json.dumps(report), flush=True)
    if out['final_failed']:
        raise SystemExit(f'{it}: the solve failed; {root} is kept for a look')
    if not a.check:
        shutil.move(str(root/'final/flow_000_vol.cgns'), it/'flow_000_vol.cgns')
        (it/'regenerated.json').write_text(json.dumps(report, indent=1)+'\n')
    shutil.rmtree(root)
